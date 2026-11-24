#!/usr/bin/env python3
"""Outcome-blind recipient-level P2 scorer, P-all decider, graded posteriors, and nulls.

Pure Python standard library on all acyclic pathway rules; optional lazy scipy/HiGHS
MILP fallback when cyclic metabolite conversion rules are supplied (runs under `python3 -I`).
Implements protocols/P2-plexitrophy-recipient.md with cross-family fixes:
  - Trap 9 guard: prokaryotic ('Bacteria'/'Archaea') + single-consensus ('dastool-*')
    + 95%-ANI dereplicated representative filtering.
  - per_mag/*.tsv parser: handles eggNOG `ko:K...` comma/semicolon lists (from
    merge_annotations.py) and derives contig first/last CDS strictly from complete
    contig feature order BEFORE pathway KO filtering (retaining 'unknown' edge status
    when complete order is unavailable; Astra 248964604793110528).
  - Calibrated gap guard (Astra 248968746164994048 #4): `max_false_absence_prob=0.05`
    enforced by default; uncalibrated or poor-detectability MAGs never produce clean
    recipient absences even at zero KO hits.
  - Exact P-all single-node well-founded decider (Opus 248974707227279360 & Astra
    248973173642936320 / 248974912022560768):
      * One acyclic justification per `(organism, metabolite)` and producer-only export
        (no relaying, blocking F5 `p<->q`, F7 `p1<->p2`, and F8 `D:p->D:q->B:r->D:p`).
      * Checks `is_metabolite_dep_acyclic` first: when acyclic (true for 100% of P2
        pathway rules), uses Theorem 3 + Lemma L5 (`O(|V|)` two-node `S_max` calls, zero
        bundle or support enumeration); when cyclic (coNP-complete by Theorem 3', F6),
        dispatches to an exact proof-depth MILP with single-node DAG self-certification.
  - Contig-block downsampling (`downsample_contig_blocks`, resolving #3): splits only
    within original contigs (never merging across original contig boundaries), retains
    blocks uniformly across genome positions, and records actual achieved completeness.
  - Graded biological-state posteriors (`compute_pathway_posterior` and
    `compute_fanin_state_posterior`, per grid 248970390835478528 & Astra 248970428366110720 /
    248973173642936320): continuous `P(no single partner suffices | data, E)` over latent
    biological states (supporting both conditional marginals and coupled within-donor joint
    subset PMFs `donor_joint_pmfs`), with hard calls as thresholded summaries.
"""

from __future__ import annotations

import csv
import itertools
import math
import random
import re
from dataclasses import dataclass
from typing import Dict, FrozenSet, Iterable, List, Optional, Sequence, Set, Tuple

from p2_pathway_rules import (
    C0_CONTROL_MODULES,
    PATHWAY_RULES,
    STATE_ABSENT,
    STATE_PRESENT,
    STATE_UNRESOLVED,
    GeneHit,
    GenomeAnnotations,
    PathwayRule,
    classify_coencoder_bypass,
    estimate_false_absence_prob,
    evaluate_b12_de_novo,
    evaluate_recipient_b12_auxotrophy,
    p_intact,
)

RAW_BINNER_PREFIXES = (
    "lorbin_",
    "metabat_",
    "semibin_",
    "maxbin_",
    "comebin_",
    "vamb_",
)

CALL_AUTONOMOUS = "autonomous"
CALL_NOT_ELIGIBLE = "not_eligible"
CALL_D0_AMBIENT = "d0_ambient"
CALL_D1_OR_CHAIN = "d1_or_chain"
CALL_PLEXITROPHIC = "plexitrophic"
CALL_UNRESOLVED_AMBIG = "unresolved_ambig_donor"
CALL_UNRESOLVED_NO_DONOR = "unresolved_no_donor"

ALL_P2_ENDPOINTS: Tuple[str, ...] = (
    "B12",
    "THI",
    "RED_S",
    "FIXED_N",
    *C0_CONTROL_MODULES,
)


# ---------------------------------------------------------------------------
# 1. Trap 9 Deduplication & Prokaryotic Scope Guards
# ---------------------------------------------------------------------------


def parse_domain_from_gtdb(lineage: str) -> str:
    """Return 'Bacteria', 'Archaea', 'Eukaryota', or 'U' from a GTDB lineage string."""
    if not lineage:
        return "U"
    first = lineage.split(";", 1)[0].strip()
    if first in ("d__Bacteria", "Bacteria"):
        return "Bacteria"
    if first in ("d__Archaea", "Archaea"):
        return "Archaea"
    if first.startswith("d__Euk") or first == "Eukaryota":
        return "Eukaryota"
    return "U"


def filter_derep_prokaryotic_mags(
    quality_rows: Sequence[Dict[str, str]],
    derep_representatives: Optional[Set[str]] = None,
    taxonomy_by_mag: Optional[Dict[str, str]] = None,
    consensus_prefix: str = "dastool-",
    min_completeness: float = 0.0,
    max_contamination: float = 100.0,
) -> Tuple[List[str], Dict[str, int]]:
    """Two-step Trap 9 deduplication + prokaryotic domain filter."""
    counts = {
        "raw_quality_rows": len(quality_rows),
        "consensus_bins": 0,
        "derep_bins": 0,
        "prokaryotic_bins": 0,
        "quality_passing_bins": 0,
    }
    kept: List[str] = []
    seen: Set[str] = set()
    for row in quality_rows:
        name = (row.get("Name") or row.get("bin_id") or row.get("mag_id") or "").strip()
        if not name or name in seen:
            continue
        seen.add(name)
        if not name.startswith(consensus_prefix):
            continue
        counts["consensus_bins"] += 1
        if derep_representatives is not None and name not in derep_representatives:
            continue
        counts["derep_bins"] += 1
        if taxonomy_by_mag is not None:
            dom = parse_domain_from_gtdb(taxonomy_by_mag.get(name, ""))
            if dom not in ("Bacteria", "Archaea"):
                continue
        counts["prokaryotic_bins"] += 1
        comp = float(row.get("Completeness", 0.0))
        cont = float(row.get("Contamination", 100.0))
        if comp >= min_completeness and cont <= max_contamination:
            counts["quality_passing_bins"] += 1
            kept.append(name)
    return sorted(kept), counts


def assert_non_redundant_mag_ids(mag_ids: Iterable[str], consensus_prefix: str = "dastool-") -> None:
    """Fail loudly if raw multi-binner IDs are passed alongside consensus bins (Trap 9)."""
    ids = list(mag_ids)
    has_consensus = any(m.startswith(consensus_prefix) for m in ids)
    raw_hits = [m for m in ids if m.startswith(RAW_BINNER_PREFIXES)]
    if has_consensus and raw_hits:
        raise ValueError(
            f"Trap 9 violation: raw component binner IDs ({raw_hits[:3]}) mixed with "
            f"{consensus_prefix}* consensus bins."
        )


# ---------------------------------------------------------------------------
# 2. GFF & per_mag/*.tsv Parser (Order-Invariant Contig-End & KO Extraction)
# ---------------------------------------------------------------------------

EC_FULL_RE = re.compile(r"^\d+\.\d+\.\d+\.\d+$")


def parse_ko_field(raw_ko: str) -> FrozenSet[str]:
    """Extract clean Kxxxxx IDs from KofamScan or eggNOG fallback ('ko:K...,ko:K...')."""
    if not raw_ko or raw_ko == "-":
        return frozenset()
    out: Set[str] = set()
    for tok in re.split(r"[,;\s]+", raw_ko.strip()):
        k = tok.strip().replace("ko:", "")
        if len(k) == 6 and k.startswith("K") and k[1:].isdigit():
            out.add(k)
    return frozenset(out)


def parse_ec_field(raw_ec: str) -> FrozenSet[str]:
    if not raw_ec or raw_ec == "-":
        return frozenset()
    out: Set[str] = set()
    for tok in re.split(r"[,;\s]+", raw_ec.strip()):
        e = tok.strip().replace("EC:", "")
        if EC_FULL_RE.match(e):
            out.add(e)
    return frozenset(out)


def parse_gff_cds_order(gff_lines: Iterable[str]) -> Tuple[Dict[str, Tuple[str, ...]], Dict[str, int]]:
    """Derive complete per-contig CDS feature order and contig lengths from GFF3."""
    by_contig: Dict[str, List[Tuple[int, int, str]]] = {}
    contig_lengths: Dict[str, int] = {}
    for line in gff_lines:
        if line.startswith("##sequence-region"):
            parts = line.strip().split()
            if len(parts) >= 4 and parts[3].isdigit():
                contig_lengths[parts[1]] = int(parts[3])
            continue
        if not line or line.startswith("#"):
            continue
        parts = line.rstrip("\n").split("\t")
        if len(parts) < 9 or parts[2] != "CDS":
            continue
        contig = parts[0]
        start, end = int(parts[3]), int(parts[4])
        attrs = parts[8]
        pid = None
        m_id = re.search(r"(?:^|;)ID=([^;]+)", attrs)
        m_lt = re.search(r"(?:^|;)locus_tag=([^;]+)", attrs)
        if m_id:
            pid = m_id.group(1)
        elif m_lt:
            pid = m_lt.group(1)
        if pid:
            by_contig.setdefault(contig, []).append((start, end, pid))
    order: Dict[str, Tuple[str, ...]] = {}
    for contig, items in by_contig.items():
        items_sorted = sorted(items, key=lambda x: (x[0], x[1], x[2]))
        order[contig] = tuple(x[2] for x in items_sorted)
    return order, contig_lengths


def _extract_contig_seq_index(protein_id: str, contig_id: str) -> Optional[int]:
    """Parse `<contig_id>_<positive_int>` CDS index if strictly matching contig_id."""
    if not protein_id or not contig_id:
        return None
    prefix = contig_id + "_"
    if protein_id.startswith(prefix):
        tail = protein_id[len(prefix) :]
        if tail.isdigit():
            val = int(tail)
            if val >= 1:
                return val
    return None


def parse_per_mag_rows(
    rows: Sequence[Dict[str, str]],
    mag_id: str,
    domain: str = "Bacteria",
    completeness: float = 100.0,
    contamination: float = 0.0,
    genome_size_bp: int = 0,
    contig_lengths: Sequence[int] = (),
    gff_cds_order: Optional[Dict[str, Tuple[str, ...]]] = None,
    complete_cds_table: bool = False,
) -> Tuple[GenomeAnnotations, List[GeneHit]]:
    """Parse per_mag TSV rows into GenomeAnnotations + GeneHit list (row-order invariant)."""
    by_contig_raw: Dict[str, List[Tuple[str, Optional[int], Dict[str, str]]]] = {}
    for row in rows:
        pid = (row.get("protein_id") or "").strip()
        cid = (row.get("contig_id") or "").strip()
        if not pid:
            continue
        if not cid:
            parts = pid.rsplit("_", 1)
            cid = parts[0] if len(parts) == 2 and parts[1].isdigit() else pid
        idx = _extract_contig_seq_index(pid, cid)
        by_contig_raw.setdefault(cid, []).append((pid, idx, row))

    edge_by_pid: Dict[Tuple[str, str], str] = {}

    for cid, items in by_contig_raw.items():
        if gff_cds_order is not None:
            order = gff_cds_order.get(cid)
            if not order:
                for pid, _, _ in items:
                    edge_by_pid[(cid, pid)] = "unknown"
                continue
            pos_map = {p: i for i, p in enumerate(order)}
            n_tot = len(order)
            for pid, _, _ in items:
                pos = pos_map.get(pid)
                if pos is None:
                    edge_by_pid[(cid, pid)] = "unknown"
                elif n_tot >= 3 and 0 < pos < n_tot - 1:
                    edge_by_pid[(cid, pid)] = "internal"
                else:
                    edge_by_pid[(cid, pid)] = "end"
        else:
            indices = [idx for _, idx, _ in items]
            if any(i is None for i in indices):
                for pid, _, _ in items:
                    edge_by_pid[(cid, pid)] = "unknown"
                continue
            int_indices = [int(i) for i in indices if i is not None]
            if len(set(int_indices)) != len(int_indices):
                for pid, _, _ in items:
                    edge_by_pid[(cid, pid)] = "unknown"
                continue
            lo, hi = min(int_indices), max(int_indices)
            is_complete_consecutive = (
                complete_cds_table
                and lo == 1
                and set(int_indices) == set(range(1, hi + 1))
            )
            for pid, idx, _ in items:
                assert idx is not None
                if lo < idx < hi:
                    edge_by_pid[(cid, pid)] = "internal"
                elif idx == 1:
                    edge_by_pid[(cid, pid)] = "end"
                elif is_complete_consecutive and idx == hi:
                    edge_by_pid[(cid, pid)] = "end"
                else:
                    edge_by_pid[(cid, pid)] = "unknown"

    hits: List[GeneHit] = []
    all_kos: Set[str] = set()
    internal_kos: Set[str] = set()
    end_seen_kos: Set[str] = set()
    unknown_seen_kos: Set[str] = set()

    for cid, items in sorted(by_contig_raw.items()):
        # Sort hits canonically by parsed CDS index (if available) then protein_id
        items_canonical = sorted(
            items, key=lambda x: (0 if x[1] is not None else 1, x[1] or 0, x[0])
        )
        for pid, _, row in items_canonical:
            kos = parse_ko_field(row.get("KO", ""))
            ecs = parse_ec_field(row.get("EC", ""))
            estatus = edge_by_pid[(cid, pid)]
            hits.append(
                GeneHit(
                    protein_id=pid,
                    contig_id=cid,
                    kos=kos,
                    ecs=ecs,
                    edge_status=estatus,
                )
            )
            all_kos |= kos
            if estatus == "internal":
                internal_kos |= kos
            elif estatus == "end":
                end_seen_kos |= kos
            else:
                unknown_seen_kos |= kos

    end_only_kos = end_seen_kos - internal_kos
    unknown_only_kos = unknown_seen_kos - internal_kos - end_only_kos

    genome = GenomeAnnotations(
        mag_id=mag_id,
        domain=domain,
        completeness=completeness,
        contamination=contamination,
        genome_size_bp=genome_size_bp or sum(contig_lengths),
        contig_lengths=tuple(contig_lengths),
        kos=frozenset(all_kos),
        end_kos=frozenset(end_only_kos),
        unknown_edge_kos=frozenset(unknown_only_kos),
    )
    return genome, hits


def parse_per_mag_tsv(
    tsv_path: str,
    mag_id: Optional[str] = None,
    domain: str = "Bacteria",
    completeness: float = 100.0,
    contamination: float = 0.0,
    genome_size_bp: int = 0,
    contig_lengths: Sequence[int] = (),
    gff_cds_order: Optional[Dict[str, Tuple[str, ...]]] = None,
    complete_cds_table: bool = False,
) -> Tuple[GenomeAnnotations, List[GeneHit]]:
    with open(tsv_path, "r", encoding="utf-8") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        rows = list(reader)
    inferred_id = mag_id or tsv_path.rsplit("/", 1)[-1].replace(".tsv", "")
    return parse_per_mag_rows(
        rows=rows,
        mag_id=inferred_id,
        domain=domain,
        completeness=completeness,
        contamination=contamination,
        genome_size_bp=genome_size_bp,
        contig_lengths=contig_lengths,
        gff_cds_order=gff_cds_order,
        complete_cds_table=complete_cds_table,
    )


# ---------------------------------------------------------------------------
# 3. Well-Founded Horn Semantics, P-All Polynomial Decider & Graph Invariants
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class OrganismProfile:
    """Metabolic exchange profile of one coexisting living organism in event e."""
    org_id: str
    de_novo: FrozenSet[str] = frozenset()
    conversions: Tuple[Tuple[FrozenSet[str], str], ...] = ()
    required: FrozenSet[str] = frozenset()
    exportable: Optional[FrozenSet[str]] = None

    def rules(self) -> List[Tuple[FrozenSet[str], str]]:
        out: List[Tuple[FrozenSet[str], str]] = [
            (frozenset(), m) for m in sorted(self.de_novo)
        ]
        out.extend(self.conversions)
        return out

    def can_export(self, compound: str) -> bool:
        if self.exportable is None:
            return True
        return compound in self.exportable


def horn_closure_org(
    org: OrganismProfile,
    initial_avail: Iterable[str],
) -> Tuple[Set[str], Set[str]]:
    """Single-organism Horn least fixpoint starting from `initial_avail`.

    Returns `(avail, rule_produced)` where a rule `(ins -> p)` fires only when
    `ins <= avail`, adding `p` to `rule_produced` and `avail`.
    """
    avail: Set[str] = set(initial_avail)
    produced: Set[str] = set()
    rules = org.rules()
    changed = True
    while changed:
        changed = False
        for ins, prod in rules:
            if prod not in produced and ins.issubset(avail):
                produced.add(prod)
                if prod not in avail:
                    avail.add(prod)
                changed = True
    return avail, produced


def can_jointly_produce_and_live(
    org: OrganismProfile,
    exported_subset: FrozenSet[str],
    imported_pool: FrozenSet[str],
    ambient: FrozenSet[str],
) -> bool:
    """Exact well-foundedness test for simultaneous rule-production of `exported_subset`.

    Theorem: `org` has a single acyclic well-founded derivation (with one justification
    per `(org, m)` pair) in which EVERY `m in exported_subset` is justified by a rule of
    `org` (never imported, blocking both 1-metabolite `p -> q -> p` and 2-metabolite
    `p1 <-> p2` laundering; Astra 248969335645061120) AND `org.required` is met IFF
    `exported_subset <= rule_produced` and `org.required <= avail` when starting from
    `(ambient | imported_pool) - exported_subset`.
    """
    if any(not org.can_export(m) for m in exported_subset):
        return False
    init = (set(ambient) | set(imported_pool)) - set(exported_subset)
    avail, rule_prods = horn_closure_org(org, init)
    return exported_subset.issubset(rule_prods) and org.required.issubset(avail)


def is_metabolite_dep_acyclic(
    profiles: Dict[str, OrganismProfile],
    ambient: FrozenSet[str] = frozenset(),
) -> bool:
    """Return True iff the global metabolite rule-dependency graph is a DAG (Lemma L5).

    Directed edge `a -> p` for every rule `(ins -> p)` in any organism and `a in ins \\ ambient`.
    All P2 pathway rules (`COB_RING, DMB -> B12`, `THI_HMP, THI_HET -> THI`, and de novo
    endpoints) have topological depth <= 1 and are strictly acyclic.
    """
    adj: Dict[str, Set[str]] = {}
    for org in profiles.values():
        for ins, prod in org.rules():
            for a in ins - ambient:
                adj.setdefault(a, set()).add(prod)
    state: Dict[str, int] = {}

    def _dfs(u: str) -> bool:
        state[u] = 1
        for v in adj.get(u, ()):
            st_v = state.get(v, 0)
            if st_v == 1 or (st_v == 0 and not _dfs(v)):
                return False
        state[u] = 2
        return True

    return all(_dfs(u) for u in sorted(adj) if u not in state)


def two_node_closure(
    profiles: Dict[str, OrganismProfile],
    members: Iterable[str],
    ambient: FrozenSet[str],
    allow_edge: Optional[object] = None,
    extra_imports: Optional[Dict[str, FrozenSet[str]]] = None,
) -> Tuple[Dict[str, Set[str]], Dict[str, Set[str]]]:
    """Two-node Horn least fixpoint (`avail` and `prod` separate; producer-only exports).

    Exact for community viability on ALL instances (Lemma L1'), and exact for restricted
    supply-graph questions whenever `is_metabolite_dep_acyclic(profiles, ambient)` holds
    (Lemma L5, `discussion/deep/plexitrophy-computation.md`).
    """
    mem_list = sorted(members)
    extra = extra_imports or {}
    avail_by: Dict[str, Set[str]] = {
        y: set(ambient) | set(extra.get(y, ())) for y in mem_list
    }
    prod_by: Dict[str, Set[str]] = {y: set() for y in mem_list}

    changed = True
    while changed:
        changed = False
        for y in mem_list:
            for ins, prod in profiles[y].rules():
                if prod not in prod_by[y] and ins.issubset(avail_by[y]):
                    prod_by[y].add(prod)
                    avail_by[y].add(prod)
                    changed = True
        for y in mem_list:
            for z in mem_list:
                if z == y:
                    continue
                org_z = profiles[z]
                for m in prod_by[z] - avail_by[y]:
                    if org_z.can_export(m) and (
                        allow_edge is None or allow_edge(z, y, m)  # type: ignore[operator]
                    ):
                        avail_by[y].add(m)
                        changed = True
    return avail_by, prod_by


def two_node_smax(
    profiles: Dict[str, OrganismProfile],
    members: Iterable[str],
    ambient: FrozenSet[str],
    allow_edge: Optional[object] = None,
    extra_imports: Optional[Dict[str, FrozenSet[str]]] = None,
) -> Tuple[Set[str], Dict[str, Set[str]], Dict[str, Set[str]]]:
    """Greatest feasible sub-community S_max (Lemma L2)."""
    cur = set(members)
    while True:
        avail_by, prod_by = two_node_closure(
            profiles, cur, ambient, allow_edge=allow_edge, extra_imports=extra_imports
        )
        dead = {y for y in cur if not profiles[y].required.issubset(avail_by[y])}
        if not dead:
            return cur, avail_by, prod_by
        cur -= dead


def _decide_nonplexi_acyclic_poly(
    focal_id: str,
    profiles: Dict[str, OrganismProfile],
    amb: FrozenSet[str],
    cyc_mode: str,
) -> Tuple[str, Optional[Tuple[str, ...]]]:
    """Exact polynomial P-all decider when `is_metabolite_dep_acyclic` holds (Theorem 3).

    Performs at most `2|V| + 2` `two_node_smax` calls with ZERO bundle or support enumeration.
    """
    all_ids = set(profiles.keys())
    v_alive, _, _ = two_node_smax(profiles, all_ids, amb)
    if focal_id not in v_alive:
        return ("no_donor", None)

    if focal_id in two_node_smax(profiles, {focal_id}, amb)[0]:
        return ("d0", ("d0",))

    others = sorted(v_alive - {focal_id})

    # Pure sink check (d_X = 1, out-degree(X) = 0; valid across 'fanin_only', 'cyc', 'scc')
    c_set, _, prod_c = two_node_smax(profiles, others, amb)
    for d_id in others:
        if d_id in c_set:
            pool_d = frozenset(m for m in prod_c[d_id] if profiles[d_id].can_export(m))
            if (
                focal_id
                in two_node_smax(
                    profiles, {focal_id}, amb, extra_imports={focal_id: pool_d}
                )[0]
            ):
                return ("not_plexi", ("sink", d_id))

    if cyc_mode == "fanin_only":
        for d_id in others:
            allow_f1 = lambda z, y, m, d=d_id: (y != focal_id or z == d)
            if focal_id in two_node_smax(profiles, v_alive, amb, allow_edge=allow_f1)[0]:
                return ("not_plexi", ("fanin1", d_id))
        return ("plexi", None)

    if cyc_mode == "cyc":
        for d_id in others:
            allow_priv = lambda z, y, m, d=d_id: (y != focal_id or z == d) and (
                z != focal_id or y == d
            )
            if focal_id in two_node_smax(profiles, v_alive, amb, allow_edge=allow_priv)[0]:
                return ("not_plexi", ("private_dyad", d_id))
        return ("plexi", None)

    if cyc_mode == "scc":
        for d_id in others:
            u_set, _, prod_u = two_node_smax(profiles, set(others) - {d_id}, amb)
            pool_u = frozenset(
                m for z in u_set for m in prod_u[z] if profiles[z].can_export(m)
            )
            s2, _, _ = two_node_smax(
                profiles, {focal_id, d_id}, amb, extra_imports={d_id: pool_u}
            )
            if focal_id in s2 and d_id in s2:
                return ("not_plexi", ("dyad_scc", d_id))
        return ("plexi", None)

    raise ValueError(f"Unknown cyc_mode: {cyc_mode}")


def _decide_nonplexi_cyclic_milp(
    focal_id: str,
    profiles: Dict[str, OrganismProfile],
    amb: FrozenSet[str],
    cyc_mode: str,
    time_limit: Optional[float] = 60.0,
) -> Tuple[str, Optional[Tuple[str, ...]]]:
    """Exact single-node S1 counterexample MILP for cyclic metabolite dependencies (Theorem 3').

    Encodes strict proof-depth variables `t_{Y, m}` per `(organism, metabolite)` node and
    producer-only exports (`x_{Z, Y, m} <= sum_k j_{Z, k}`), then self-certifies any returned
    witness as an explicit acyclic single-justification proof DAG.
    """
    if cyc_mode not in ("fanin_only", "cyc", "scc"):
        raise ValueError(f"Unknown cyc_mode: {cyc_mode}")

    # Lemma L1': unrestricted community viability under S1 equals two_node_smax
    v_alive, _, _ = two_node_smax(profiles, set(profiles.keys()), amb)
    if focal_id not in v_alive:
        return ("no_donor", None)
    if focal_id in two_node_smax(profiles, {focal_id}, amb)[0]:
        return ("d0", ("d0",))

    try:
        import numpy as np
        from scipy.optimize import Bounds, LinearConstraint, milp
        from scipy.sparse import coo_matrix
    except ImportError:
        return ("unresolved_cyclic_no_milp", None)

    v_list = sorted(v_alive)
    selfcl = {y: two_node_closure(profiles, [y], amb)[0][y] for y in v_list}
    uses: Dict[str, Set[str]] = {}
    for y in v_list:
        u_set = set(profiles[y].required)
        for ins, _ in profiles[y].rules():
            u_set |= set(ins)
        uses[y] = u_set - set(amb)

    rules_list: List[Tuple[str, int, FrozenSet[str], str]] = []
    prods: Dict[str, Set[str]] = {y: set() for y in v_list}
    for y in v_list:
        for k, (ins, prod) in enumerate(profiles[y].rules()):
            if prod in amb or prod in ins:
                continue
            clean_ins = frozenset(ins - amb)
            rules_list.append((y, k, clean_ins, prod))
            prods[y].add(prod)

    nodes = [(y, m) for y in v_list for m in sorted(uses[y] | prods[y])]
    big_n = float(len(nodes) + 1)

    idx: Dict[Tuple[object, ...], int] = {}
    lb: List[float] = []
    ub: List[float] = []
    integ: List[int] = []

    def _var(key: Tuple[object, ...], lo: float = 0.0, hi: float = 1.0, is_int: int = 1) -> int:
        idx[key] = len(lb)
        lb.append(lo)
        ub.append(hi)
        integ.append(is_int)
        return idx[key]

    for y in v_list:
        _var(("s", y))
    for y, m in nodes:
        _var(("a", y, m))
        _var(("t", y, m), 0.0, big_n, 0)
    for y, k, _, _ in rules_list:
        _var(("j", y, k))

    imps: List[Tuple[str, str, str]] = []
    for z in v_list:
        org_z = profiles[z]
        exp_z = {m for m in prods[z] if org_z.can_export(m)}
        for y in v_list:
            if y == z:
                continue
            for m in sorted(exp_z & uses[y]):
                if m in selfcl[y]:
                    continue
                imps.append((z, y, m))
                _var(("x", z, y, m))

    hedges = sorted({(z, y) for z, y, _ in imps})
    for z, y in hedges:
        _var(("h", z, y))
    if cyc_mode in ("cyc", "scc"):
        for z in v_list:
            if z != focal_id:
                _var(("r", z), 0.0, 1.0, 0)

    rows: List[int] = []
    cols: List[int] = []
    vals: List[float] = []
    rlo: List[float] = []
    rhi: List[float] = []

    def _con(coefs: Sequence[Tuple[Tuple[object, ...], float]], lo: float, hi: float) -> None:
        r_idx = len(rlo)
        for k_var, val in coefs:
            rows.append(r_idx)
            cols.append(idx[k_var])
            vals.append(val)
        rlo.append(lo)
        rhi.append(hi)

    inf = float(np.inf)
    _con([(("s", focal_id), 1.0)], 1.0, 1.0)

    just_vars: Dict[Tuple[str, str], List[Tuple[object, ...]]] = {n: [] for n in nodes}
    rules_by_ym: Dict[Tuple[str, str], List[Tuple[int, FrozenSet[str]]]] = {}

    for y, k, ins, prod in rules_list:
        just_vars[(y, prod)].append(("j", y, k))
        rules_by_ym.setdefault((y, prod), []).append((k, ins))
        _con([(("j", y, k), 1.0), (("s", y), -1.0)], -inf, 0.0)
        for i_met in ins:
            _con([(("j", y, k), 1.0), (("a", y, i_met), -1.0)], -inf, 0.0)
            _con(
                [(("t", y, prod), 1.0), (("t", y, i_met), -1.0), (("j", y, k), -big_n)],
                1.0 - big_n,
                inf,
            )

    imps_by_ym: Dict[Tuple[str, str], List[str]] = {}
    for z, y, m in imps:
        just_vars[(y, m)].append(("x", z, y, m))
        imps_by_ym.setdefault((y, m), []).append(z)
        z_rules = [("j", z, k) for (k, _) in rules_by_ym.get((z, m), ())]
        _con([(("x", z, y, m), 1.0)] + [(zr, -1.0) for zr in z_rules], -inf, 0.0)
        _con([(("x", z, y, m), 1.0), (("s", z), -1.0)], -inf, 0.0)
        _con([(("x", z, y, m), 1.0), (("s", y), -1.0)], -inf, 0.0)
        _con(
            [(("t", y, m), 1.0), (("t", z, m), -1.0), (("x", z, y, m), -big_n)],
            1.0 - big_n,
            inf,
        )
        _con([(("h", z, y), 1.0), (("x", z, y, m), -1.0)], 0.0, inf)

    for y, m in nodes:
        _con([(("a", y, m), 1.0)] + [(jv, -1.0) for jv in just_vars[(y, m)]], -inf, 0.0)

    for y in v_list:
        for req_m in profiles[y].required - amb:
            _con([(("s", y), 1.0), (("a", y, req_m), -1.0)], -inf, 0.0)

    in_x = [z for z, y in hedges if y == focal_id]
    if in_x:
        _con([(("h", z, focal_id), 1.0) for z in in_x], -inf, 1.0)

    hs_set = set(hedges)
    if cyc_mode in ("cyc", "scc"):
        for z, y in hedges:
            if y == focal_id:
                continue
            if z == focal_id:
                if cyc_mode == "cyc":
                    extra_term = [(("h", y, focal_id), 1.0)] if y in in_x else []
                    _con([(("r", y), 1.0), (("h", focal_id, y), -1.0)] + extra_term, 0.0, inf)
                else:
                    _con([(("r", y), 1.0), (("h", focal_id, y), -1.0)], 0.0, inf)
            else:
                _con([(("r", y), 1.0), (("r", z), -1.0), (("h", z, y), -1.0)], -1.0, inf)

        for z in v_list:
            if z == focal_id:
                continue
            if cyc_mode == "cyc":
                if z in in_x:
                    _con([(("r", z), 1.0), (("h", z, focal_id), 1.0)], -inf, 1.0)
            else:
                for d_cand in in_x:
                    if d_cand != z and (z, d_cand) in hs_set:
                        _con(
                            [(("h", z, d_cand), 1.0), (("r", z), 1.0), (("h", d_cand, focal_id), 1.0)],
                            -inf,
                            2.0,
                        )

    n_vars = len(lb)
    mat = coo_matrix((vals, (rows, cols)), shape=(len(rlo), n_vars)).tocsr()
    obj = np.zeros(n_vars)
    for y in v_list:
        obj[idx[("s", y)]] = 1.0
    opts = {} if time_limit is None else {"time_limit": float(time_limit)}
    res = milp(
        obj,
        constraints=LinearConstraint(mat, rlo, rhi),
        integrality=np.array(integ),
        bounds=Bounds(lb, ub),
        options=opts,
    )

    if res.status == 2:
        return ("plexi", None)
    if res.x is None:
        return ("unresolved_timeout", None)

    sol = res.x
    # Self-certification: extract and verify an explicit single-node S1 proof DAG
    alive = {y for y in v_list if sol[idx[("s", y)]] > 0.5}
    assert focal_id in alive, "MILP witness omitted focal organism"

    chosen_just: Dict[Tuple[str, str], Tuple[object, ...]] = {}
    t_map: Dict[Tuple[str, str], float] = {}
    for y in alive:
        for m in amb:
            chosen_just[(y, m)] = ("ext",)
            t_map[(y, m)] = -1.0
    for y, m in nodes:
        if y not in alive:
            continue
        t_ym = float(sol[idx[("t", y, m)]])
        t_map[(y, m)] = t_ym
        active_rules = [
            (k, ins)
            for (k, ins) in rules_by_ym.get((y, m), ())
            if sol[idx[("j", y, k)]] > 0.5
        ]
        if active_rules:
            _, ins = active_rules[0]
            chosen_just[(y, m)] = ("rule", ins)
            continue
        active_imps = [
            z
            for z in imps_by_ym.get((y, m), ())
            if z in alive and sol[idx[("x", z, y, m)]] > 0.5
        ]
        if active_imps:
            chosen_just[(y, m)] = ("import", active_imps[0])

    # Backward cone verification from focal_id
    cone_orgs: Set[str] = {focal_id}
    seen_pairs: Set[Tuple[str, str]] = set()
    p_queue: List[Tuple[str, str]] = [
        (focal_id, r_m) for r_m in sorted(profiles[focal_id].required)
    ]
    supply_edges: Set[Tuple[str, str]] = set()

    while p_queue:
        ym = p_queue.pop()
        if ym in seen_pairs:
            continue
        assert ym in chosen_just, f"MILP witness left {ym} unjustified"
        seen_pairs.add(ym)
        y, m = ym
        j = chosen_just[ym]
        if j[0] == "ext":
            assert m in amb
        elif j[0] == "rule":
            ins = j[1]  # type: ignore[assignment]
            for a_m in ins:  # type: ignore[union-attr]
                assert t_map[(y, a_m)] <= t_map[ym] - 0.5, "Proof-depth rule violation"
                p_queue.append((y, a_m))
        elif j[0] == "import":
            z = str(j[1])
            assert z != y and z in alive and profiles[z].can_export(m)
            assert chosen_just.get((z, m), ("none",))[0] == "rule", "Relay violation in witness"
            assert t_map[(z, m)] <= t_map[ym] - 0.5, "Proof-depth import violation"
            supply_edges.add((z, y))
            p_queue.append((z, m))
            if z not in cone_orgs:
                cone_orgs.add(z)
                for req_z in sorted(profiles[z].required):
                    p_queue.append((z, req_z))

    d_suppliers = sorted({u for (u, v) in supply_edges if v == focal_id and u != focal_id})
    assert len(d_suppliers) <= 1, f"MILP witness had d_X={len(d_suppliers)} > 1"
    if cyc_mode == "scc":
        sccs = tarjan_sccs(cone_orgs, supply_edges)
        c_x = next((len(comp) for comp in sccs if focal_id in comp), 1)
        assert c_x <= 2, f"MILP witness violated |SCC(X)|<=2: {c_x}"
    elif cyc_mode == "cyc":
        assert not on_simple_cycle_ge_3(cone_orgs, supply_edges, focal_id), (
            "MILP witness violated simple-cycle kappa(X)==0"
        )

    if not d_suppliers:
        return ("d0", ("d0",))
    d_id = d_suppliers[0]
    has_out = any(u == focal_id for (u, _) in supply_edges)
    return ("not_plexi", ("sink" if not has_out else "milp_d1", d_id))


def decide_nonplexi_p_all(
    focal_id: str,
    profiles: Dict[str, OrganismProfile],
    ambient: Iterable[str] = (),
    cyc_mode: str = "cyc",
    time_limit: Optional[float] = 60.0,
) -> Tuple[str, Optional[Tuple[str, ...]]]:
    """Exact P-all single-node well-founded decider (Opus 248974707227279360 & Astra 248974912022560768).

    - Checks `is_metabolite_dep_acyclic(profiles, amb)` first:
      * If True (holds for 100% of P2 pathway rules): runs `_decide_nonplexi_acyclic_poly`
        (Theorem 3 + Lemma L5; `O(|V|)` `two_node_smax` calls, pure Python stdlib, zero
        bundle or support enumeration).
      * If False (cyclic metabolite conversion rules; coNP-complete by Theorem 3'): runs
        `_decide_nonplexi_cyclic_milp` (exact proof-depth MILP with single-node S1 DAG
        self-certification).
    """
    amb = frozenset(ambient)
    if focal_id not in profiles:
        raise KeyError(f"Unknown focal organism {focal_id}")
    if is_metabolite_dep_acyclic(profiles, amb):
        return _decide_nonplexi_acyclic_poly(focal_id, profiles, amb, cyc_mode)
    return _decide_nonplexi_cyclic_milp(focal_id, profiles, amb, cyc_mode, time_limit=time_limit)


@dataclass(frozen=True)
class SupportWiring:
    """Feasible minimal support + directed living-partner wiring (for small-graph inspection)."""
    focal_id: str
    members: FrozenSet[str]
    partner_edges: FrozenSet[Tuple[str, str]]
    transfers: FrozenSet[Tuple[str, str, str]] = frozenset()

    def invariants(self, target_id: Optional[str] = None) -> Tuple[int, int, int]:
        """Return (d, m, c) where c = |SCC(target_id)|."""
        x = target_id or self.focal_id
        if x not in self.members:
            raise ValueError(f"{x} not in support {self.members}")
        d_x = len({u for (u, v) in self.partner_edges if v == x and u != x})
        m_x = len(self.members)
        sccs = tarjan_sccs(self.members, self.partner_edges)
        c_x = next((len(comp) for comp in sccs if x in comp), 1)
        return d_x, m_x, c_x

    def kappa(self, target_id: Optional[str] = None) -> int:
        """Return 1 iff target_id lies on a simple directed cycle of length >= 3."""
        x = target_id or self.focal_id
        return 1 if on_simple_cycle_ge_3(self.members, self.partner_edges, x) else 0

    def is_plexitrophic(self, target_id: Optional[str] = None, cyc_mode: str = "cyc") -> bool:
        d_x, _, c_x = self.invariants(target_id)
        if d_x >= 2:
            return True
        if cyc_mode == "scc":
            return c_x >= 3
        if cyc_mode == "cyc":
            return self.kappa(target_id) == 1
        return False


def on_simple_cycle_ge_3(
    nodes: Iterable[str],
    edges: Iterable[Tuple[str, str]],
    focal_id: str,
) -> bool:
    """Return True iff `focal_id` lies on a simple directed cycle of length >= 3."""
    node_set = set(nodes)
    adj: Dict[str, Set[str]] = {n: set() for n in node_set}
    for u, v in edges:
        if u != v and u in node_set and v in node_set:
            adj[u].add(v)

    def _dfs(cur: str, visited: Set[str]) -> bool:
        for nxt in adj.get(cur, ()):
            if nxt == focal_id and len(visited) >= 3:
                return True
            if nxt not in visited and _dfs(nxt, visited | {nxt}):
                return True
        return False

    return _dfs(focal_id, {focal_id})


def tarjan_sccs(
    nodes: Iterable[str],
    edges: Iterable[Tuple[str, str]],
) -> List[FrozenSet[str]]:
    """Exact Tarjan strongly connected components on a directed graph."""
    adj: Dict[str, List[str]] = {n: [] for n in nodes}
    for u, v in edges:
        if u in adj and v in adj:
            adj[u].append(v)

    index = 0
    indices: Dict[str, int] = {}
    lowlink: Dict[str, int] = {}
    stack: List[str] = []
    on_stack: Set[str] = set()
    sccs: List[FrozenSet[str]] = []

    def strongconnect(v: str) -> None:
        nonlocal index
        indices[v] = index
        lowlink[v] = index
        index += 1
        stack.append(v)
        on_stack.add(v)
        for w in adj[v]:
            if w not in indices:
                strongconnect(w)
                lowlink[v] = min(lowlink[v], lowlink[w])
            elif w in on_stack:
                lowlink[v] = min(lowlink[v], indices[w])
        if lowlink[v] == indices[v]:
            comp: Set[str] = set()
            while True:
                w = stack.pop()
                on_stack.remove(w)
                comp.add(w)
                if w == v:
                    break
            sccs.append(frozenset(comp))

    for n in sorted(adj):
        if n not in indices:
            strongconnect(n)
    return sccs


def is_wiring_feasible(
    focal_id: str,
    members: FrozenSet[str],
    transfers: FrozenSet[Tuple[str, str, str]],
    profiles: Dict[str, OrganismProfile],
    ambient: FrozenSet[str],
) -> bool:
    """Verify strict single-node S1 groundedness, producer-only exports, viability, and reachability."""
    if focal_id not in members:
        return False
    exports_from: Dict[str, Set[str]] = {u: set() for u in members}
    imports_into: Dict[str, Set[str]] = {u: set() for u in members}
    for u, v, cmpd in transfers:
        if u == v or u not in members or v not in members or cmpd in ambient:
            return False
        if not profiles[u].can_export(cmpd):
            return False
        exports_from[u].add(cmpd)
        imports_into[v].add(cmpd)

    # Under S1 (one justification per (u, m)), no organism can both import and export m
    for u in members:
        if exports_from[u] & imports_into[u]:
            return False

    avail: Dict[str, Set[str]] = {u: set(ambient) for u in members}
    rule_prod: Dict[str, Set[str]] = {u: set() for u in members}
    changed = True
    while changed:
        changed = False
        for u in members:
            for ins, prod in profiles[u].rules():
                if (
                    prod not in imports_into[u]
                    and prod not in rule_prod[u]
                    and ins.issubset(avail[u])
                ):
                    rule_prod[u].add(prod)
                    avail[u].add(prod)
                    changed = True
        for u, v, cmpd in transfers:
            if cmpd in rule_prod[u] and cmpd not in avail[v]:
                avail[v].add(cmpd)
                changed = True

    for u, _, cmpd in transfers:
        if cmpd not in rule_prod[u]:
            return False
    for u in members:
        if not profiles[u].required.issubset(avail[u]):
            return False

    rev_adj: Dict[str, Set[str]] = {u: set() for u in members}
    for u, v, _ in transfers:
        rev_adj[v].add(u)
    seen = {focal_id}
    queue = [focal_id]
    while queue:
        cur = queue.pop()
        for prev in rev_adj[cur]:
            if prev not in seen:
                seen.add(prev)
                queue.append(prev)
    return seen == set(members)


def enumerate_minimal_wirings(
    focal_id: str,
    profiles: Dict[str, OrganismProfile],
    ambient: Iterable[str] = (),
    max_support_size: Optional[int] = None,
) -> List[SupportWiring]:
    """Enumerate inclusion-minimal wirings (for small fixture inspection; NOT used for P2 decisions)."""
    amb = frozenset(ambient)
    solo = frozenset({focal_id})
    if is_wiring_feasible(focal_id, solo, frozenset(), profiles, amb):
        return [SupportWiring(focal_id=focal_id, members=solo, partner_edges=frozenset())]

    potentially_produced: Dict[str, Set[str]] = {}
    for uid, p in profiles.items():
        prods = set(p.de_novo) | {out for _, out in p.conversions}
        potentially_produced[uid] = {c for c in prods if p.can_export(c)}

    potentially_needed: Dict[str, Set[str]] = {}
    for uid, p in profiles.items():
        needs = set(p.required)
        for req_in, _ in p.conversions:
            needs |= set(req_in)
        potentially_needed[uid] = needs - amb

    cand_transfers_by_pair: Dict[Tuple[str, str], List[str]] = {}
    for u in profiles:
        for v in profiles:
            if u != v:
                shared = sorted(potentially_produced[u] & potentially_needed[v])
                if shared:
                    cand_transfers_by_pair[(u, v)] = shared

    other_ids = sorted(uid for uid in profiles if uid != focal_id)
    cap = len(other_ids) if max_support_size is None else min(len(other_ids), max_support_size - 1)
    found_wirings: List[SupportWiring] = []

    for k in range(1, cap + 1):
        for comb in itertools.combinations(other_ids, k):
            members = frozenset((focal_id, *comb))
            pair_transfers = [
                (u, v, c)
                for u in members
                for v in members
                if u != v and (u, v) in cand_transfers_by_pair
                for c in cand_transfers_by_pair[(u, v)]
            ]
            if not pair_transfers:
                continue
            minimal_t_sets: List[FrozenSet[Tuple[str, str, str]]] = []
            for r in range(len(members) - 1, len(pair_transfers) + 1):
                for t_comb in itertools.combinations(pair_transfers, r):
                    t_set = frozenset(t_comb)
                    if any(prev.issubset(t_set) for prev in minimal_t_sets):
                        continue
                    if is_wiring_feasible(focal_id, members, t_set, profiles, amb):
                        minimal_t_sets.append(t_set)
            for t_set in minimal_t_sets:
                p_edges = frozenset((u, v) for (u, v, _) in t_set)
                found_wirings.append(
                    SupportWiring(
                        focal_id=focal_id,
                        members=members,
                        partner_edges=p_edges,
                        transfers=t_set,
                    )
                )

    unique_map: Dict[Tuple[FrozenSet[str], FrozenSet[Tuple[str, str]]], SupportWiring] = {}
    for w in found_wirings:
        key = (w.members, w.partner_edges)
        if key not in unique_map or len(w.transfers) < len(unique_map[key].transfers):
            unique_map[key] = w

    all_w = list(unique_map.values())
    return [
        w
        for w in all_w
        if not any(
            other is not w
            and other.members.issubset(w.members)
            and other.partner_edges.issubset(w.partner_edges)
            and (other.members != w.members or other.partner_edges != w.partner_edges)
            for other in all_w
        )
    ]


# ---------------------------------------------------------------------------
# 4. Genome -> P2 Metabolic Profile & Exact P-All MAG Entry Point
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class SupplyScenario:
    name: str
    ambient_compounds: FrozenSet[str] = frozenset()


SCENARIO_E_MINIMAL = SupplyScenario("E_minimal", frozenset())
SCENARIO_E_PLUS_DIN = SupplyScenario("E_plus_DIN", frozenset({"FIXED_N"}))
SCENARIO_E_PLUS_DMSP = SupplyScenario("E_plus_DMSP", frozenset({"RED_S"}))
SCENARIO_E_PLUS_VITAMINS = SupplyScenario("E_plus_vitamins", frozenset({"B12", "THI"}))


def recipient_endpoint_status(
    g: GenomeAnnotations,
    endpoint: str,
    max_false_absence_prob: Optional[float] = 0.05,
) -> str:
    """Evaluate focal genome status for `endpoint` BEFORE inspecting donors.

    Requires calibrated detectability (`p_false_absence < max_false_absence_prob`,
    default 0.05) before calling any clean `STATE_ABSENT` gap (Astra 248968746164994048 #4).
    """
    if endpoint == "B12":
        return evaluate_recipient_b12_auxotrophy(g, max_false_absence_prob=max_false_absence_prob)
    if endpoint == "RED_S":
        return PATHWAY_RULES["ASSIM_SO4"].evaluate(g, max_false_absence_prob=max_false_absence_prob).state
    if endpoint == "FIXED_N":
        nif = PATHWAY_RULES["A_NIF"].evaluate(g, max_false_absence_prob=max_false_absence_prob).state
        no3 = PATHWAY_RULES["ASSIM_NO3"].evaluate(g, max_false_absence_prob=max_false_absence_prob).state
        if nif == STATE_PRESENT or no3 == STATE_PRESENT:
            return STATE_PRESENT
        if nif == STATE_ABSENT and no3 == STATE_ABSENT:
            return STATE_ABSENT
        return STATE_UNRESOLVED
    if endpoint == "THI":
        return PATHWAY_RULES["C_THI_SPLIT"].evaluate(g, max_false_absence_prob=max_false_absence_prob).state
    if endpoint in PATHWAY_RULES:
        return PATHWAY_RULES[endpoint].evaluate(g, max_false_absence_prob=max_false_absence_prob).state
    raise KeyError(f"Unknown endpoint {endpoint}")


def build_mag_organism_profile(
    g: GenomeAnnotations,
    required_override: Optional[FrozenSet[str]] = None,
    modeled_donor_needs: Sequence[str] = (),
    optimistic_capabilities: bool = False,
    max_false_absence_prob: Optional[float] = 0.05,
    include_c0: bool = True,
    exportable_override: Optional[FrozenSet[str]] = None,
) -> OrganismProfile:
    """Construct `OrganismProfile` for conservative (`I^-`) or optimistic (`I^+`) completion."""
    target_states = {STATE_PRESENT, STATE_UNRESOLVED} if optimistic_capabilities else {STATE_PRESENT}
    de_novo: Set[str] = set()
    conversions: List[Tuple[FrozenSet[str], str]] = []

    aer_ring = PATHWAY_RULES["P_B12_AER"].evaluate(g, max_false_absence_prob).state
    ana_ring = PATHWAY_RULES["A_B12_ANA"].evaluate(g, max_false_absence_prob).state
    blub_full = PATHWAY_RULES["P_BLUB"].evaluate(g, max_false_absence_prob).state
    blub_dmb = PATHWAY_RULES["P_BLUB_DMB"].evaluate(g, max_false_absence_prob).state
    bza = PATHWAY_RULES["A_BZA"].evaluate(g, max_false_absence_prob).state
    attach = PATHWAY_RULES["B12_LOWER_ATTACH"].evaluate(g, max_false_absence_prob).state

    ring_ok = (aer_ring in target_states) or (ana_ring in target_states)
    lig_ok = (blub_full in target_states) or (bza in target_states)

    if ring_ok:
        de_novo.add("COB_RING")
    if blub_dmb in target_states:
        de_novo.add("DMB")
    if ring_ok and lig_ok:
        de_novo.add("B12")
    if attach in target_states:
        conversions.append((frozenset({"COB_RING", "DMB"}), "B12"))

    hmp = PATHWAY_RULES["THI_HMP"].evaluate(g, max_false_absence_prob).state
    het = PATHWAY_RULES["THI_HET"].evaluate(g, max_false_absence_prob).state
    tcpl = PATHWAY_RULES["THI_COUPLE"].evaluate(g, max_false_absence_prob).state
    if hmp in target_states:
        de_novo.add("THI_HMP")
    if het in target_states:
        de_novo.add("THI_HET")
    if hmp in target_states and het in target_states and tcpl in target_states:
        de_novo.add("THI")
    elif tcpl in target_states:
        conversions.append((frozenset({"THI_HMP", "THI_HET"}), "THI"))

    if PATHWAY_RULES["A_DSR"].evaluate(g, max_false_absence_prob).state in target_states:
        de_novo.add("RED_S")
    if PATHWAY_RULES["A_NIF"].evaluate(g, max_false_absence_prob).state in target_states:
        de_novo.add("FIXED_N")

    if include_c0:
        for mod_id in C0_CONTROL_MODULES:
            if PATHWAY_RULES[mod_id].evaluate(g, max_false_absence_prob).state in target_states:
                de_novo.add(mod_id)

    if required_override is not None:
        reqs = required_override
    else:
        need_set: Set[str] = set()
        for ep in modeled_donor_needs:
            st = recipient_endpoint_status(g, ep, max_false_absence_prob=max_false_absence_prob)
            if optimistic_capabilities:
                # In I^+ (fewer needs), only confident STATE_ABSENT counts as a donor need
                if st == STATE_ABSENT:
                    need_set.add(ep)
            else:
                # In I^- (more needs), both STATE_ABSENT and STATE_UNRESOLVED count as needs
                if st in (STATE_ABSENT, STATE_UNRESOLVED):
                    need_set.add(ep)
        reqs = frozenset(need_set)

    return OrganismProfile(
        org_id=g.mag_id,
        de_novo=frozenset(de_novo),
        conversions=tuple(conversions),
        required=reqs,
        exportable=exportable_override,
    )


def _build_relaxed_pool_pseudo_organism(
    pool_genome: GenomeAnnotations,
    org_id: str,
    exportable_override: Optional[FrozenSet[str]] = None,
) -> OrganismProfile:
    """Proposition 9 relaxed pseudo-organism `P*` for `I^+` only (never produces resolved d=1)."""
    u_opt = build_mag_organism_profile(
        pool_genome,
        required_override=frozenset(),
        optimistic_capabilities=True,
        max_false_absence_prob=None,
    )
    u_prods = set(u_opt.de_novo) | {out for _, out in u_opt.conversions}
    if "COB_RING" in u_prods and "DMB" in u_prods:
        u_prods.add("B12")
    if "THI_HMP" in u_prods and "THI_HET" in u_prods:
        u_prods.add("THI")
    if exportable_override is not None:
        u_prods &= set(exportable_override)
    return OrganismProfile(
        org_id=org_id,
        de_novo=frozenset(u_prods),
        conversions=(),
        required=frozenset(),
        exportable=exportable_override,
    )


def evaluate_recipient_requirements_in_event(
    recipient: GenomeAnnotations,
    required_endpoints: Sequence[str],
    co_occurring_mags: Sequence[GenomeAnnotations],
    scenario: SupplyScenario = SCENARIO_E_MINIMAL,
    unbinned_genome: Optional[GenomeAnnotations] = None,
    mode: str = "fanin_only",
    donor_need_endpoints: Sequence[str] = (),
    max_false_absence_prob: Optional[float] = 0.05,
    eukaryotic_donors: Sequence[GenomeAnnotations] = (),
    unassigned_euk_genome: Optional[GenomeAnnotations] = None,
    euk_exportable_endpoints: Optional[FrozenSet[str]] = None,
    euk_annotated_routes_only: bool = False,
) -> Tuple[str, Optional[Tuple[str, ...]]]:
    """Score prokaryotic HQ recipient X in event e using two-point monotone P-all search.

    Modes:
      - `mode='fanin_only'` (default one-level P2 scope): donors have no modeled needs;
        tests exact fan-in (`d >= 2`) across ALL supports of ANY size (`O(|V|)` time).
      - `mode='consortium_cyc'`: focal X exports its own synthesized endpoints and
        co-occurring prokaryotic MAGs carry their own auxotrophies over `donor_need_endpoints`,
        testing `d >= 2 OR kappa(X) == 1` (simple cycle of length >= 3).

    Eukaryotic donor class (Opus 248976762356547584 & Astra 248976870271795200 / 248978604545196032):
      - Focal recipients remain strictly prokaryotic (`Bacteria`/`Archaea`).
      - `eukaryotic_donors`: organism-resolved eukaryotic MAGs/bins in event `e`.
        * In `I^-`: exports only routes with `STATE_PRESENT` in that single eukaryotic organism.
        * In `I^+` (default `euk_annotated_routes_only=False`): uses `max_false_absence_prob`
          (`0.05` default). Biologically excluded prokaryote-only routes (`P_B12_AER`, `A_B12_ANA`,
          `P_BLUB*`, `A_BZA`, `A_NIF`, `A_DSR`, etc.) evaluate to `STATE_ABSENT` (`p_fa = 0`),
          whereas eukaryote-eligible routes (`THI`, `C0` amino acids/cofactors, `B12_LOWER_ATTACH`)
          with uncalibrated/fragmentary eukaryotic detectability evaluate to `STATE_UNRESOLVED`
          and are RETAINED in `I^+`, downgrading putative fan-in to `unresolved_ambig_donor`
          unless eukaryotic absence detectability is explicitly calibrated (`false_absence_probs < 0.05`)
          or `euk_annotated_routes_only=True` (annotated-routes-only conditional bound).
      - `unassigned_euk_genome`: pooled multi-taxon MetaEuk/MarFERReT eukaryotic hits in `e`.
        Placed STRICTLY in `I^+` (`_euk_unassigned_P_star`), never manufacturing a resolved `d = 1`.
    """
    if recipient.domain not in ("Bacteria", "Archaea"):
        return CALL_NOT_ELIGIBLE, None

    reqs = tuple(dict.fromkeys(required_endpoints))
    statuses = [
        recipient_endpoint_status(recipient, ep, max_false_absence_prob=max_false_absence_prob)
        for ep in reqs
    ]
    if any(s == "not_dependent" for s in statuses):
        return CALL_NOT_ELIGIBLE, None
    if all(s == STATE_PRESENT for s in statuses):
        return CALL_AUTONOMOUS, None
    if not all(s == STATE_ABSENT for s in statuses):
        return CALL_NOT_ELIGIBLE, None

    cyc_mode = "fanin_only" if mode == "fanin_only" else "cyc"
    active_donor_needs: Sequence[str] = (
        ()
        if mode == "fanin_only"
        else (donor_need_endpoints if donor_need_endpoints else ALL_P2_ENDPOINTS)
    )

    prof_minus: Dict[str, OrganismProfile] = {}
    prof_plus: Dict[str, OrganismProfile] = {}

    prof_minus[recipient.mag_id] = build_mag_organism_profile(
        recipient,
        required_override=frozenset(reqs),
        optimistic_capabilities=False,
        max_false_absence_prob=max_false_absence_prob,
    )
    prof_plus[recipient.mag_id] = build_mag_organism_profile(
        recipient,
        required_override=frozenset(reqs),
        optimistic_capabilities=True,
        max_false_absence_prob=max_false_absence_prob,
    )

    for donor in co_occurring_mags:
        if donor.mag_id == recipient.mag_id:
            continue
        prof_minus[donor.mag_id] = build_mag_organism_profile(
            donor,
            modeled_donor_needs=active_donor_needs,
            optimistic_capabilities=False,
            max_false_absence_prob=max_false_absence_prob,
        )
        prof_plus[donor.mag_id] = build_mag_organism_profile(
            donor,
            modeled_donor_needs=active_donor_needs,
            optimistic_capabilities=True,
            max_false_absence_prob=max_false_absence_prob,
        )

    for euk in eukaryotic_donors:
        if euk.mag_id == recipient.mag_id:
            continue
        euk_g = (
            euk
            if euk.domain == "Eukaryota"
            else GenomeAnnotations(
                mag_id=euk.mag_id,
                domain="Eukaryota",
                completeness=euk.completeness,
                contamination=euk.contamination,
                genome_size_bp=euk.genome_size_bp,
                contig_lengths=euk.contig_lengths,
                kos=euk.kos,
                end_kos=euk.end_kos,
                unknown_edge_kos=euk.unknown_edge_kos,
                confirmed_absent_modules=euk.confirmed_absent_modules,
                false_absence_probs=euk.false_absence_probs,
            )
        )
        euk_fap_threshold = None if euk_annotated_routes_only else max_false_absence_prob
        prof_minus[euk_g.mag_id] = build_mag_organism_profile(
            euk_g,
            required_override=frozenset(),
            optimistic_capabilities=False,
            max_false_absence_prob=euk_fap_threshold,
            exportable_override=euk_exportable_endpoints,
        )
        prof_plus[euk_g.mag_id] = build_mag_organism_profile(
            euk_g,
            required_override=frozenset(),
            optimistic_capabilities=True,
            max_false_absence_prob=euk_fap_threshold,
            exportable_override=euk_exportable_endpoints,
        )

    if unbinned_genome is not None:
        prof_plus["_unbinned_P_star"] = _build_relaxed_pool_pseudo_organism(
            unbinned_genome, "_unbinned_P_star"
        )

    if unassigned_euk_genome is not None:
        unass_euk_g = (
            unassigned_euk_genome
            if unassigned_euk_genome.domain == "Eukaryota"
            else GenomeAnnotations(
                mag_id=unassigned_euk_genome.mag_id,
                domain="Eukaryota",
                completeness=unassigned_euk_genome.completeness,
                contamination=unassigned_euk_genome.contamination,
                genome_size_bp=unassigned_euk_genome.genome_size_bp,
                contig_lengths=unassigned_euk_genome.contig_lengths,
                kos=unassigned_euk_genome.kos,
                end_kos=unassigned_euk_genome.end_kos,
                unknown_edge_kos=unassigned_euk_genome.unknown_edge_kos,
                confirmed_absent_modules=unassigned_euk_genome.confirmed_absent_modules,
                false_absence_probs=unassigned_euk_genome.false_absence_probs,
            )
        )
        prof_plus["_euk_unassigned_P_star"] = _build_relaxed_pool_pseudo_organism(
            unass_euk_g,
            "_euk_unassigned_P_star",
            exportable_override=euk_exportable_endpoints,
        )

    st_minus, wit_minus = decide_nonplexi_p_all(
        recipient.mag_id, prof_minus, scenario.ambient_compounds, cyc_mode=cyc_mode
    )
    st_plus, wit_plus = decide_nonplexi_p_all(
        recipient.mag_id, prof_plus, scenario.ambient_compounds, cyc_mode=cyc_mode
    )

    if st_minus == "d0":
        return CALL_D0_AMBIENT, wit_minus
    if st_minus == "not_plexi":
        return CALL_D1_OR_CHAIN, wit_minus
    if st_plus == "no_donor":
        return CALL_UNRESOLVED_NO_DONOR, None
    if st_minus == "plexi" and st_plus == "plexi":
        return CALL_PLEXITROPHIC, None
    return CALL_UNRESOLVED_AMBIG, wit_plus


def compare_calls_with_eukaryotic_donors(
    evaluations: Sequence[
        Tuple[
            GenomeAnnotations,
            Sequence[str],
            Sequence[GenomeAnnotations],
            Sequence[GenomeAnnotations],
            Optional[GenomeAnnotations],
        ]
    ],
    scenario: SupplyScenario = SCENARIO_E_MINIMAL,
    euk_exportable_endpoints: Optional[FrozenSet[str]] = None,
) -> Dict[str, object]:
    """Matched-event transition audit with vs. without eukaryotic donor evidence (Astra 248976870271795200).

    Each tuple in `evaluations` is `(recipient, req_endpoints, prok_mags, euk_donors, unassigned_euk)`.
    Returns `{'n_evaluated': int, 'n_changed': int, 'transitions': {(prok_call, with_euk_call): count}}`.
    """
    transitions: Dict[Tuple[str, str], int] = {}
    n_eval = 0
    n_changed = 0
    for rec, reqs, prok_mags, euk_donors, unass_euk in evaluations:
        c_prok, _ = evaluate_recipient_requirements_in_event(
            rec, reqs, prok_mags, scenario=scenario
        )
        c_euk, _ = evaluate_recipient_requirements_in_event(
            rec,
            reqs,
            prok_mags,
            scenario=scenario,
            eukaryotic_donors=euk_donors,
            unassigned_euk_genome=unass_euk,
            euk_exportable_endpoints=euk_exportable_endpoints,
        )
        n_eval += 1
        if c_prok != c_euk:
            n_changed += 1
        transitions[(c_prok, c_euk)] = transitions.get((c_prok, c_euk), 0) + 1
    return {
        "n_evaluated": n_eval,
        "n_changed": n_changed,
        "transitions": transitions,
    }


# ---------------------------------------------------------------------------
# 5. Graded / Probabilistic Biological-State Layer (grid 248970390835478528)
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class GradedFanInPosterior:
    """Continuous biological-state summary per (recipient, endpoint_set, event, E).

    Computed over the latent biological capability state space Omega (NEVER by uniform
    weighting over enumerated routes, so adding redundant/duplicate wiring paths cannot
    distort probabilities; Astra 248970428366110720).
    """
    p_recipient_needs_all: float
    p_d0_ambient: float
    p_d1_single_supplier: float
    p_fanin_ge_2: float  # P(covered AND no single supplier suffices | X needs all, data, E)
    p_uncovered: float
    joint_p_plexitrophic: float  # p_recipient_needs_all * p_fanin_ge_2


def compute_pathway_posterior(
    rule: PathwayRule,
    g: GenomeAnnotations,
    prior_present: float = 0.5,
    fp_per_step: float = 0.02,
) -> float:
    """Smooth posterior P(pathway truly encoded in g | observed steps, contig_lengths, c, contam).

    Replaces brittle step-fraction cliffs with Bayes' rule over step observations
    weighted by intact-locus detectability and CheckM2 contamination.
    """
    ev = rule.evaluate(g, max_false_absence_prob=None)
    n = max(1, ev.n_steps)
    k = ev.n_present
    c = g.completeness / 100.0 if g.completeness > 1.0 else g.completeness
    c = min(max(float(c), 0.05), 0.999)
    contam = min(max(float(g.contamination) / 100.0, 0.0), 0.50)

    # Effective per-step detection sensitivity when truly present
    f_gene = p_intact(1200.0, g.contig_lengths, completeness=1.0) if g.contig_lengths else 0.85
    sens = min(max(c * (0.5 * f_gene + 0.5), 0.05), 0.995)
    # Effective per-step false-positive rate from annotation noise + bin contamination
    eps = min(max(fp_per_step + 0.5 * contam, 0.005), 0.45)

    log_lik_pres = k * math.log(sens) + (n - k) * math.log(1.0 - sens)
    if ev.key_ok and rule.key_indices:
        log_lik_pres += 0.25
    elif rule.key_indices and not ev.key_ok and k > 0:
        log_lik_pres -= 0.75
    if ev.end_truncated and k < n:
        # Contig-end truncation raises likelihood of missing steps under true presence
        log_lik_pres += 0.5 * (n - k) * math.log(max((1.0 - sens * 0.5) / (1.0 - sens), 1.0))

    log_lik_abs = k * math.log(eps) + (n - k) * math.log(1.0 - eps)
    p0 = min(max(prior_present, 1e-4), 1.0 - 1e-4)
    logit_post = math.log(p0 / (1.0 - p0)) + (log_lik_pres - log_lik_abs)
    logit_post = min(max(logit_post, -30.0), 30.0)
    return 1.0 / (1.0 + math.exp(-logit_post))


def compute_fanin_state_posterior(
    recipient_need_probs: Dict[str, float],
    donor_supply_probs: Sequence[Tuple[float, Dict[str, float]]] = (),
    ambient: Iterable[str] = (),
    donor_joint_pmfs: Optional[Sequence[Dict[FrozenSet[str], float]]] = None,
) -> GradedFanInPosterior:
    """Exact closed-form posterior over biological states (Proposition 1 marginalization).

    Supports BOTH:
      - `donor_supply_probs`: sequence of `(p_active_i, {ep: P(donor i exports ep | active)})`
        under conditional within-donor independence across `ep` given `active_i`; OR
      - `donor_joint_pmfs`: sequence of exact joint subset PMFs `{exported_subset: prob}`
        per donor `i != X` (sum of probs <= 1.0, residual assigned to `frozenset()`),
        handling arbitrary within-donor route coupling (e.g., mutual exclusion `p11 = 0`
        or shared precursor routes; Astra 248973173642936320).
    """
    if donor_supply_probs and donor_joint_pmfs is not None:
        raise ValueError("Pass either donor_supply_probs or donor_joint_pmfs, not both.")

    amb = set(ambient)
    reqs = list(recipient_need_probs.keys())
    p_need_all = 1.0
    for ep in reqs:
        p_need_all *= min(max(recipient_need_probs[ep], 0.0), 1.0)

    missing = [ep for ep in reqs if ep not in amb]
    if not missing:
        return GradedFanInPosterior(
            p_recipient_needs_all=round(p_need_all, 6),
            p_d0_ambient=1.0,
            p_d1_single_supplier=0.0,
            p_fanin_ge_2=0.0,
            p_uncovered=0.0,
            joint_p_plexitrophic=0.0,
        )

    missing_set = frozenset(missing)
    k = len(missing)

    if donor_joint_pmfs is not None:
        # Normalize/project each donor's joint PMF onto subsets of `missing_set`
        donor_az: List[Tuple[float, Dict[Tuple[str, ...], float]]] = []
        for pmf in donor_joint_pmfs:
            proj: Dict[FrozenSet[str], float] = {}
            tot = 0.0
            for s_exp, pr in pmf.items():
                p_clamped = max(0.0, float(pr))
                tot += p_clamped
                key = frozenset(s_exp) & missing_set
                proj[key] = proj.get(key, 0.0) + p_clamped
            if tot < 1.0:
                proj[frozenset()] = proj.get(frozenset(), 0.0) + (1.0 - tot)
            elif tot > 1.0 + 1e-9:
                for key in list(proj.keys()):
                    proj[key] /= tot
            a_i = min(max(proj.get(missing_set, 0.0), 0.0), 1.0)
            z_map: Dict[Tuple[str, ...], float] = {}
            for s_size in range(1, k + 1):
                for sub in itertools.combinations(missing, s_size):
                    sub_set = frozenset(sub)
                    z_map[sub] = min(
                        max(sum(pr for k_exp, pr in proj.items() if not (k_exp & sub_set)), 0.0),
                        1.0,
                    )
            donor_az.append((a_i, z_map))
    else:
        donor_az = []
        for p_act, cap_map in donor_supply_probs:
            pa = min(max(p_act, 0.0), 1.0)
            a_i = pa
            for ep in missing:
                a_i *= min(max(cap_map.get(ep, 0.0), 0.0), 1.0)
            z_map = {}
            for s_size in range(1, k + 1):
                for sub in itertools.combinations(missing, s_size):
                    p_none_given_act = 1.0
                    for ep in sub:
                        p_none_given_act *= 1.0 - min(max(cap_map.get(ep, 0.0), 0.0), 1.0)
                    z_map[sub] = (1.0 - pa) + pa * p_none_given_act
            donor_az.append((a_i, z_map))

    prod_not_all = 1.0
    for a_i, _ in donor_az:
        prod_not_all *= 1.0 - a_i

    p_uncov = 0.0
    for s_size in range(1, k + 1):
        for sub in itertools.combinations(missing, s_size):
            prod_none_sub = 1.0
            for _, z_map in donor_az:
                prod_none_sub *= z_map[sub]
            sign = 1.0 if (s_size % 2 == 1) else -1.0
            p_uncov += sign * prod_none_sub

    p_uncov = min(max(p_uncov, 0.0), 1.0)
    p_d1 = min(max(1.0 - prod_not_all, 0.0), 1.0)
    p_fanin = min(max(prod_not_all - p_uncov, 0.0), 1.0)

    return GradedFanInPosterior(
        p_recipient_needs_all=round(p_need_all, 6),
        p_d0_ambient=0.0,
        p_d1_single_supplier=round(p_d1, 6),
        p_fanin_ge_2=round(p_fanin, 6),
        p_uncovered=round(p_uncov, 6),
        joint_p_plexitrophic=round(p_need_all * p_fanin, 6),
    )


# ---------------------------------------------------------------------------
# 6. Contig-Block Downsampling (Null 1), Curveball, Bounds & Decision Bands
# ---------------------------------------------------------------------------


def downsample_contig_blocks(
    complete_hits: Sequence[GeneHit],
    target_completeness: float,
    target_contig_count: int,
    rng: random.Random,
    mag_id: str = "downsampled",
    domain: str = "Bacteria",
    orig_completeness: float = 100.0,
    total_genome_bp: int = 3_000_000,
) -> GenomeAnnotations:
    """Null 1: Contig-block fragmentation + completeness downsampling (Astra 248968746164994048 #3).

    Fixes all three defects in earlier draft:
      1. Splits ONLY within original contigs — never merges genes across two original contigs.
      2. Retains blocks uniformly across genome positions (shuffles BEFORE budget selection,
         never re-sorting indices prior to truncation).
      3. Computes and records the ACTUAL achieved completeness (`orig_completeness * kept / total`)
         and proportional simulated contig lengths in bp.
    """
    ordered = list(complete_hits)
    n_genes = len(ordered)
    if n_genes == 0:
        return GenomeAnnotations.from_kos(
            mag_id, (), domain=domain, completeness=0.0, contig_lengths=()
        )

    # Group genes by original contig preserving within-contig feature order
    orig_contigs: List[List[GeneHit]] = []
    cur_cid: Optional[str] = None
    for h in ordered:
        if h.contig_id != cur_cid:
            orig_contigs.append([h])
            cur_cid = h.contig_id
        else:
            orig_contigs[-1].append(h)

    # Allocate extra internal cuts across contigs that have >= 2 genes
    n_extra_cuts = max(0, min(target_contig_count - len(orig_contigs), n_genes - len(orig_contigs)))
    possible_cuts: List[Tuple[int, int]] = [
        (c_idx, g_pos)
        for c_idx, c_genes in enumerate(orig_contigs)
        for g_pos in range(1, len(c_genes))
    ]
    chosen_cuts = set(rng.sample(possible_cuts, n_extra_cuts)) if n_extra_cuts > 0 else set()

    blocks: List[List[GeneHit]] = []
    for c_idx, c_genes in enumerate(orig_contigs):
        cuts_here = sorted(g_pos for (ci, g_pos) in chosen_cuts if ci == c_idx)
        bounds = [0] + cuts_here + [len(c_genes)]
        for i in range(len(bounds) - 1):
            blocks.append(c_genes[bounds[i] : bounds[i + 1]])

    c_target_frac = target_completeness / 100.0 if target_completeness > 1.0 else target_completeness
    c_orig_frac = orig_completeness / 100.0 if orig_completeness > 1.0 else orig_completeness
    retain_ratio = min(max(c_target_frac / max(c_orig_frac, 1e-6), 0.0), 1.0)
    target_keep = max(1, int(round(retain_ratio * n_genes)))

    perm_indices = list(range(len(blocks)))
    rng.shuffle(perm_indices)

    selected_indices: List[int] = []
    kept_genes = 0
    for b_i in perm_indices:
        b_len = len(blocks[b_i])
        if not selected_indices:
            selected_indices.append(b_i)
            kept_genes += b_len
        elif abs((kept_genes + b_len) - target_keep) <= abs(kept_genes - target_keep):
            selected_indices.append(b_i)
            kept_genes += b_len
        if kept_genes >= target_keep:
            break

    bp_per_gene = max(500.0, float(total_genome_bp) / float(n_genes))
    kept_rows: List[Dict[str, str]] = []
    sim_contig_lengths: List[int] = []

    for out_idx, b_i in enumerate(sorted(selected_indices), start=1):
        blk = blocks[b_i]
        cid = f"simcontig_{out_idx:04d}"
        sim_contig_lengths.append(int(round(len(blk) * bp_per_gene)))
        for g_pos, h in enumerate(blk, start=1):
            kept_rows.append(
                {
                    "protein_id": f"{cid}_{g_pos:05d}",
                    "contig_id": cid,
                    "KO": ",".join(sorted(h.kos)),
                    "EC": ",".join(sorted(h.ecs)),
                }
            )

    orig_pct = orig_completeness if orig_completeness > 1.0 else orig_completeness * 100.0
    achieved_comp = round(orig_pct * (kept_genes / float(n_genes)), 4)

    genome, _ = parse_per_mag_rows(
        kept_rows,
        mag_id=mag_id,
        domain=domain,
        completeness=achieved_comp,
        genome_size_bp=sum(sim_contig_lengths),
        contig_lengths=sim_contig_lengths,
        complete_cds_table=True,
    )
    return genome


@dataclass
class StratumCounts:
    """Exhaustive 5-state partition of N_cand_aux (Section 4.1)."""
    n_auto: int = 0
    n_d0: int = 0
    n_d1: int = 0
    n_plexi: int = 0
    n_ambig: int = 0
    n_no_donor: int = 0

    @property
    def n_cand_aux(self) -> int:
        return self.n_d0 + self.n_d1 + self.n_plexi + self.n_ambig + self.n_no_donor

    @property
    def n_unresolved(self) -> int:
        return self.n_ambig + self.n_no_donor

    def record(self, call: str) -> None:
        if call == CALL_AUTONOMOUS:
            self.n_auto += 1
        elif call == CALL_D0_AMBIENT:
            self.n_d0 += 1
        elif call == CALL_D1_OR_CHAIN:
            self.n_d1 += 1
        elif call == CALL_PLEXITROPHIC:
            self.n_plexi += 1
        elif call == CALL_UNRESOLVED_AMBIG:
            self.n_ambig += 1
        elif call == CALL_UNRESOLVED_NO_DONOR:
            self.n_no_donor += 1

    def prevalence_bounds(self) -> Tuple[float, float]:
        n = self.n_cand_aux
        if n == 0:
            return (0.0, 1.0)
        return (self.n_plexi / n, (self.n_plexi + self.n_unresolved) / n)


def observed_minus_null_id_bounds(
    obs: StratumCounts,
    null: StratumCounts,
) -> Tuple[float, float]:
    a, b = obs.prevalence_bounds()
    c, d = null.prevalence_bounds()
    return (a - d, b - c)


def manski_log_odds_bounds(
    inc: StratumCounts,
    comp: StratumCounts,
    pseudocount: float = 0.5,
) -> Tuple[float, float]:
    n_i = inc.n_cand_aux
    n_c = comp.n_cand_aux
    if n_i == 0 or n_c == 0:
        return (-math.inf, math.inf)

    def _logit(k: int, n: int) -> float:
        if 0 < k < n:
            return math.log(k / (n - k))
        return math.log((k + pseudocount) / (n - k + pseudocount))

    low = _logit(inc.n_plexi, n_i) - _logit(comp.n_plexi + comp.n_unresolved, n_c)
    high = _logit(inc.n_plexi + inc.n_unresolved, n_i) - _logit(comp.n_plexi, n_c)
    return (low, high)


def event_bootstrap_q1_bounds(
    per_event_obs: Sequence[StratumCounts],
    per_event_null: Sequence[StratumCounts],
    n_boot: int = 1000,
    seed: int = 42,
    alpha: float = 0.05,
) -> Tuple[float, float]:
    if len(per_event_obs) != len(per_event_null) or not per_event_obs:
        return (-1.0, 1.0)
    rng = random.Random(seed)
    m = len(per_event_obs)
    lows: List[float] = []
    highs: List[float] = []
    for _ in range(n_boot):
        idxs = [rng.randrange(m) for _ in range(m)]
        obs_tot = StratumCounts(
            n_auto=sum(per_event_obs[i].n_auto for i in idxs),
            n_d0=sum(per_event_obs[i].n_d0 for i in idxs),
            n_d1=sum(per_event_obs[i].n_d1 for i in idxs),
            n_plexi=sum(per_event_obs[i].n_plexi for i in idxs),
            n_ambig=sum(per_event_obs[i].n_ambig for i in idxs),
            n_no_donor=sum(per_event_obs[i].n_no_donor for i in idxs),
        )
        nul_tot = StratumCounts(
            n_auto=sum(per_event_null[i].n_auto for i in idxs),
            n_d0=sum(per_event_null[i].n_d0 for i in idxs),
            n_d1=sum(per_event_null[i].n_d1 for i in idxs),
            n_plexi=sum(per_event_null[i].n_plexi for i in idxs),
            n_ambig=sum(per_event_null[i].n_ambig for i in idxs),
            n_no_donor=sum(per_event_null[i].n_no_donor for i in idxs),
        )
        lo, hi = observed_minus_null_id_bounds(obs_tot, nul_tot)
        lows.append(lo)
        highs.append(hi)
    lows.sort()
    highs.sort()
    q_lo = max(0, int(math.floor((alpha / 2.0) * n_boot)))
    q_hi = min(n_boot - 1, int(math.ceil((1.0 - alpha / 2.0) * n_boot)) - 1)
    return (lows[q_lo], highs[q_hi])


def curveball_permute_binary_matrix(
    matrix: Sequence[Sequence[int]],
    rng: random.Random,
    n_steps: int = 50,
    forbid_row_fn: Optional[object] = None,
) -> List[List[int]]:
    """Row- and column-sum-preserving Curveball (supports constrained N2c or GI incidence)."""
    n_rows = len(matrix)
    n_cols = len(matrix[0]) if n_rows else 0
    adj: List[Set[int]] = [{j for j, val in enumerate(row) if val} for row in matrix]
    if n_rows < 2:
        return [list(row) for row in matrix]

    for _ in range(n_steps * n_rows):
        i, k = rng.sample(range(n_rows), 2)
        si, sk = adj[i], adj[k]
        only_i = list(si - sk)
        only_k = list(sk - si)
        if not only_i or not only_k:
            continue
        pool = only_i + only_k
        rng.shuffle(pool)
        shared = si & sk
        ni = shared | set(pool[: len(only_i)])
        nk = shared | set(pool[len(only_i) :])
        if forbid_row_fn is not None and (forbid_row_fn(ni) or forbid_row_fn(nk)):  # type: ignore[operator]
            continue
        adj[i], adj[k] = ni, nk

    return [[1 if j in adj[i] else 0 for j in range(n_cols)] for i in range(n_rows)]


def check_data_adequacy_gates(
    n_independent_events: int,
    unbinned_donor_ko_share: float,
    unresolved_fraction: float,
    min_events: int = 10,
    max_unbinned_share: float = 0.50,
    max_unresolved_frac: float = 0.50,
) -> Optional[str]:
    if n_independent_events < min_events:
        return f"Inconclusive (data-adequacy gate: n_events={n_independent_events} < {min_events})"
    if unbinned_donor_ko_share > max_unbinned_share:
        return (
            f"Inconclusive (data-adequacy gate: unbinned_donor_ko_share="
            f"{unbinned_donor_ko_share:.3f} > {max_unbinned_share:.2f})"
        )
    if unresolved_fraction > max_unresolved_frac:
        return (
            f"Inconclusive (data-adequacy gate: unresolved_fraction="
            f"{unresolved_fraction:.3f} > {max_unresolved_frac:.2f})"
        )
    return None


def decide_p2_q1(
    ci_null1: Tuple[float, float],
    ci_null2: Tuple[float, float],
    delta_p: float = 0.05,
    gate_failure: Optional[str] = None,
) -> str:
    if gate_failure is not None:
        return gate_failure
    l1, u1 = ci_null1
    l2, u2 = ci_null2
    if u1 < delta_p or u2 < delta_p:
        return "Falsified"
    if l1 > delta_p and l2 > delta_p:
        return "Supported"
    return "Inconclusive"


def decide_p2_q2(
    ci_log_or: Tuple[float, float],
    delta_or: float = math.log(2.0),
    gate_failure: Optional[str] = None,
) -> str:
    if gate_failure is not None:
        return gate_failure
    lo, hi = ci_log_or
    if hi < delta_or:
        return "Falsified"
    if lo > delta_or:
        return "Supported"
    return "Inconclusive"
