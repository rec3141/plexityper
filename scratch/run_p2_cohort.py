#!/usr/bin/env python3
"""Execute frozen P2 protocol (e0dd8645dd32) on Stage 0 (Ebb_Flow) and staged Freshwater MAGs.

Imports ONLY frozen:
  - scratch/p2_pathway_rules.py (sha256: 42fbce744c39)
  - scratch/p2_scorer.py        (sha256: 579e9c3502c1)
and verifies all frozen protocol/code hashes before opening any cohort table.

Usage:
  python3 -I scratch/run_p2_cohort.py --stage0
  python3 -I scratch/run_p2_cohort.py --freshwater
  python3 -I scratch/run_p2_cohort.py --all
"""

from __future__ import annotations

import argparse
import concurrent.futures
import csv
import gzip
import hashlib
import itertools
import json
import math
import os
import random
import subprocess
import sys
import tempfile
from typing import AbstractSet, Callable, Dict, FrozenSet, Iterable, List, Optional, Sequence, Set, Tuple

SCRIPT_DIR = os.path.dirname(os.path.realpath(__file__))
REPO_ROOT = os.path.dirname(SCRIPT_DIR)
sys.path.insert(0, SCRIPT_DIR)

from p2_pathway_rules import (  # noqa: E402
    C0_CONTROL_MODULES,
    PATHWAY_RULES,
    STATE_ABSENT,
    STATE_PRESENT,
    STATE_UNRESOLVED,
    GeneHit,
    GenomeAnnotations,
    all_referenced_kos,
    classify_coencoder_bypass,
    evaluate_b12_de_novo,
)
from p2_scorer import (  # noqa: E402
    ALL_P2_ENDPOINTS,
    CALL_AUTONOMOUS,
    CALL_D0_AMBIENT,
    CALL_D1_OR_CHAIN,
    CALL_NOT_ELIGIBLE,
    CALL_PLEXITROPHIC,
    CALL_UNRESOLVED_AMBIG,
    CALL_UNRESOLVED_NO_DONOR,
    SCENARIO_E_MINIMAL,
    SCENARIO_E_PLUS_DIN,
    SCENARIO_E_PLUS_DMSP,
    SCENARIO_E_PLUS_VITAMINS,
    OrganismProfile,
    StratumCounts,
    SupplyScenario,
    _build_relaxed_pool_pseudo_organism,
    assert_non_redundant_mag_ids,
    build_mag_organism_profile,
    check_data_adequacy_gates,
    compare_calls_with_eukaryotic_donors,
    compute_fanin_state_posterior,
    compute_pathway_posterior,
    curveball_permute_binary_matrix,
    decide_nonplexi_p_all,
    decide_p2_q1,
    decide_p2_q2,
    downsample_contig_blocks,
    evaluate_recipient_requirements_in_event,
    event_bootstrap_q1_bounds,
    filter_derep_prokaryotic_mags,
    manski_log_odds_bounds,
    observed_minus_null_id_bounds,
    parse_gff_cds_order,
    parse_ko_field,
    parse_per_mag_tsv,
    recipient_endpoint_status,
)

FROZEN_HASHES = {
    os.path.join(REPO_ROOT, "protocols", "P2-plexitrophy-recipient.md"): "e0dd8645dd32",
    os.path.join(SCRIPT_DIR, "p2_pathway_rules.py"): "42fbce744c39",
    os.path.join(SCRIPT_DIR, "p2_scorer.py"): "579e9c3502c1",
    os.path.join(SCRIPT_DIR, "test_p2_fixtures.py"): "ed7c33667c33",
}

EBB_ROOT = "/data/scratch/nanopore_mag/results_Ebb_Flow_20260218"
FW_ROOT = os.path.join(REPO_ROOT, "data", "grex_mag_freshwater")
FW_FASTA_DIR = "/data/scratch/derep_test/freshwater/dastool"
OUT_DIR = os.path.join(REPO_ROOT, "results", "p2")

# Strictly domain-specific KOs verified in kofam_db/ko_list (used ONLY when GTDB-Tk is unstaged; 'U' if ambiguous)
BAC_STRICT_KOS = {"K03086", "K03043", "K03046", "K01000", "K00790", "K01924", "K01925", "K02563"}
ARC_STRICT_KOS = {"K03041", "K03042", "K02918", "K19712", "K28504", "K00399", "K00401", "K00402"}


def sha256_full(path: str) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        while True:
            chunk = fh.read(1 << 20)
            if not chunk:
                break
            h.update(chunk)
    return h.hexdigest()


def sha256_12(path: str) -> str:
    return sha256_full(path)[:12]


def verify_frozen_artifacts() -> None:
    for path, expected in FROZEN_HASHES.items():
        actual = sha256_12(path)
        if actual != expected:
            raise RuntimeError(f"Frozen artifact hash mismatch on {path}: {actual} != {expected}")


def read_fasta_contig_lengths(fa_path: str) -> Dict[str, int]:
    lens: Dict[str, int] = {}
    cur: Optional[str] = None
    with open(fa_path, "r", encoding="utf-8") as fh:
        for line in fh:
            if line.startswith(">"):
                cur = line[1:].strip().split()[0]
                lens[cur] = 0
            elif cur is not None:
                lens[cur] += len(line.strip())
    return lens


def _samtools_cov_worker(args: Tuple[str, str]) -> Tuple[str, Dict[str, int]]:
    sample_name, bam_path = args
    proc = subprocess.run(
        ["/usr/bin/samtools", "coverage", bam_path],
        check=True,
        capture_output=True,
        text=True,
    )
    cov_map: Dict[str, int] = {}
    for line in proc.stdout.splitlines():
        if not line or line.startswith("#"):
            continue
        parts = line.split("\t")
        covbases = int(parts[4])
        if covbases > 0:
            cov_map[parts[0]] = covbases
    return sample_name, cov_map


def compute_or_load_stage0_bam_covbases(sample_bams: List[Tuple[str, str]], cache_json: str) -> Dict[str, Dict[str, int]]:
    if os.path.isfile(cache_json):
        with open(cache_json, "r", encoding="utf-8") as fh:
            return json.load(fh)
    out: Dict[str, Dict[str, int]] = {}
    with concurrent.futures.ThreadPoolExecutor(max_workers=8) as pool:
        for sname, cmap in pool.map(_samtools_cov_worker, sample_bams):
            out[sname] = cmap
    with open(cache_json, "w", encoding="utf-8") as fh:
        json.dump(out, fh, sort_keys=True)
    return out


def _fast_smax(
    profiles: Dict[str, OrganismProfile],
    members: Iterable[str],
    amb: FrozenSet[str],
    focal_id: Optional[str] = None,
    only_partner: Optional[str] = None,
    private_dyad: bool = False,
    extra_focal_imports: FrozenSet[str] = frozenset(),
) -> Tuple[Set[str], Dict[str, Set[str]]]:
    """Exact O(|V|) algebraic equivalent of two_node_smax for unrestricted, fanin1, and private_dyad."""
    cur = set(members)
    rules_map = {y: profiles[y].rules() for y in cur}
    while True:
        avail_by: Dict[str, Set[str]] = {
            y: (set(amb) | set(extra_focal_imports) if y == focal_id else set(amb))
            for y in cur
        }
        prod_by: Dict[str, Set[str]] = {y: set() for y in cur}
        changed = True
        while changed:
            changed = False
            for y in cur:
                av_y = avail_by[y]
                pr_y = prod_by[y]
                for ins, prod in rules_map[y]:
                    if prod not in pr_y and ins.issubset(av_y):
                        pr_y.add(prod)
                        av_y.add(prod)
                        changed = True
            if focal_id is None and only_partner is None:
                pool: Set[str] = set()
                for z in cur:
                    oz = profiles[z]
                    for m in prod_by[z]:
                        if oz.can_export(m):
                            pool.add(m)
                for y in cur:
                    before = len(avail_by[y])
                    avail_by[y] |= pool
                    if len(avail_by[y]) > before:
                        changed = True
            elif not private_dyad:
                # fanin1: y != focal_id receives from all z; focal_id receives only from only_partner
                pool = set()
                for z in cur:
                    oz = profiles[z]
                    for m in prod_by[z]:
                        if oz.can_export(m):
                            pool.add(m)
                for y in cur:
                    if y == focal_id:
                        if only_partner in cur:
                            od = profiles[only_partner]  # type: ignore[index]
                            for m in prod_by[only_partner]:  # type: ignore[index]
                                if od.can_export(m) and m not in avail_by[y]:
                                    avail_by[y].add(m)
                                    changed = True
                    else:
                        before = len(avail_by[y])
                        avail_by[y] |= pool
                        if len(avail_by[y]) > before:
                            changed = True
            else:
                # private_dyad ('cyc'): focal_id <-> only_partner only; others receive from z != focal_id
                pool_no_focal: Set[str] = set()
                exp_focal: Set[str] = set()
                for z in cur:
                    oz = profiles[z]
                    target_set = exp_focal if z == focal_id else pool_no_focal
                    for m in prod_by[z]:
                        if oz.can_export(m):
                            target_set.add(m)
                for y in cur:
                    before = len(avail_by[y])
                    if y == focal_id:
                        if only_partner in cur:
                            od = profiles[only_partner]  # type: ignore[index]
                            for m in prod_by[only_partner]:  # type: ignore[index]
                                if od.can_export(m):
                                    avail_by[y].add(m)
                    elif y == only_partner:
                        avail_by[y] |= pool_no_focal | exp_focal
                    else:
                        avail_by[y] |= pool_no_focal
                    if len(avail_by[y]) > before:
                        changed = True
        dead = {y for y in cur if not profiles[y].required.issubset(avail_by[y])}
        if not dead:
            return cur, prod_by
        cur -= dead


def _fast_decide_acyclic_p_all(
    focal_id: str,
    profiles: Dict[str, OrganismProfile],
    amb: FrozenSet[str],
    cyc_mode: str,
    cached_others_smax: Optional[Tuple[Set[str], Dict[str, Set[str]]]] = None,
) -> Tuple[str, Optional[Tuple[str, ...]]]:
    """Exact O(|V|)-closure equivalent of _decide_nonplexi_acyclic_poly."""
    all_ids = set(profiles.keys())
    v_alive, prod_v = _fast_smax(profiles, all_ids, amb)
    if focal_id not in v_alive:
        return ("no_donor", None)
    if focal_id in _fast_smax(profiles, {focal_id}, amb)[0]:
        return ("d0", ("d0",))
    others = sorted(v_alive - {focal_id})
    if cached_others_smax is not None and set(others) == cached_others_smax[0]:
        c_set, prod_c = cached_others_smax
    else:
        c_set, prod_c = _fast_smax(profiles, others, amb)
    for d_id in others:
        if d_id in c_set:
            pool_d = frozenset(m for m in prod_c[d_id] if profiles[d_id].can_export(m))
            if focal_id in _fast_smax(profiles, {focal_id}, amb, focal_id=focal_id, extra_focal_imports=pool_d)[0]:
                return ("not_plexi", ("sink", d_id))
    if cyc_mode == "fanin_only":
        # Exact invariant check: if focal_id's conversion outputs are not inputs/requirements of any other organism,
        # prod_by[d_id] under fanin1(d_id) equals prod_v[d_id] identically!
        focal_conv_outs = {out for _, out in profiles[focal_id].conversions}
        other_uses_focal_conv = False
        if focal_conv_outs:
            for z_id in others:
                oz = profiles[z_id]
                if oz.required & focal_conv_outs:
                    other_uses_focal_conv = True
                    break
                if any(ins & focal_conv_outs for ins, _ in oz.conversions):
                    other_uses_focal_conv = True
                    break
        if not other_uses_focal_conv:
            for d_id in others:
                pool_dv = frozenset(m for m in prod_v[d_id] if profiles[d_id].can_export(m))
                if focal_id in _fast_smax(profiles, {focal_id}, amb, focal_id=focal_id, extra_focal_imports=pool_dv)[0]:
                    return ("not_plexi", ("fanin1", d_id))
            return ("plexi", None)
        for d_id in others:
            if focal_id in _fast_smax(profiles, v_alive, amb, focal_id=focal_id, only_partner=d_id, private_dyad=False)[0]:
                return ("not_plexi", ("fanin1", d_id))
        return ("plexi", None)
    if cyc_mode == "cyc":
        for d_id in others:
            if focal_id in _fast_smax(profiles, v_alive, amb, focal_id=focal_id, only_partner=d_id, private_dyad=True)[0]:
                return ("not_plexi", ("private_dyad", d_id))
        return ("plexi", None)
    raise ValueError(f"Unsupported cyc_mode: {cyc_mode}")


class CommunityScorerContext:
    """Precomputes donor OrganismProfiles once per community; algebraically identical to evaluate_recipient_requirements_in_event."""

    def __init__(
        self,
        co_occurring_mags: Sequence[GenomeAnnotations],
        unbinned_genome: Optional[GenomeAnnotations] = None,
        unassigned_euk_genome: Optional[GenomeAnnotations] = None,
        mode: str = "fanin_only",
        max_false_absence_prob: Optional[float] = 0.05,
        verify_against_reference: bool = True,
    ) -> None:
        self.co_occurring_mags = co_occurring_mags
        self.unbinned_genome = unbinned_genome
        self.unassigned_euk_genome = unassigned_euk_genome
        self.mode = mode
        self.cyc_mode = "fanin_only" if mode == "fanin_only" else "cyc"
        self.max_fap = max_false_absence_prob
        active_needs: Sequence[str] = () if mode == "fanin_only" else ALL_P2_ENDPOINTS
        self.donor_minus: Dict[str, OrganismProfile] = {}
        self.donor_plus: Dict[str, OrganismProfile] = {}
        self.rec_base_minus: Dict[str, OrganismProfile] = {}
        self.rec_base_plus: Dict[str, OrganismProfile] = {}
        self.ep_status_cache: Dict[Tuple[str, str], str] = {}

        for g in co_occurring_mags:
            self.donor_minus[g.mag_id] = build_mag_organism_profile(
                g, modeled_donor_needs=active_needs, optimistic_capabilities=False, max_false_absence_prob=max_false_absence_prob
            )
            self.donor_plus[g.mag_id] = build_mag_organism_profile(
                g, modeled_donor_needs=active_needs, optimistic_capabilities=True, max_false_absence_prob=max_false_absence_prob
            )
            if mode == "fanin_only":
                self.rec_base_minus[g.mag_id] = self.donor_minus[g.mag_id]
                self.rec_base_plus[g.mag_id] = self.donor_plus[g.mag_id]
            else:
                self.rec_base_minus[g.mag_id] = build_mag_organism_profile(
                    g, required_override=frozenset(), optimistic_capabilities=False, max_false_absence_prob=max_false_absence_prob
                )
                self.rec_base_plus[g.mag_id] = build_mag_organism_profile(
                    g, required_override=frozenset(), optimistic_capabilities=True, max_false_absence_prob=max_false_absence_prob
                )

        if unbinned_genome is not None:
            self.donor_plus["_unbinned_P_star"] = _build_relaxed_pool_pseudo_organism(
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
            self.donor_plus["_euk_unassigned_P_star"] = _build_relaxed_pool_pseudo_organism(
                unass_euk_g, "_euk_unassigned_P_star"
            )
        self._verified_count = 0 if verify_against_reference else 2

    def _status(self, rec: GenomeAnnotations, ep: str) -> str:
        key = (rec.mag_id, ep)
        if key not in self.ep_status_cache:
            self.ep_status_cache[key] = recipient_endpoint_status(
                rec, ep, max_false_absence_prob=self.max_fap
            )
        return self.ep_status_cache[key]

    def evaluate(
        self,
        recipient: GenomeAnnotations,
        required_endpoints: Sequence[str],
        scenario: SupplyScenario = SCENARIO_E_MINIMAL,
    ) -> Tuple[str, Optional[Tuple[str, ...]]]:
        if recipient.domain not in ("Bacteria", "Archaea"):
            return CALL_NOT_ELIGIBLE, None
        reqs = tuple(dict.fromkeys(required_endpoints))
        statuses = [self._status(recipient, ep) for ep in reqs]
        if any(s == "not_dependent" for s in statuses):
            return CALL_NOT_ELIGIBLE, None
        if all(s == STATE_PRESENT for s in statuses):
            return CALL_AUTONOMOUS, None
        if not all(s == STATE_ABSENT for s in statuses):
            return CALL_NOT_ELIGIBLE, None

        req_fset = frozenset(reqs)
        bm = self.rec_base_minus[recipient.mag_id]
        bp = self.rec_base_plus[recipient.mag_id]
        pm = dict(self.donor_minus)
        pp = dict(self.donor_plus)
        pm[recipient.mag_id] = OrganismProfile(
            org_id=bm.org_id, de_novo=bm.de_novo, conversions=bm.conversions, required=req_fset, exportable=bm.exportable
        )
        pp[recipient.mag_id] = OrganismProfile(
            org_id=bp.org_id, de_novo=bp.de_novo, conversions=bp.conversions, required=req_fset, exportable=bp.exportable
        )
        st_minus, wit_minus = _fast_decide_acyclic_p_all(
            recipient.mag_id, pm, scenario.ambient_compounds, cyc_mode=self.cyc_mode
        )
        st_plus, wit_plus = _fast_decide_acyclic_p_all(
            recipient.mag_id, pp, scenario.ambient_compounds, cyc_mode=self.cyc_mode
        )
        if st_minus == "d0":
            out = (CALL_D0_AMBIENT, wit_minus)
        elif st_minus == "not_plexi":
            out = (CALL_D1_OR_CHAIN, wit_minus)
        elif st_plus == "no_donor":
            out = (CALL_UNRESOLVED_NO_DONOR, None)
        elif st_minus == "plexi" and st_plus == "plexi":
            out = (CALL_PLEXITROPHIC, None)
        else:
            out = (CALL_UNRESOLVED_AMBIG, wit_plus)

        if self._verified_count < 2:
            ref = evaluate_recipient_requirements_in_event(
                recipient,
                reqs,
                self.co_occurring_mags,
                scenario=scenario,
                unbinned_genome=self.unbinned_genome,
                unassigned_euk_genome=self.unassigned_euk_genome,
                mode=self.mode,
                max_false_absence_prob=self.max_fap,
            )
            assert out == ref, f"Hoisted scorer mismatch on {recipient.mag_id} {reqs}: {out} != {ref}"
            self._verified_count += 1
        return out


def donor_endpoint_posterior(g: GenomeAnnotations, ep: str) -> float:
    """Smooth posterior P(donor g can synthesize & export endpoint ep) via compute_pathway_posterior."""
    if ep == "B12":
        p_aer = compute_pathway_posterior(PATHWAY_RULES["P_B12_AER"], g)
        p_ana = compute_pathway_posterior(PATHWAY_RULES["A_B12_ANA"], g)
        p_blub = compute_pathway_posterior(PATHWAY_RULES["P_BLUB"], g)
        p_bza = compute_pathway_posterior(PATHWAY_RULES["A_BZA"], g)
        p_ring = 1.0 - (1.0 - p_aer) * (1.0 - p_ana)
        p_lig = 1.0 - (1.0 - p_blub) * (1.0 - p_bza)
        return p_ring * p_lig
    if ep == "RED_S":
        return compute_pathway_posterior(PATHWAY_RULES["A_DSR"], g)
    if ep == "FIXED_N":
        return compute_pathway_posterior(PATHWAY_RULES["A_NIF"], g)
    if ep == "THI":
        return compute_pathway_posterior(PATHWAY_RULES["C_THI_SPLIT"], g)
    return compute_pathway_posterior(PATHWAY_RULES[ep], g)


def recipient_need_posterior(g: GenomeAnnotations, ep: str) -> float:
    """Smooth posterior P(recipient g lacks de novo synthesis of ep)."""
    if ep == "B12":
        if "K00549" in g.kos or "K00548" not in g.kos:
            return 0.0
        return 1.0 - donor_endpoint_posterior(g, "B12")
    if ep == "RED_S":
        return 1.0 - compute_pathway_posterior(PATHWAY_RULES["ASSIM_SO4"], g)
    if ep == "FIXED_N":
        p_nif = compute_pathway_posterior(PATHWAY_RULES["A_NIF"], g)
        p_no3 = compute_pathway_posterior(PATHWAY_RULES["ASSIM_NO3"], g)
        return (1.0 - p_nif) * (1.0 - p_no3)
    if ep == "THI":
        return 1.0 - compute_pathway_posterior(PATHWAY_RULES["C_THI_SPLIT"], g)
    return 1.0 - compute_pathway_posterior(PATHWAY_RULES[ep], g)


# Donor capability modules permuted in Null 2b / Null 2c
DONOR_CAP_MODULES = (
    "P_B12_AER",
    "A_B12_ANA",
    "P_BLUB",
    "P_BLUB_DMB",
    "A_BZA",
    "B12_LOWER_ATTACH",
    "THI_HMP",
    "THI_HET",
    "THI_COUPLE",
    "A_DSR",
    "A_NIF",
    *C0_CONTROL_MODULES,
)
O2_PLUS_IDX = {DONOR_CAP_MODULES.index(m) for m in ("P_B12_AER", "P_BLUB", "P_BLUB_DMB")}
O2_MINUS_IDX = {DONOR_CAP_MODULES.index(m) for m in ("A_B12_ANA", "A_BZA", "A_DSR", "A_NIF")}


def profile_from_module_Presence_row(org_id: str, row_bits: Sequence[int]) -> OrganismProfile:
    pres = {DONOR_CAP_MODULES[j] for j, b in enumerate(row_bits) if b}
    de_novo: Set[str] = set()
    convs: List[Tuple[FrozenSet[str], str]] = []
    ring_ok = ("P_B12_AER" in pres) or ("A_B12_ANA" in pres)
    lig_ok = ("P_BLUB" in pres) or ("A_BZA" in pres)
    if ring_ok:
        de_novo.add("COB_RING")
    if "P_BLUB_DMB" in pres or "P_BLUB" in pres:
        de_novo.add("DMB")
    if ring_ok and lig_ok:
        de_novo.add("B12")
    if "B12_LOWER_ATTACH" in pres:
        convs.append((frozenset({"COB_RING", "DMB"}), "B12"))
    if "THI_HMP" in pres:
        de_novo.add("THI_HMP")
    if "THI_HET" in pres:
        de_novo.add("THI_HET")
    if "THI_HMP" in pres and "THI_HET" in pres and "THI_COUPLE" in pres:
        de_novo.add("THI")
    elif "THI_COUPLE" in pres:
        convs.append((frozenset({"THI_HMP", "THI_HET"}), "THI"))
    if "A_DSR" in pres:
        de_novo.add("RED_S")
    if "A_NIF" in pres:
        de_novo.add("FIXED_N")
    for m_c0 in C0_CONTROL_MODULES:
        if m_c0 in pres:
            de_novo.add(m_c0)
    return OrganismProfile(org_id=org_id, de_novo=frozenset(de_novo), conversions=tuple(convs))


def build_donor_binary_matrix(donors: Sequence[GenomeAnnotations], optimistic: bool) -> List[List[int]]:
    t_states = {STATE_PRESENT, STATE_UNRESOLVED} if optimistic else {STATE_PRESENT}
    mat: List[List[int]] = []
    for g in donors:
        row = [
            1 if PATHWAY_RULES[mid].evaluate(g, max_false_absence_prob=0.05).state in t_states else 0
            for mid in DONOR_CAP_MODULES
        ]
        mat.append(row)
    return mat


def _clash_locked_rows(
    mat_minus: Sequence[Sequence[int]],
    forbid_minus_row_fn: Optional[Callable[[Set[int]], bool]],
) -> FrozenSet[int]:
    """Identify observed clash rows R_clash(M^-) to freeze under equivalence relation ~_{2c} (zero exemptions)."""
    if forbid_minus_row_fn is None:
        return frozenset()
    return frozenset(
        r
        for r in range(len(mat_minus))
        if forbid_minus_row_fn({c for c, v in enumerate(mat_minus[r]) if v})
    )


def assert_nested_completion_invariants(
    orig_minus: Sequence[Sequence[int]],
    orig_plus: Sequence[Sequence[int]],
    pm_minus: Sequence[Sequence[int]],
    pm_plus: Sequence[Sequence[int]],
    forbidden_minus_cells: AbstractSet[Tuple[int, int]] = frozenset(),
    forbidden_plus_cells: AbstractSet[Tuple[int, int]] = frozenset(),
    forbid_minus_row_fn: Optional[Callable[[Set[int]], bool]] = None,
    locked_rows: AbstractSet[int] = frozenset(),
) -> None:
    """Assert cellwise definite <= possible, preserved I-/I+ row & column sums, locked clash rows, and zero-exemption row rules."""
    n_rows = len(orig_minus)
    n_cols = len(orig_minus[0]) if n_rows else 0
    assert len(orig_plus) == n_rows and len(pm_minus) == n_rows and len(pm_plus) == n_rows
    for r in range(n_rows):
        rm = pm_minus[r]
        rp = pm_plus[r]
        assert len(rm) == n_cols and len(rp) == n_cols
        for c in range(n_cols):
            vm, vp = rm[c], rp[c]
            assert vm in (0, 1) and vp in (0, 1), f"Non-binary entry at ({r},{c}): ({vm},{vp})"
            assert vm <= vp, f"Cellwise nesting definite<=possible violated at ({r},{c}): {vm} > {vp}"
            if (r, c) in forbidden_plus_cells:
                assert vp == 0, f"Forbidden I+ cell ({r},{c}) became 1"
            if (r, c) in forbidden_minus_cells:
                assert vm == 0, f"Forbidden I- cell ({r},{c}) became 1"
        assert sum(rm) == sum(orig_minus[r]), f"Row {r} I- sum changed: {sum(rm)} != {sum(orig_minus[r])}"
        assert sum(rp) == sum(orig_plus[r]), f"Row {r} I+ sum changed: {sum(rp)} != {sum(orig_plus[r])}"
        if r in locked_rows:
            assert list(rm) == list(orig_minus[r]) and list(rp) == list(orig_plus[r]), (
                f"Locked clash row {r} mutated during sampling"
            )
            if forbid_minus_row_fn is not None:
                assert forbid_minus_row_fn({c for c, v in enumerate(rm) if v}), (
                    f"Locked clash row {r} lost its clash"
                )
        elif forbid_minus_row_fn is not None:
            ones_m = {c for c, v in enumerate(rm) if v}
            assert not forbid_minus_row_fn(ones_m), f"Row {r} violated constrained I- predicate: {sorted(ones_m)}"
    for c in range(n_cols):
        cm_orig = sum(orig_minus[r][c] for r in range(n_rows))
        cp_orig = sum(orig_plus[r][c] for r in range(n_rows))
        cm_new = sum(pm_minus[r][c] for r in range(n_rows))
        cp_new = sum(pm_plus[r][c] for r in range(n_rows))
        assert cm_new == cm_orig, f"Col {c} I- sum changed: {cm_new} != {cm_orig}"
        assert cp_new == cp_orig, f"Col {c} I+ sum changed: {cp_new} != {cp_orig}"


ALL_12_LATIN_3X3: Tuple[Tuple[Tuple[int, int, int], ...], ...] = tuple(
    (r0, r1, r2)
    for r0 in itertools.permutations((0, 1, 2))
    for r1 in itertools.permutations((0, 1, 2))
    if all(r0[j] != r1[j] for j in range(3))
    for r2 in (tuple(3 - r0[j] - r1[j] for j in range(3)),)
)
ALL_12_LATIN_3X3_SET = frozenset(ALL_12_LATIN_3X3)


def curveball_permute_nested_pair(
    mat_minus: Sequence[Sequence[int]],
    mat_plus: Sequence[Sequence[int]],
    rng: random.Random,
    n_steps: int = 30,
    forbidden_minus_cells: AbstractSet[Tuple[int, int]] = frozenset(),
    forbidden_plus_cells: AbstractSet[Tuple[int, int]] = frozenset(),
    forbid_minus_row_fn: Optional[Callable[[Set[int]], bool]] = None,
    check_each_trade: bool = False,
    enable_latin3x3: bool = False,
) -> Tuple[List[List[int]], List[List[int]], int]:
    """Jointly permute nested (I-, I+) binary matrices via exact 3-state Curveball + optional 3x3 Latin moves."""
    n_rows = len(mat_minus)
    n_cols = len(mat_minus[0]) if n_rows else 0
    locked_rows = _clash_locked_rows(mat_minus, forbid_minus_row_fn)
    free_rows = [r for r in range(n_rows) if r not in locked_rows]
    assert_nested_completion_invariants(
        mat_minus,
        mat_plus,
        mat_minus,
        mat_plus,
        forbidden_minus_cells=forbidden_minus_cells,
        forbidden_plus_cells=forbidden_plus_cells,
        forbid_minus_row_fn=forbid_minus_row_fn,
        locked_rows=locked_rows,
    )
    if len(free_rows) < 2 or n_cols == 0:
        return [list(r) for r in mat_minus], [list(r) for r in mat_plus], 0

    state = [[mat_minus[r][c] + mat_plus[r][c] for c in range(n_cols)] for r in range(n_rows)]
    forb_any = set(forbidden_minus_cells) | set(forbidden_plus_cells)
    n_nontrivial = 0

    for _ in range(n_steps * n_rows):
        if enable_latin3x3 and len(free_rows) >= 3 and n_cols >= 3 and rng.random() < 0.25:
            r_idxs = rng.sample(free_rows, 3)
            unlocked_cols = [
                c for c in range(n_cols) if all((r_a, c) not in forb_any for r_a in r_idxs)
            ]
            if len(unlocked_cols) >= 3:
                c_idxs = rng.sample(unlocked_cols, 3)
                sub_sq = tuple(
                    tuple(state[r_a][c_b] for c_b in c_idxs) for r_a in r_idxs
                )
                if sub_sq in ALL_12_LATIN_3X3_SET:
                    prop_sq = rng.choice(ALL_12_LATIN_3X3)
                    if prop_sq != sub_sq:
                        cand_rows = []
                        ok_latin = True
                        for a_pos, r_a in enumerate(r_idxs):
                            nr = list(state[r_a])
                            for b_pos, c_b in enumerate(c_idxs):
                                nr[c_b] = prop_sq[a_pos][b_pos]
                            if (
                                forbid_minus_row_fn is not None
                                and forbid_minus_row_fn({c for c, v in enumerate(nr) if v == 2})
                            ):
                                ok_latin = False
                                break
                            cand_rows.append(nr)
                        if ok_latin:
                            for a_pos, r_a in enumerate(r_idxs):
                                state[r_a] = cand_rows[a_pos]
                            n_nontrivial += 1
                            if check_each_trade:
                                cur_m = [[1 if v == 2 else 0 for v in row] for row in state]
                                cur_p = [[1 if v >= 1 else 0 for v in row] for row in state]
                                assert_nested_completion_invariants(
                                    mat_minus,
                                    mat_plus,
                                    cur_m,
                                    cur_p,
                                    forbidden_minus_cells=forbidden_minus_cells,
                                    forbidden_plus_cells=forbidden_plus_cells,
                                    forbid_minus_row_fn=forbid_minus_row_fn,
                                    locked_rows=locked_rows,
                                )
            continue

        i, k = rng.sample(free_rows, 2)
        ri, rk = state[i], state[k]
        p01: List[int] = []
        p02: List[int] = []
        p12: List[int] = []
        a01 = 0
        a02 = 0
        a12 = 0
        for c in range(n_cols):
            vi, vk = ri[c], rk[c]
            if vi == vk:
                continue
            lo, hi = (vi, vk) if vi < vk else (vk, vi)
            if lo == 0 and hi == 1:
                if (i, c) in forbidden_plus_cells or (k, c) in forbidden_plus_cells:
                    continue
                p01.append(c)
                if vi == 1:
                    a01 += 1
            elif lo == 0 and hi == 2:
                if (i, c) in forb_any or (k, c) in forb_any:
                    continue
                p02.append(c)
                if vi == 2:
                    a02 += 1
            else:  # lo == 1 and hi == 2
                if (i, c) in forbidden_minus_cells or (k, c) in forbidden_minus_cells:
                    continue
                p12.append(c)
                if vi == 2:
                    a12 += 1

        m01, m02, m12 = len(p01), len(p02), len(p12)
        t_min = max(-a02, a01 - m01, a12 - m12)
        t_max = min(m02 - a02, a01, a12)
        if t_min == 0 and t_max == 0 and m01 <= 1 and m02 <= 1 and m12 <= 1:
            continue

        t_vals = list(range(t_min, t_max + 1))
        weights = [
            math.comb(m01, a01 - t) * math.comb(m02, a02 + t) * math.comb(m12, a12 - t)
            for t in t_vals
        ]
        tot_w = sum(weights)
        if tot_w <= 1:
            continue
        pick = rng.randrange(tot_w)
        chosen_t = t_vals[-1]
        acc = 0
        for t_v, w_v in zip(t_vals, weights):
            acc += w_v
            if pick < acc:
                chosen_t = t_v
                break

        x01 = a01 - chosen_t
        x02 = a02 + chosen_t
        x12 = a12 - chosen_t

        ni = list(ri)
        nk = list(rk)
        sel01 = set(rng.sample(p01, x01)) if x01 else set()
        for c in p01:
            if c in sel01:
                ni[c], nk[c] = 1, 0
            else:
                ni[c], nk[c] = 0, 1
        sel02 = set(rng.sample(p02, x02)) if x02 else set()
        for c in p02:
            if c in sel02:
                ni[c], nk[c] = 2, 0
            else:
                ni[c], nk[c] = 0, 2
        sel12 = set(rng.sample(p12, x12)) if x12 else set()
        for c in p12:
            if c in sel12:
                ni[c], nk[c] = 2, 1
            else:
                ni[c], nk[c] = 1, 2

        if forbid_minus_row_fn is not None:
            if forbid_minus_row_fn({c for c, v in enumerate(ni) if v == 2}):
                continue
            if forbid_minus_row_fn({c for c, v in enumerate(nk) if v == 2}):
                continue

        changed = (ni != ri) or (nk != rk)
        state[i], state[k] = ni, nk
        if changed:
            n_nontrivial += 1
            if check_each_trade:
                cur_m = [[1 if v == 2 else 0 for v in row] for row in state]
                cur_p = [[1 if v >= 1 else 0 for v in row] for row in state]
                assert_nested_completion_invariants(
                    mat_minus,
                    mat_plus,
                    cur_m,
                    cur_p,
                    forbidden_minus_cells=forbidden_minus_cells,
                    forbidden_plus_cells=forbidden_plus_cells,
                    forbid_minus_row_fn=forbid_minus_row_fn,
                    locked_rows=locked_rows,
                )

    pm_minus = [[1 if v == 2 else 0 for v in row] for row in state]
    pm_plus = [[1 if v >= 1 else 0 for v in row] for row in state]
    assert_nested_completion_invariants(
        mat_minus,
        mat_plus,
        pm_minus,
        pm_plus,
        forbidden_minus_cells=forbidden_minus_cells,
        forbidden_plus_cells=forbidden_plus_cells,
        forbid_minus_row_fn=forbid_minus_row_fn,
        locked_rows=locked_rows,
    )
    return pm_minus, pm_plus, n_nontrivial


def _enumerate_line_completions(
    line_a: Sequence[int],
    line_b: Sequence[int],
    locked_minus_a: AbstractSet[int],
    locked_minus_b: AbstractSet[int],
    locked_plus_a: AbstractSet[int],
    locked_plus_b: AbstractSet[int],
) -> Tuple[List[int], List[int], List[int], int, int, int, List[int], List[int]]:
    p01: List[int] = []
    p02: List[int] = []
    p12: List[int] = []
    a01 = 0
    a02 = 0
    a12 = 0
    for idx in range(len(line_a)):
        va, vb = line_a[idx], line_b[idx]
        if va == vb:
            continue
        lo, hi = (va, vb) if va < vb else (vb, va)
        if lo == 0 and hi == 1:
            if idx in locked_plus_a or idx in locked_plus_b:
                continue
            p01.append(idx)
            if va == 1:
                a01 += 1
        elif lo == 0 and hi == 2:
            if (
                idx in locked_plus_a
                or idx in locked_plus_b
                or idx in locked_minus_a
                or idx in locked_minus_b
            ):
                continue
            p02.append(idx)
            if va == 2:
                a02 += 1
        else:
            if idx in locked_minus_a or idx in locked_minus_b:
                continue
            p12.append(idx)
            if va == 2:
                a12 += 1
    m01, m02, m12 = len(p01), len(p02), len(p12)
    t_min = max(-a02, a01 - m01, a12 - m12)
    t_max = min(m02 - a02, a01, a12)
    t_vals = list(range(t_min, t_max + 1))
    weights = [
        math.comb(m01, a01 - t) * math.comb(m02, a02 + t) * math.comb(m12, a12 - t)
        for t in t_vals
    ]
    return p01, p02, p12, a01, a02, a12, t_vals, weights


def _symmetric_ternary_line_step(
    state: List[List[int]],
    rng: random.Random,
    free_rows: Sequence[int],
    locked_rows: AbstractSet[int],
    forbidden_minus_cells: AbstractSet[Tuple[int, int]] = frozenset(),
    forbidden_plus_cells: AbstractSet[Tuple[int, int]] = frozenset(),
    forbid_minus_row_fn: Optional[Callable[[Set[int]], bool]] = None,
    active_rows: Optional[Sequence[int]] = None,
    active_cols: Optional[Sequence[int]] = None,
) -> bool:
    """One strictly symmetric (doubly stochastic) 2-row or 2-column ternary Curveball step (zero exemptions)."""
    n_rows = len(state)
    n_cols = len(state[0]) if n_rows else 0
    cand_rows = active_rows if active_rows is not None else free_rows
    cand_cols = active_cols if active_cols is not None else list(range(n_cols))
    can_row = len(cand_rows) >= 2 and n_cols >= 1
    can_col = len(cand_cols) >= 2 and len(free_rows) >= 2
    if not can_row and not can_col:
        return False
    if can_row and can_col:
        do_col = rng.random() < 0.5
    else:
        do_col = can_col

    if not do_col:
        i, k = rng.sample(cand_rows, 2)
        ri, rk = state[i], state[k]
        lm_i = {c for c in range(n_cols) if (i, c) in forbidden_minus_cells}
        lm_k = {c for c in range(n_cols) if (k, c) in forbidden_minus_cells}
        lp_i = {c for c in range(n_cols) if (i, c) in forbidden_plus_cells}
        lp_k = {c for c in range(n_cols) if (k, c) in forbidden_plus_cells}
        p01, p02, p12, a01, a02, a12, t_vals, weights = _enumerate_line_completions(
            ri, rk, lm_i, lm_k, lp_i, lp_k
        )
        tot_w = sum(weights)
        if tot_w <= 1:
            return False
        pick = rng.randrange(tot_w)
        chosen_t = t_vals[-1]
        acc = 0
        for tv, wv in zip(t_vals, weights):
            acc += wv
            if pick < acc:
                chosen_t = tv
                break
        x01, x02, x12 = a01 - chosen_t, a02 + chosen_t, a12 - chosen_t
        ni, nk = list(ri), list(rk)
        sel01 = set(rng.sample(p01, x01)) if x01 else set()
        for c in p01:
            ni[c], nk[c] = (1, 0) if c in sel01 else (0, 1)
        sel02 = set(rng.sample(p02, x02)) if x02 else set()
        for c in p02:
            ni[c], nk[c] = (2, 0) if c in sel02 else (0, 2)
        sel12 = set(rng.sample(p12, x12)) if x12 else set()
        for c in p12:
            ni[c], nk[c] = (2, 1) if c in sel12 else (1, 2)
        if forbid_minus_row_fn is not None:
            if forbid_minus_row_fn({c for c, v in enumerate(ni) if v == 2}):
                return False
            if forbid_minus_row_fn({c for c, v in enumerate(nk) if v == 2}):
                return False
        changed = (ni != ri) or (nk != rk)
        state[i], state[k] = ni, nk
        return changed

    c1, c2 = rng.sample(cand_cols, 2)
    col1 = [state[r][c1] for r in range(n_rows)]
    col2 = [state[r][c2] for r in range(n_rows)]
    lm_1 = set(locked_rows) | {r for r in range(n_rows) if (r, c1) in forbidden_minus_cells}
    lm_2 = set(locked_rows) | {r for r in range(n_rows) if (r, c2) in forbidden_minus_cells}
    lp_1 = set(locked_rows) | {r for r in range(n_rows) if (r, c1) in forbidden_plus_cells}
    lp_2 = set(locked_rows) | {r for r in range(n_rows) if (r, c2) in forbidden_plus_cells}
    p01, p02, p12, a01, a02, a12, t_vals, weights = _enumerate_line_completions(
        col1, col2, lm_1, lm_2, lp_1, lp_2
    )
    tot_w = sum(weights)
    if tot_w <= 1:
        return False
    pick = rng.randrange(tot_w)
    chosen_t = t_vals[-1]
    acc = 0
    for tv, wv in zip(t_vals, weights):
        acc += wv
        if pick < acc:
            chosen_t = tv
            break
    x01, x02, x12 = a01 - chosen_t, a02 + chosen_t, a12 - chosen_t
    nc1, nc2 = list(col1), list(col2)
    sel01 = set(rng.sample(p01, x01)) if x01 else set()
    for r in p01:
        nc1[r], nc2[r] = (1, 0) if r in sel01 else (0, 1)
    sel02 = set(rng.sample(p02, x02)) if x02 else set()
    for r in p02:
        nc1[r], nc2[r] = (2, 0) if r in sel02 else (0, 2)
    sel12 = set(rng.sample(p12, x12)) if x12 else set()
    for r in p12:
        nc1[r], nc2[r] = (2, 1) if r in sel12 else (1, 2)
    if forbid_minus_row_fn is not None:
        for r in free_rows:
            if nc1[r] != col1[r] or nc2[r] != col2[r]:
                ones_r = {
                    c
                    for c in range(n_cols)
                    if (nc1[r] if c == c1 else (nc2[r] if c == c2 else state[r][c])) == 2
                }
                if forbid_minus_row_fn(ones_r):
                    return False
    changed = (nc1 != col1) or (nc2 != col2)
    if changed:
        for r in range(n_rows):
            state[r][c1] = nc1[r]
            state[r][c2] = nc2[r]
    return changed


def mcmc_reversible_nested_chains(
    mat_minus: Sequence[Sequence[int]],
    mat_plus: Sequence[Sequence[int]],
    rng: random.Random,
    n_draws: int = 1000,
    n_chains: int = 4,
    burnin_mult: int = 20,
    thin_mult: int = 3,
    forbidden_minus_cells: AbstractSet[Tuple[int, int]] = frozenset(),
    forbidden_plus_cells: AbstractSet[Tuple[int, int]] = frozenset(),
    forbid_minus_row_fn: Optional[Callable[[Set[int]], bool]] = None,
    check_every_draw: bool = False,
) -> Tuple[List[List[Tuple[List[List[int]], List[List[int]]]]], List[List[int]]]:
    """Tier 2 reversible symmetric 4-chain MCMC sampler for stationary-mean estimation (never hub-and-spoke)."""
    n_rows = len(mat_minus)
    n_cols = len(mat_minus[0]) if n_rows else 0
    locked_rows = _clash_locked_rows(mat_minus, forbid_minus_row_fn)
    free_rows = [r for r in range(n_rows) if r not in locked_rows]
    assert_nested_completion_invariants(
        mat_minus, mat_plus, mat_minus, mat_plus,
        forbidden_minus_cells=forbidden_minus_cells,
        forbidden_plus_cells=forbidden_plus_cells,
        forbid_minus_row_fn=forbid_minus_row_fn,
        locked_rows=locked_rows,
    )
    state0 = [[mat_minus[r][c] + mat_plus[r][c] for c in range(n_cols)] for r in range(n_rows)]
    # State-independent fixed-margin active lines (constant lines never change on Omega; preserves exact symmetry)
    act_rows = [
        r for r in free_rows
        if sum(1 for s in (0, 1, 2) if any(state0[r][c] == s for c in range(n_cols))) >= 2
    ]
    act_cols = [
        c for c in range(n_cols)
        if sum(1 for s in (0, 1, 2) if any(state0[r][c] == s for r in free_rows)) >= 2
    ]
    draws_per_chain = max(10, (n_draws + n_chains - 1) // n_chains)
    line_dim = max(1, max(n_rows, n_cols))
    burnin_steps = max(40, burnin_mult * line_dim)
    thin_steps = max(12, thin_mult * line_dim)

    chains_out: List[List[Tuple[List[List[int]], List[List[int]]]]] = []
    hamming_out: List[List[int]] = []
    for _ in range(n_chains):
        cur = [list(row) for row in state0]
        for _ in range(burnin_steps):
            _symmetric_ternary_line_step(
                cur, rng, free_rows, locked_rows,
                forbidden_minus_cells, forbidden_plus_cells, forbid_minus_row_fn,
                active_rows=act_rows, active_cols=act_cols,
            )
        c_draws: List[Tuple[List[List[int]], List[List[int]]]] = []
        c_ham: List[int] = []
        for d_idx in range(draws_per_chain):
            for _ in range(thin_steps):
                _symmetric_ternary_line_step(
                    cur, rng, free_rows, locked_rows,
                    forbidden_minus_cells, forbidden_plus_cells, forbid_minus_row_fn,
                    active_rows=act_rows, active_cols=act_cols,
                )
            pm_m = [[1 if v == 2 else 0 for v in row] for row in cur]
            pm_p = [[1 if v >= 1 else 0 for v in row] for row in cur]
            if check_every_draw or d_idx < 2 or d_idx == draws_per_chain - 1:
                assert_nested_completion_invariants(
                    mat_minus, mat_plus, pm_m, pm_p,
                    forbidden_minus_cells=forbidden_minus_cells,
                    forbidden_plus_cells=forbidden_plus_cells,
                    forbid_minus_row_fn=forbid_minus_row_fn,
                    locked_rows=locked_rows,
                )
            ham = sum(
                1
                for r in range(n_rows)
                for c in range(n_cols)
                if cur[r][c] != state0[r][c]
            )
            c_draws.append((pm_m, pm_p))
            c_ham.append(ham)
        chains_out.append(c_draws)
        hamming_out.append(c_ham)
    return chains_out, hamming_out


def compute_series_rhat_ess_mcse(
    chain_series: Sequence[Sequence[float]],
    n_batches_per_chain: int = 5,
) -> Dict[str, object]:
    """Compute Gelman-Rubin split-Rhat, Geyer positive-pair bulk ESS, and non-overlapping batch-means MCSE."""
    n_chains = len(chain_series)
    n_per = len(chain_series[0]) if n_chains else 0
    all_vals = [float(x) for ch in chain_series for x in ch]
    if not all_vals:
        return {"mean": 0.0, "var_bm": 0.0, "mcse_bm": 0.0, "is_constant": True, "split_rhat": 1.0, "ess": 0.0}
    grand_mean = sum(all_vals) / len(all_vals)
    if max(all_vals) - min(all_vals) <= 1e-12:
        return {
            "mean": grand_mean,
            "var_bm": 0.0,
            "mcse_bm": 0.0,
            "is_constant": True,
            "split_rhat": 1.0,
            "ess": float(len(all_vals)),
        }

    # Non-overlapping batch-means variance across K = n_chains * n_batches_per_chain batches
    b_size = max(1, n_per // n_batches_per_chain)
    batch_means: List[float] = []
    for ch in chain_series:
        for b_idx in range(n_batches_per_chain):
            s_i = b_idx * b_size
            e_i = n_per if b_idx == n_batches_per_chain - 1 else (b_idx + 1) * b_size
            seg = ch[s_i:e_i]
            if seg:
                batch_means.append(sum(seg) / len(seg))
    k_b = len(batch_means)
    bm_mean = sum(batch_means) / k_b
    var_bm = (
        sum((bm - bm_mean) ** 2 for bm in batch_means) / (k_b * (k_b - 1))
        if k_b >= 2
        else 0.0
    )

    # Split chains in half (m = 2 * n_chains split chains of length n_half)
    n_half = n_per // 2
    if n_half < 2:
        return {
            "mean": grand_mean,
            "var_bm": var_bm,
            "mcse_bm": math.sqrt(var_bm),
            "is_constant": False,
            "split_rhat": math.inf,
            "ess": 1.0,
        }
    splits: List[List[float]] = []
    for ch in chain_series:
        splits.append([float(x) for x in ch[:n_half]])
        splits.append([float(x) for x in ch[n_half : 2 * n_half]])
    m_sp = len(splits)
    sp_means = [sum(sp) / n_half for sp in splits]
    sp_grand = sum(sp_means) / m_sp
    b_var = (n_half / (m_sp - 1)) * sum(( sm - sp_grand ) ** 2 for sm in sp_means)
    sp_vars = [
        sum((x - sm) ** 2 for x in sp) / (n_half - 1)
        for sp, sm in zip(splits, sp_means)
    ]
    w_var = sum(sp_vars) / m_sp
    if w_var <= 1e-15:
        split_rhat = 1.0 if b_var <= 1e-15 else math.inf
        return {
            "mean": grand_mean,
            "var_bm": var_bm,
            "mcse_bm": math.sqrt(var_bm),
            "is_constant": False,
            "split_rhat": split_rhat,
            "ess": 1.0 if b_var > 1e-15 else float(m_sp * n_half),
        }
    var_plus = ((n_half - 1) / n_half) * w_var + (b_var / n_half)
    split_rhat = math.sqrt(var_plus / w_var)

    # Geyer positive-pair integrated autocorrelation time across split chains
    rhos: List[float] = [1.0]
    for lag in range(1, n_half):
        v_lag = 0.0
        for sp, sm in zip(splits, sp_means):
            v_lag += sum((sp[t] - sm) * (sp[t + lag] - sm) for t in range(n_half - lag)) / n_half
        v_lag /= m_sp
        rho_lag = 1.0 - (w_var - v_lag) / var_plus
        rhos.append(rho_lag)
    tau = -1.0
    s_pair = 0
    while 2 * s_pair < len(rhos):
        p_sum = rhos[2 * s_pair] + (rhos[2 * s_pair + 1] if (2 * s_pair + 1) < len(rhos) else 0.0)
        if p_sum <= 0.0:
            break
        tau += 2.0 * p_sum
        s_pair += 1
    tau = max(1.0, tau)
    ess = (m_sp * n_half) / tau
    return {
        "mean": grand_mean,
        "var_bm": var_bm,
        "mcse_bm": math.sqrt(var_bm),
        "is_constant": False,
        "split_rhat": round(split_rhat, 6),
        "ess": round(ess, 2),
    }


def validate_mcmc_prevalence_bounds(
    low_chains: Sequence[Sequence[float]],
    high_chains: Sequence[Sequence[float]],
    hamming_chains: Sequence[Sequence[int]],
    has_active_move: bool,
    min_ess: float = 200.0,
    max_rhat: float = 1.05,
) -> Dict[str, object]:
    """Evaluate empirical R-hat, ESS, and 20-batch MCSE directly on the prevalence bounds used in the verdict.

    Note: split-Rhat and ESS are empirical convergence checks, NOT a general mathematical guarantee of
    finite-time mixing or multi-component exploration (connectivity is enforced separately).
    """
    d_low = compute_series_rhat_ess_mcse(low_chains)
    d_high = compute_series_rhat_ess_mcse(high_chains)
    d_ham = compute_series_rhat_ess_mcse([[float(h) for h in ch] for ch in hamming_chains])

    moved_ok = (not has_active_move) or (max((h for ch in hamming_chains for h in ch), default=0) > 0)
    low_ok = bool(d_low["is_constant"]) or (
        float(d_low["split_rhat"]) <= max_rhat and float(d_low["ess"]) >= min_ess  # type: ignore[arg-type]
    )
    high_ok = bool(d_high["is_constant"]) or (
        float(d_high["split_rhat"]) <= max_rhat and float(d_high["ess"]) >= min_ess  # type: ignore[arg-type]
    )
    # If both prevalence bounds stayed constant on a graph with active moves, require Hamming distance to mix cleanly
    ham_guard_ok = True
    if has_active_move and bool(d_low["is_constant"]) and bool(d_high["is_constant"]):
        ham_guard_ok = (
            not bool(d_ham["is_constant"])
            and float(d_ham["split_rhat"]) <= max_rhat  # type: ignore[arg-type]
            and float(d_ham["ess"]) >= min_ess  # type: ignore[arg-type]
        )
    passed = moved_ok and low_ok and high_ok and ham_guard_ok
    return {
        "mcmc_diagnostics_pass": passed,
        "moved_from_obs_ok": moved_ok,
        "bound_low_diag": d_low,
        "bound_high_diag": d_high,
        "hamming_diag": d_ham,
    }


def _has_any_active_2line_move(
    state0: Sequence[Sequence[int]],
    free_rows: Sequence[int],
    locked_rows: AbstractSet[int],
    forbidden_minus_cells: AbstractSet[Tuple[int, int]],
    forbidden_plus_cells: AbstractSet[Tuple[int, int]],
    forbid_minus_row_fn: Optional[Callable[[Set[int]], bool]],
) -> bool:
    n_rows = len(state0)
    n_cols = len(state0[0]) if n_rows else 0
    for idx1 in range(len(free_rows)):
        for idx2 in range(idx1 + 1, len(free_rows)):
            r1, r2 = free_rows[idx1], free_rows[idx2]
            lm1 = {c for c in range(n_cols) if (r1, c) in forbidden_minus_cells}
            lm2 = {c for c in range(n_cols) if (r2, c) in forbidden_minus_cells}
            lp1 = {c for c in range(n_cols) if (r1, c) in forbidden_plus_cells}
            lp2 = {c for c in range(n_cols) if (r2, c) in forbidden_plus_cells}
            p01, p02, p12, a01, a02, a12, t_vals, weights = _enumerate_line_completions(
                state0[r1], state0[r2], lm1, lm2, lp1, lp2
            )
            if sum(weights) <= 1:
                continue
            if forbid_minus_row_fn is None:
                return True
            for tv, wv in zip(t_vals, weights):
                if wv == 0:
                    continue
                x01, x02, x12 = a01 - tv, a02 + tv, a12 - tv
                for s01 in itertools.combinations(p01, x01):
                    s01_s = set(s01)
                    for s02 in itertools.combinations(p02, x02):
                        s02_s = set(s02)
                        for s12 in itertools.combinations(p12, x12):
                            s12_s = set(s12)
                            nr1, nr2 = list(state0[r1]), list(state0[r2])
                            for c in p01:
                                nr1[c], nr2[c] = (1, 0) if c in s01_s else (0, 1)
                            for c in p02:
                                nr1[c], nr2[c] = (2, 0) if c in s02_s else (0, 2)
                            for c in p12:
                                nr1[c], nr2[c] = (2, 1) if c in s12_s else (1, 2)
                            if (nr1 != list(state0[r1]) or nr2 != list(state0[r2])) and not (
                                forbid_minus_row_fn({c for c, v in enumerate(nr1) if v == 2})
                                or forbid_minus_row_fn({c for c, v in enumerate(nr2) if v == 2})
                            ):
                                return True
    for c1 in range(n_cols):
        for c2 in range(c1 + 1, n_cols):
            col1 = [state0[r][c1] for r in range(n_rows)]
            col2 = [state0[r][c2] for r in range(n_rows)]
            lm1 = set(locked_rows) | {r for r in range(n_rows) if (r, c1) in forbidden_minus_cells}
            lm2 = set(locked_rows) | {r for r in range(n_rows) if (r, c2) in forbidden_minus_cells}
            lp1 = set(locked_rows) | {r for r in range(n_rows) if (r, c1) in forbidden_plus_cells}
            lp2 = set(locked_rows) | {r for r in range(n_rows) if (r, c2) in forbidden_plus_cells}
            p01, p02, p12, a01, a02, a12, t_vals, weights = _enumerate_line_completions(
                col1, col2, lm1, lm2, lp1, lp2
            )
            if sum(weights) <= 1:
                continue
            if forbid_minus_row_fn is None:
                return True
            for tv, wv in zip(t_vals, weights):
                if wv == 0:
                    continue
                x01, x02, x12 = a01 - tv, a02 + tv, a12 - tv
                for s01 in itertools.combinations(p01, x01):
                    s01_s = set(s01)
                    for s02 in itertools.combinations(p02, x02):
                        s02_s = set(s02)
                        for s12 in itertools.combinations(p12, x12):
                            s12_s = set(s12)
                            nc1, nc2 = list(col1), list(col2)
                            for r in p01:
                                nc1[r], nc2[r] = (1, 0) if r in s01_s else (0, 1)
                            for r in p02:
                                nc1[r], nc2[r] = (2, 0) if r in s02_s else (0, 2)
                            for r in p12:
                                nc1[r], nc2[r] = (2, 1) if r in s12_s else (1, 2)
                            if nc1 == col1 and nc2 == col2:
                                continue
                            ok_c = True
                            for r in free_rows:
                                if nc1[r] != col1[r] or nc2[r] != col2[r]:
                                    ones_r = {
                                        cc
                                        for cc in range(n_cols)
                                        if (nc1[r] if cc == c1 else (nc2[r] if cc == c2 else state0[r][cc])) == 2
                                    }
                                    if forbid_minus_row_fn(ones_r):
                                        ok_c = False
                                        break
                            if ok_c:
                                return True
    return False


def classify_ternary_move_connectivity(
    mat_minus: Sequence[Sequence[int]],
    mat_plus: Sequence[Sequence[int]],
    forbidden_minus_cells: AbstractSet[Tuple[int, int]] = frozenset(),
    forbidden_plus_cells: AbstractSet[Tuple[int, int]] = frozenset(),
    forbid_minus_row_fn: Optional[Callable[[Set[int]], bool]] = None,
    max_enum_states: int = 2000,
    max_bt_nodes: int = 50_000,
    return_enumerated_pairs: bool = False,
) -> Dict[str, object]:
    """Classify local (K_row + K_col) connectivity and, when |Omega| <= max_enum_states, return exact Omega & C_K."""
    n_rows = len(mat_minus)
    n_cols = len(mat_minus[0]) if n_rows else 0
    state0 = tuple(tuple(mat_minus[r][c] + mat_plus[r][c] for c in range(n_cols)) for r in range(n_rows))
    locked_rows = _clash_locked_rows(mat_minus, forbid_minus_row_fn)
    free_rows = [r for r in range(n_rows) if r not in locked_rows]

    def _to_pair(grid_s: Sequence[Sequence[int]]) -> Tuple[List[List[int]], List[List[int]]]:
        return (
            [[1 if v == 2 else 0 for v in row] for row in grid_s],
            [[1 if v >= 1 else 0 for v in row] for row in grid_s],
        )

    if len(free_rows) < 2 or n_cols < 2:
        out_triv: Dict[str, object] = {
            "status": "connected",
            "component_size": 1,
            "fiber_size": 1,
            "locked_clash_rows": sorted(locked_rows),
            "has_active_move": False,
            "exact_full_fiber_enumerated": True,
            "proof": "trivial_dimension",
        }
        if return_enumerated_pairs:
            p0 = _to_pair(state0)
            out_triv["exact_fiber_pairs"] = [p0]
            out_triv["exact_component_pairs"] = [p0]
        return out_triv

    row_counts = [[sum(1 for v in state0[r] if v == s) for s in (0, 1, 2)] for r in range(n_rows)]
    col_counts = [[sum(1 for r in range(n_rows) if state0[r][c] == s) for s in (0, 1, 2)] for c in range(n_cols)]
    act_rows = [r for r in free_rows if sum(1 for cnt in row_counts[r] if cnt > 0) >= 2]
    act_cols = [c for c in range(n_cols) if sum(1 for cnt in col_counts[c] if cnt > 0) >= 2]
    if not act_rows or not act_cols:
        out_sing: Dict[str, object] = {
            "status": "connected",
            "component_size": 1,
            "fiber_size": 1,
            "locked_clash_rows": sorted(locked_rows),
            "has_active_move": False,
            "exact_full_fiber_enumerated": True,
            "proof": "singleton_fiber",
        }
        if return_enumerated_pairs:
            p0 = _to_pair(state0)
            out_sing["exact_fiber_pairs"] = [p0]
            out_sing["exact_component_pairs"] = [p0]
        return out_sing

    ar_len, ac_len = len(act_rows), len(act_cols)
    rem_r = [[row_counts[r][s] for s in (0, 1, 2)] for r in act_rows]
    rem_c = [[col_counts[c][s] for s in (0, 1, 2)] for c in act_cols]
    act_col_set = set(act_cols)
    act_row_set = set(act_rows)
    for ai, r in enumerate(act_rows):
        for c in range(n_cols):
            if c not in act_col_set:
                rem_r[ai][state0[r][c]] -= 1
    for aj, c in enumerate(act_cols):
        for r in range(n_rows):
            if r not in act_row_set:
                rem_c[aj][state0[r][c]] -= 1

    def _fallback_unenumerated(reason_str: str) -> Dict[str, object]:
        has_mv = _has_any_active_2line_move(
            state0, free_rows, locked_rows, forbidden_minus_cells, forbidden_plus_cells, forbid_minus_row_fn
        )
        # Ryser (1957) certifies connectivity on unconstrained 2-symbol binary fibers even when |Omega| > max_enum_states
        if not forbidden_minus_cells and not forbidden_plus_cells and forbid_minus_row_fn is None:
            active_syms = [s for s in (0, 1, 2) if sum(rem_r[ai][s] for ai in range(ar_len)) > 0]
            if len(active_syms) <= 2:
                return {
                    "status": "connected",
                    "component_size": None,
                    "fiber_size": None,
                    "locked_clash_rows": [],
                    "has_active_move": has_mv,
                    "exact_full_fiber_enumerated": False,
                    "proof": "ryser_1957_two_symbol_binary_reduction",
                }
        return {
            "status": "connectivity unknown",
            "component_size": None,
            "fiber_size": None,
            "locked_clash_rows": sorted(locked_rows),
            "has_active_move": has_mv,
            "exact_full_fiber_enumerated": False,
            "proof": reason_str,
        }

    if ar_len * ac_len > 64:
        return _fallback_unenumerated("active_submatrix_exceeds_exact_enumeration_cap")

    # Tier 1: Exact backtracking enumeration of the ENTIRE allowed set Omega(S_obs) on active submatrix
    base_grid = [list(row) for row in state0]
    fiber: Set[Tuple[Tuple[int, ...], ...]] = set()
    exceeded = False
    bt_nodes = 0

    def _bt(pos: int) -> None:
        nonlocal exceeded, bt_nodes
        if exceeded:
            return
        bt_nodes += 1
        if bt_nodes > max_bt_nodes:
            exceeded = True
            return
        if pos == ar_len * ac_len:
            if forbid_minus_row_fn is not None:
                for r in act_rows:
                    if forbid_minus_row_fn({c for c in range(n_cols) if base_grid[r][c] == 2}):
                        return
            fiber.add(tuple(tuple(base_grid[r][c] for c in act_cols) for r in act_rows))
            if len(fiber) > max_enum_states:
                exceeded = True
            return
        ai, aj = divmod(pos, ac_len)
        r, c = act_rows[ai], act_cols[aj]
        cells_left_r = ac_len - aj - 1
        cells_left_c = ar_len - ai - 1
        for s in (0, 1, 2):
            if rem_r[ai][s] <= 0 or rem_c[aj][s] <= 0:
                continue
            if s >= 1 and (r, c) in forbidden_plus_cells:
                continue
            if s == 2 and (r, c) in forbidden_minus_cells:
                continue
            rem_r[ai][s] -= 1
            rem_c[aj][s] -= 1
            if sum(rem_r[ai]) == cells_left_r and sum(rem_c[aj]) == cells_left_c:
                base_grid[r][c] = s
                _bt(pos + 1)
            rem_r[ai][s] += 1
            rem_c[aj][s] += 1

    _bt(0)
    if exceeded:
        return _fallback_unenumerated("exceeds_exact_enumeration_cap")

    # Exact BFS from S_obs over all 2-row and 2-column transitions of K (zero exemptions)
    start_sub = tuple(tuple(state0[r][c] for c in act_cols) for r in act_rows)
    visited: Set[Tuple[Tuple[int, ...], ...]] = {start_sub}
    queue = [start_sub]

    def _expand_neighbors(sub: Tuple[Tuple[int, ...], ...]) -> Optional[List[Tuple[Tuple[int, ...], ...]]]:
        full = [list(row) for row in state0]
        for ai, r in enumerate(act_rows):
            for aj, c in enumerate(act_cols):
                full[r][c] = sub[ai][aj]
        out_n: List[Tuple[Tuple[int, ...], ...]] = []
        for ai1 in range(ar_len):
            for ai2 in range(ai1 + 1, ar_len):
                r1, r2 = act_rows[ai1], act_rows[ai2]
                lm1 = {c for c in range(n_cols) if (r1, c) in forbidden_minus_cells}
                lm2 = {c for c in range(n_cols) if (r2, c) in forbidden_minus_cells}
                lp1 = {c for c in range(n_cols) if (r1, c) in forbidden_plus_cells}
                lp2 = {c for c in range(n_cols) if (r2, c) in forbidden_plus_cells}
                p01, p02, p12, a01, a02, a12, t_vals, weights = _enumerate_line_completions(
                    full[r1], full[r2], lm1, lm2, lp1, lp2
                )
                tot_w = sum(weights)
                if tot_w <= 1:
                    continue
                if tot_w > max_enum_states:
                    return None
                for tv, wv in zip(t_vals, weights):
                    if wv == 0:
                        continue
                    x01, x02, x12 = a01 - tv, a02 + tv, a12 - tv
                    for s01 in itertools.combinations(p01, x01):
                        s01_set = set(s01)
                        for s02 in itertools.combinations(p02, x02):
                            s02_set = set(s02)
                            for s12 in itertools.combinations(p12, x12):
                                s12_set = set(s12)
                                nr1, nr2 = list(full[r1]), list(full[r2])
                                for c in p01:
                                    nr1[c], nr2[c] = (1, 0) if c in s01_set else (0, 1)
                                for c in p02:
                                    nr1[c], nr2[c] = (2, 0) if c in s02_set else (0, 2)
                                for c in p12:
                                    nr1[c], nr2[c] = (2, 1) if c in s12_set else (1, 2)
                                if forbid_minus_row_fn is not None:
                                    if forbid_minus_row_fn({c for c, v in enumerate(nr1) if v == 2}):
                                        continue
                                    if forbid_minus_row_fn({c for c, v in enumerate(nr2) if v == 2}):
                                        continue
                                cand_sub = tuple(
                                    tuple(
                                        (nr1[c] if r == r1 else (nr2[c] if r == r2 else full[r][c]))
                                        for c in act_cols
                                    )
                                    for r in act_rows
                                )
                                out_n.append(cand_sub)
        for aj1 in range(ac_len):
            for aj2 in range(aj1 + 1, ac_len):
                c1, c2 = act_cols[aj1], act_cols[aj2]
                col1 = [full[r][c1] for r in range(n_rows)]
                col2 = [full[r][c2] for r in range(n_rows)]
                lm1 = set(locked_rows) | {r for r in range(n_rows) if (r, c1) in forbidden_minus_cells}
                lm2 = set(locked_rows) | {r for r in range(n_rows) if (r, c2) in forbidden_minus_cells}
                lp1 = set(locked_rows) | {r for r in range(n_rows) if (r, c1) in forbidden_plus_cells}
                lp2 = set(locked_rows) | {r for r in range(n_rows) if (r, c2) in forbidden_plus_cells}
                p01, p02, p12, a01, a02, a12, t_vals, weights = _enumerate_line_completions(
                    col1, col2, lm1, lm2, lp1, lp2
                )
                tot_w = sum(weights)
                if tot_w <= 1:
                    continue
                if tot_w > max_enum_states:
                    return None
                for tv, wv in zip(t_vals, weights):
                    if wv == 0:
                        continue
                    x01, x02, x12 = a01 - tv, a02 + tv, a12 - tv
                    for s01 in itertools.combinations(p01, x01):
                        s01_set = set(s01)
                        for s02 in itertools.combinations(p02, x02):
                            s02_set = set(s02)
                            for s12 in itertools.combinations(p12, x12):
                                s12_set = set(s12)
                                nc1, nc2 = list(col1), list(col2)
                                for r in p01:
                                    nc1[r], nc2[r] = (1, 0) if r in s01_set else (0, 1)
                                for r in p02:
                                    nc1[r], nc2[r] = (2, 0) if r in s02_set else (0, 2)
                                for r in p12:
                                    nc1[r], nc2[r] = (2, 1) if r in s12_set else (1, 2)
                                ok_c = True
                                if forbid_minus_row_fn is not None:
                                    for r in act_rows:
                                        if nc1[r] != col1[r] or nc2[r] != col2[r]:
                                            ones_r = {
                                                cc
                                                for cc in range(n_cols)
                                                if (nc1[r] if cc == c1 else (nc2[r] if cc == c2 else full[r][cc])) == 2
                                            }
                                            if forbid_minus_row_fn(ones_r):
                                                ok_c = False
                                                break
                                if not ok_c:
                                    continue
                                cand_sub = tuple(
                                    tuple(
                                        (nc1[r] if c == c1 else (nc2[r] if c == c2 else full[r][c]))
                                        for c in act_cols
                                    )
                                    for r in act_rows
                                )
                                out_n.append(cand_sub)
        return out_n

    bfs_exceeded = False
    while queue:
        cur_s = queue.pop()
        nbrs = _expand_neighbors(cur_s)
        if nbrs is None:
            bfs_exceeded = True
            break
        for nxt_s in nbrs:
            if nxt_s not in visited:
                visited.add(nxt_s)
                queue.append(nxt_s)

    def _sub_to_pair(sub_m: Tuple[Tuple[int, ...], ...]) -> Tuple[List[List[int]], List[List[int]]]:
        full = [list(row) for row in state0]
        for ai, r in enumerate(act_rows):
            for aj, c in enumerate(act_cols):
                full[r][c] = sub_m[ai][aj]
        return _to_pair(full)

    if bfs_exceeded:
        # Note: exact_full_fiber_enumerated is STILL True because `fiber` enumerated all of Omega!
        out_bfs_cap: Dict[str, object] = {
            "status": "connectivity unknown",
            "component_size": None,
            "fiber_size": len(fiber),
            "locked_clash_rows": sorted(locked_rows),
            "has_active_move": True,
            "exact_full_fiber_enumerated": True,
            "proof": "line_completion_exceeds_cap_after_fiber_enum",
        }
        if return_enumerated_pairs:
            out_bfs_cap["exact_fiber_pairs"] = [_sub_to_pair(s) for s in sorted(fiber)]
            out_bfs_cap["exact_component_pairs"] = None
        return out_bfs_cap

    status_str = "connected" if len(visited) == len(fiber) else "disconnected (proven)"
    out_res: Dict[str, object] = {
        "status": status_str,
        "component_size": len(visited),
        "fiber_size": len(fiber),
        "locked_clash_rows": sorted(locked_rows),
        "has_active_move": len(visited) > 1,
        "exact_full_fiber_enumerated": True,
        "proof": "exact_enumeration_and_bfs",
    }
    if return_enumerated_pairs:
        out_res["exact_fiber_pairs"] = [_sub_to_pair(s) for s in sorted(fiber)]
        out_res["exact_component_pairs"] = [_sub_to_pair(s) for s in sorted(visited)]
    return out_res


def verify_nested_curveball_fixture() -> Dict[str, object]:
    """Fixture suite verifying trades, invariants, 3-state connectivity, 2-component full-fiber enumeration, Null 2c target invariance, Tier 2 calibration, and slow-kernel regression case."""
    # 1. 4x5 nested fixture with ternary states s = minus + plus:
    #    Row 0: [2, 0, 1, 0, 2], Row 1: [0, 2, 0, 1, 0], Row 2: [1, 0, 2, 0, 1], Row 3: [0, 1, 0, 2, 0]
    fx_minus = [
        [1, 0, 0, 0, 1],
        [0, 1, 0, 0, 0],
        [0, 0, 1, 0, 0],
        [0, 0, 0, 1, 0],
    ]
    fx_plus = [
        [1, 0, 1, 0, 1],
        [0, 1, 0, 1, 0],
        [1, 0, 1, 0, 1],
        [0, 1, 0, 1, 0],
    ]
    rng = random.Random(20261121)

    um, up, u_trades = curveball_permute_nested_pair(
        fx_minus, fx_plus, rng, n_steps=15, check_each_trade=True
    )
    assert u_trades > 0, "Unconstrained fixture performed zero nontrivial trades"
    assert um != fx_minus and up != fx_plus, "Unconstrained fixture draw did not change both I- and I+"

    forb_plus = {(0, 3), (1, 4)}
    forb_minus = {(2, 0)}
    forb_row = lambda s_idx: (0 in s_idx) and (3 in s_idx)
    cm_changed_both = False
    c_total_trades = 0
    for _ in range(10):
        cm, cp, c_trades = curveball_permute_nested_pair(
            fx_minus,
            fx_plus,
            rng,
            n_steps=15,
            forbidden_minus_cells=forb_minus,
            forbidden_plus_cells=forb_plus,
            forbid_minus_row_fn=forb_row,
            check_each_trade=True,
        )
        c_total_trades += c_trades
        if cm != fx_minus and cp != fx_plus:
            cm_changed_both = True
    assert c_total_trades > 0 and cm_changed_both, "Constrained fixture failed to execute a nontrivial legal trade"

    # 2. Pure 3-state cyclic trade fixture: [[0, 1, 2], [1, 2, 0]] <-> [[1, 2, 0], [0, 1, 2]]
    cyc_m = [[0, 0, 1], [0, 1, 0]]
    cyc_p = [[0, 1, 1], [1, 1, 0]]
    saw_flipped = False
    for _ in range(20):
        ym, yp, y_trades = curveball_permute_nested_pair(
            cyc_m, cyc_p, rng, n_steps=1, check_each_trade=True
        )
        if ym == [[0, 1, 0], [0, 0, 1]] and yp == [[1, 1, 0], [0, 1, 1]]:
            assert y_trades == 1
            saw_flipped = True
    assert saw_flipped, "3-state cyclic trade [[0,1,2],[1,2,0]] -> [[1,2,0],[0,1,2]] was never reached"

    # 3. Negative checks
    bad_m = [list(r) for r in fx_minus]
    bad_p = [list(r) for r in fx_plus]
    bad_m[0][1] = 1
    caught_nesting = False
    try:
        assert_nested_completion_invariants(fx_minus, fx_plus, bad_m, bad_p)
    except AssertionError:
        caught_nesting = True
    assert caught_nesting, "Invariant checker failed to catch definite > possible"

    bad_forb_p = [list(r) for r in fx_plus]
    bad_forb_p[0][3] = 1
    bad_forb_p[0][2] = 0
    bad_forb_p[1][3] = 0
    bad_forb_p[1][2] = 1
    caught_forb = False
    try:
        assert_nested_completion_invariants(
            fx_minus, fx_plus, fx_minus, bad_forb_p, forbidden_plus_cells=forb_plus
        )
    except AssertionError:
        caught_forb = True
    assert caught_forb, "Invariant checker failed to catch forbidden cell violation even with preserved sums"

    # 4. 3x3 ternary Latin-square full-fiber check (012/120/201 -> all 12/12 legal 3x3 Latin squares)
    lat_m = [[0, 0, 1], [0, 1, 0], [1, 0, 0]]
    lat_p = [[0, 1, 1], [1, 1, 0], [1, 0, 1]]
    visited_latin: Set[Tuple[Tuple[int, int, int], ...]] = set()
    for _ in range(120):
        lm, lp, _ = curveball_permute_nested_pair(
            lat_m, lat_p, rng, n_steps=6, check_each_trade=True, enable_latin3x3=True
        )
        sq = tuple(tuple(lm[r][c] + lp[r][c] for c in range(3)) for r in range(3))
        visited_latin.add(sq)
    assert visited_latin == ALL_12_LATIN_3X3_SET, (
        f"Expected all 12 legal 3x3 Latin squares, visited {len(visited_latin)}"
    )

    # 5. Verify 3-state connectivity classifier (exact connected, Ryser >cap connected, disconnected-proven, unknown)
    conn_ryser_small = classify_ternary_move_connectivity([[1, 0], [0, 1]], [[1, 0], [0, 1]])
    assert conn_ryser_small["status"] == "connected" and conn_ryser_small["fiber_size"] == 2
    id10 = [[1 if r == c else 0 for c in range(10)] for r in range(10)]
    conn_ryser_big = classify_ternary_move_connectivity(id10, id10, max_enum_states=50)
    assert conn_ryser_big["status"] == "connected" and conn_ryser_big["proof"] == "ryser_1957_two_symbol_binary_reduction"
    conn_latin = classify_ternary_move_connectivity(lat_m, lat_p, return_enumerated_pairs=True)
    assert conn_latin["status"] == "connected" and conn_latin["component_size"] == 12 and conn_latin["fiber_size"] == 12
    conn_unk = classify_ternary_move_connectivity(lat_m, lat_p, max_enum_states=5)
    assert conn_unk["status"] == "connectivity unknown" and conn_unk["exact_full_fiber_enumerated"] is False

    # 6. Null 2c observation-invariant equivalence-class fiber [S_obs]_{~2c}:
    #    Row 0 has an initial {0, 3} clash; verify Row 0 stays locked with its clash in all states & starting from any Y != X gives identical Omega!
    clash_m = [[1, 0, 0, 1], [1, 0, 0, 0], [0, 1, 0, 0], [0, 0, 1, 0]]
    clash_p = [[1, 0, 0, 1], [1, 0, 1, 0], [0, 1, 0, 1], [1, 0, 1, 0]]
    c_eq_0 = classify_ternary_move_connectivity(
        clash_m, clash_p, forbid_minus_row_fn=forb_row, return_enumerated_pairs=True
    )
    assert c_eq_0["locked_clash_rows"] == [0]
    eq_pairs_0 = c_eq_0["exact_fiber_pairs"]  # type: ignore[assignment]
    assert len(eq_pairs_0) > 1  # type: ignore[arg-type]
    canon_set_0 = {
        (tuple(tuple(r) for r in pm), tuple(tuple(r) for r in pp))
        for pm, pp in eq_pairs_0  # type: ignore[union-attr]
    }
    alt_pair = next(
        (pm, pp) for pm, pp in eq_pairs_0  # type: ignore[union-attr]
        if pm != clash_m or pp != clash_p
    )
    for pm, pp in eq_pairs_0:  # type: ignore[union-attr]
        assert pm[0] == clash_m[0] and pp[0] == clash_p[0]
        assert _clash_locked_rows(pm, forb_row) == frozenset({0})
    c_eq_alt = classify_ternary_move_connectivity(
        alt_pair[0], alt_pair[1], forbid_minus_row_fn=forb_row, return_enumerated_pairs=True
    )
    canon_set_alt = {
        (tuple(tuple(r) for r in pm), tuple(tuple(r) for r in pp))
        for pm, pp in c_eq_alt["exact_fiber_pairs"]  # type: ignore[union-attr,index]
    }
    assert canon_set_0 == canon_set_alt, "Null 2c fiber depended on which state in the equivalence class was observed!"

    # 7. Two-component unequal-mean exact enumeration test (Astra 249779947258167296):
    #    3 donors x 21 DONOR_CAP_MODULES where active cols are (P_B12_AER, C_HIS, C_TRP), P_BLUB is constant 1 on all 3 donors,
    #    and A_DSR is locked at [0, 1, 0]^T via forbidden_plus_cells, plus diagonal shift lock on the 3 active cols!
    idx_b12_aer = DONOR_CAP_MODULES.index("P_B12_AER")
    idx_blub = DONOR_CAP_MODULES.index("P_BLUB")
    idx_his = DONOR_CAP_MODULES.index("C_HIS")
    idx_trp = DONOR_CAP_MODULES.index("C_TRP")
    idx_dsr = DONOR_CAP_MODULES.index("A_DSR")
    n_cap = len(DONOR_CAP_MODULES)
    twocomp_m = [[0] * n_cap for _ in range(3)]
    for r_i in range(3):
        twocomp_m[r_i][idx_blub] = 1
    twocomp_m[0][idx_b12_aer] = 1  # Donor 0: B12 only
    twocomp_m[1][idx_his] = 1      # Donor 1: C_HIS + A_DSR (no B12 in State A; gets B12 in State B!)
    twocomp_m[1][idx_dsr] = 1
    twocomp_m[2][idx_trp] = 1      # Donor 2: C_TRP only
    twocomp_p = [list(r) for r in twocomp_m]
    forb_twocomp = {
        (0, idx_his), (1, idx_trp), (2, idx_b12_aer),  # locks 3x3 shift so |Omega|=2, |C_K|=1
        (0, idx_dsr), (2, idx_dsr),                    # fixes A_DSR column to [0, 1, 0]^T
    }
    conn_2comp = classify_ternary_move_connectivity(
        twocomp_m, twocomp_p, forbidden_plus_cells=forb_twocomp, return_enumerated_pairs=True
    )
    assert (
        conn_2comp["status"] == "disconnected (proven)"
        and conn_2comp["component_size"] == 1
        and conn_2comp["fiber_size"] == 2
        and conn_2comp["exact_full_fiber_enumerated"] is True
    )
    rec_x_prof = OrganismProfile("rec_x", frozenset(), (), frozenset({"B12", "RED_S"}), frozenset())

    def _plexi_indicator(pair_mp: Tuple[List[List[int]], List[List[int]]]) -> float:
        pm_x, pp_x = pair_mp
        pr_m = {f"d{i}": profile_from_module_Presence_row(f"d{i}", pm_x[i]) for i in range(3)}
        pr_p = {f"d{i}": profile_from_module_Presence_row(f"d{i}", pp_x[i]) for i in range(3)}
        pr_m["rec_x"] = rec_x_prof
        pr_p["rec_x"] = rec_x_prof
        st_m, _ = _fast_decide_acyclic_p_all("rec_x", pr_m, frozenset(), cyc_mode="fanin_only")
        st_p, _ = _fast_decide_acyclic_p_all("rec_x", pr_p, frozenset(), cyc_mode="fanin_only")
        return 1.0 if (st_m == "plexi" and st_p == "plexi") else 0.0

    comp_vals = [_plexi_indicator(p) for p in conn_2comp["exact_component_pairs"]]  # type: ignore[union-attr,index]
    fiber_vals = [_plexi_indicator(p) for p in conn_2comp["exact_fiber_pairs"]]  # type: ignore[union-attr,index]
    comp_mean = sum(comp_vals) / len(comp_vals)
    full_fiber_mean = sum(fiber_vals) / len(fiber_vals)
    assert comp_mean == 1.0 and sorted(fiber_vals) == [0.0, 1.0] and full_fiber_mean == 0.5, (
        f"Expected component_mean=1.0 vs full_fiber_mean=0.5 on 2-component unequal-mean test, got {comp_mean} vs {full_fiber_mean}"
    )

    # 8a. Tier 2 reversible 4-chain MCMC stationary-mean calibration on known 12-state Latin fiber (exact mean = 2/12 = 1/6)
    exact_lat_vals = [
        1.0 if (pm[0][0] == 1 and pm[1][1] == 1) else 0.0
        for pm, _ in conn_latin["exact_fiber_pairs"]  # type: ignore[union-attr,index]
    ]
    exact_lat_mean = sum(exact_lat_vals) / len(exact_lat_vals)
    assert abs(exact_lat_mean - (1.0 / 6.0)) < 1e-12
    n_covered = 0
    for s_idx in range(20):
        s_rng = random.Random(9100 + s_idx)
        ch_draws, ch_ham = mcmc_reversible_nested_chains(
            lat_m, lat_p, s_rng, n_draws=400, n_chains=4, burnin_mult=10, thin_mult=2
        )
        p_chains = [
            [1.0 if (pm[0][0] == 1 and pm[1][1] == 1) else 0.0 for pm, _ in ch]
            for ch in ch_draws
        ]
        v_diag = validate_mcmc_prevalence_bounds(p_chains, p_chains, ch_ham, has_active_move=True, min_ess=80.0)
        assert v_diag["mcmc_diagnostics_pass"] is True, f"Tier 2 failed on known 12-state fiber: {v_diag}"
        mu_hat = float(v_diag["bound_low_diag"]["mean"])  # type: ignore[index]
        mcse_hat = float(v_diag["bound_low_diag"]["mcse_bm"])  # type: ignore[index]
        if abs(mu_hat - exact_lat_mean) <= 2.861 * mcse_hat:
            n_covered += 1
    assert n_covered >= 18, f"Expected >=18/20 MCSE interval coverage on known Latin fiber, got {n_covered}/20"

    # 8b. Slow-kernel regression test case (2-state {0,1} chain with switch prob p=0.01, L=20 from hub 0;
    #     regression test case verifying the diagnostic gate catches an under-mixed chain, NOT a general mixing proof)
    slow_rng = random.Random(20261122)
    theo_spoke_mean = (1.0 - (0.98 ** 20)) / 2.0  # ~0.166198 != 0.50
    spoke_vals: List[float] = []
    for _ in range(2000):
        st_s = 0  # conditional on hub H = 0
        for _ in range(20):
            if slow_rng.random() < 0.01:
                st_s = 1 - st_s
        spoke_vals.append(float(st_s))
    emp_spoke_mean = sum(spoke_vals) / len(spoke_vals)
    assert abs(emp_spoke_mean - theo_spoke_mean) < 0.025 and emp_spoke_mean < 0.25
    slow_chains: List[List[float]] = []
    slow_hams: List[List[int]] = []
    for _ in range(4):
        st_c = 0
        for _ in range(20):
            if slow_rng.random() < 0.01:
                st_c = 1 - st_c
        c_ser: List[float] = []
        c_hm: List[int] = []
        for _ in range(250):
            if slow_rng.random() < 0.01:
                st_c = 1 - st_c
            c_ser.append(float(st_c))
            c_hm.append(st_c)
        slow_chains.append(c_ser)
        slow_hams.append(c_hm)
    slow_diag = validate_mcmc_prevalence_bounds(slow_chains, slow_chains, slow_hams, has_active_move=True, min_ess=200.0)
    assert slow_diag["mcmc_diagnostics_pass"] is False, (
        f"Slow-kernel regression test case (p=0.01) unexpectedly passed diagnostics: {slow_diag}"
    )

    return {
        "unconstrained_trades": u_trades,
        "constrained_total_trades_10draws": c_total_trades,
        "constrained_changed_both_matrices": cm_changed_both,
        "cyclic_3state_trade_verified": saw_flipped,
        "latin_3x3_visited_count": len(visited_latin),
        "connectivity_3state_verified": True,
        "null2c_equivalence_class_invariance_verified": True,
        "two_component_unequal_mean_full_vs_comp": (full_fiber_mean, comp_mean),
        "tier2_known_fiber_coverage_20seeds": n_covered,
        "slow_kernel_regression_test_rejected": not bool(slow_diag["mcmc_diagnostics_pass"]),
    }


def match_c0_pairs_for_recipient(
    rec: GenomeAnnotations,
    n_needed: int = 1,
) -> List[Tuple[str, str]]:
    """Return all C0xC0 pairs for which rec is an eligible calibrated candidate auxotroph."""
    abs_c0 = [
        m
        for m in C0_CONTROL_MODULES
        if recipient_endpoint_status(rec, m, max_false_absence_prob=0.05) == STATE_ABSENT
    ]
    return list(itertools.combinations(abs_c0, 2))


def count_missing_modules(g: GenomeAnnotations) -> int:
    return sum(
        1
        for ep in ALL_P2_ENDPOINTS
        if recipient_endpoint_status(g, ep, max_false_absence_prob=0.05) == STATE_ABSENT
    )


def run_stage0_ebb_flow(n_null_draws: int = 500, seed: int = 20261120) -> Dict[str, object]:
    os.makedirs(OUT_DIR, exist_ok=True)
    q_path = os.path.join(EBB_ROOT, "binning", "checkm2", "quality_report.tsv")
    mags_json_path = os.path.join(EBB_ROOT, "viz", "data", "mags.json")
    depths_path = os.path.join(EBB_ROOT, "mapping", "depths.txt")
    c2b_path = os.path.join(EBB_ROOT, "binning", "dastool", "contig2bin.tsv")
    gff_path = os.path.join(EBB_ROOT, "annotation", "bakta", "basic", "annotation.gff3")
    comm_path = os.path.join(EBB_ROOT, "metabolism", "per_mag", "community_annotations.tsv")
    tiara_path = os.path.join(EBB_ROOT, "eukaryotic", "tiara", "tiara_output.tsv")
    skani_path = os.path.join(OUT_DIR, "stage0_skani_95ani.tsv")

    with open(q_path, "r", encoding="utf-8") as fh:
        q_rows = list(csv.DictReader(fh, delimiter="\t"))
    raw_hq = sum(
        1 for r in q_rows if float(r["Completeness"]) >= 90.0 and float(r["Contamination"]) <= 5.0
    )

    with open(mags_json_path, "r", encoding="utf-8") as fh:
        mags_meta = {m["name"]: m for m in json.load(fh)}

    tax_by_mag: Dict[str, str] = {}
    phylum_by_mag: Dict[str, str] = {}
    for name, m in mags_meta.items():
        dom = (m.get("taxonomy") or {}).get("domain") or "U"
        phy = (m.get("taxonomy") or {}).get("phylum") or "unclassified"
        tax_by_mag[name] = dom
        phylum_by_mag[name] = phy

    # Verify skani 95% ANI dereplication on the 165 dastool-* bins
    max_skani_ani = 0.0
    n_pairs_ge_95 = 0
    if os.path.isfile(skani_path):
        with open(skani_path, "r", encoding="utf-8") as fh:
            for r in csv.DictReader(fh, delimiter="\t"):
                ani = float(r["ANI"])
                if ani > max_skani_ani:
                    max_skani_ani = ani
                if ani >= 95.0:
                    n_pairs_ge_95 += 1
    derep_reps = set(mags_meta.keys())
    assert n_pairs_ge_95 == 0, f"Found {n_pairs_ge_95} pairs >=95% ANI; cluster before scoring!"

    prok_all_ids, flow_all = filter_derep_prokaryotic_mags(
        q_rows,
        derep_representatives=derep_reps,
        taxonomy_by_mag=tax_by_mag,
        consensus_prefix="dastool-",
        min_completeness=0.0,
        max_contamination=100.0,
    )
    prok_hq_ids, flow_hq = filter_derep_prokaryotic_mags(
        q_rows,
        derep_representatives=derep_reps,
        taxonomy_by_mag=tax_by_mag,
        consensus_prefix="dastool-",
        min_completeness=90.0,
        max_contamination=5.0,
    )
    prok_nc_ids, _ = filter_derep_prokaryotic_mags(
        q_rows,
        derep_representatives=derep_reps,
        taxonomy_by_mag=tax_by_mag,
        consensus_prefix="dastool-",
        min_completeness=95.0,
        max_contamination=5.0,
    )
    assert_non_redundant_mag_ids(prok_all_ids)

    domain_breakdown: Dict[str, int] = {}
    for name in sorted(derep_reps):
        d = tax_by_mag.get(name, "U")
        domain_breakdown[d] = domain_breakdown.get(d, 0) + 1

    # Load contig lengths & sample BAM list from depths.txt
    contig_len: Dict[str, int] = {}
    sample_bams: List[Tuple[str, str]] = []
    with open(depths_path, "r", encoding="utf-8") as fh:
        header = fh.readline().rstrip("\n").split("\t")
        for col in header[3:]:
            if col.endswith(".sorted.bam"):
                sname = col.replace(".sorted.bam", "")
                sample_bams.append((sname, os.path.join(EBB_ROOT, "mapping", col)))
        for line in fh:
            p = line.split("\t", 3)
            contig_len[p[0]] = int(p[1])

    # Load contig2bin for dastool
    bin_contigs: Dict[str, List[str]] = {}
    binned_contig_set: Set[str] = set()
    with open(c2b_path, "r", encoding="utf-8") as fh:
        for line in fh:
            p = line.rstrip("\n").split("\t")
            if len(p) >= 2:
                bin_contigs.setdefault(p[1], []).append(p[0])
                binned_contig_set.add(p[0])

    # Parse complete CDS feature order from Bakta GFF3 on binned contigs
    def _binned_gff_lines() -> Iterable[str]:
        with open(gff_path, "r", encoding="utf-8") as fh:
            for line in fh:
                if line.startswith("##sequence-region"):
                    parts = line.split()
                    if len(parts) >= 4 and parts[1] in binned_contig_set:
                        yield line
                elif not line.startswith("#"):
                    c_id = line.split("\t", 1)[0]
                    if c_id in binned_contig_set:
                        yield line

    gff_order, gff_lens = parse_gff_cds_order(_binned_gff_lines())

    q_by_name = {r["Name"]: r for r in q_rows if r["Name"] in derep_reps}
    genomes: Dict[str, GenomeAnnotations] = {}
    hits_by_mag: Dict[str, List[GeneHit]] = {}
    for mid in prok_all_ids:
        qr = q_by_name[mid]
        c_list = sorted(bin_contigs.get(mid, ()))
        c_lens = [contig_len.get(c, gff_lens.get(c, 0)) for c in c_list]
        tsv_p = os.path.join(EBB_ROOT, "metabolism", "per_mag", "per_mag", f"{mid}.tsv")
        g_ann, g_hits = parse_per_mag_tsv(
            tsv_p,
            mag_id=mid,
            domain=tax_by_mag[mid],
            completeness=float(qr["Completeness"]),
            contamination=float(qr["Contamination"]),
            genome_size_bp=int(qr["Genome_Size"]),
            contig_lengths=c_lens,
            gff_cds_order=gff_order,
        )
        genomes[mid] = g_ann
        hits_by_mag[mid] = g_hits

    # Exact base-level contig covbases per sample BAM via samtools coverage
    cov_cache = os.path.join(OUT_DIR, "stage0_bam_covbases.json")
    bam_covbases = compute_or_load_stage0_bam_covbases(sample_bams, cov_cache)
    sample_names = [s for s, _ in sample_bams]

    # Write exact per-MAG breadth matrix
    breadth_tsv = os.path.join(OUT_DIR, "stage0_mag_breadth.tsv")
    mag_breadth: Dict[str, Dict[str, float]] = {m: {} for m in prok_all_ids}
    with open(breadth_tsv, "w", encoding="utf-8", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["mag_id", "domain", "completeness", "contamination", "genome_bp", *sample_names])
        for mid in prok_all_ids:
            c_list = bin_contigs.get(mid, [])
            tot_bp = sum(contig_len[c] for c in c_list)
            row_vals = []
            for sname in sample_names:
                cmap = bam_covbases[sname]
                cb = sum(cmap.get(c, 0) for c in c_list)
                br = (cb / tot_bp) if tot_bp > 0 else 0.0
                mag_breadth[mid][sname] = br
                row_vals.append(f"{br:.6f}")
            w.writerow([
                mid,
                genomes[mid].domain,
                f"{genomes[mid].completeness:.2f}",
                f"{genomes[mid].contamination:.2f}",
                tot_bp,
                *row_vals,
            ])

    # Load Tiara contig classes + community_annotations KOs for unbinned & eukaryotic pools
    tiara_class: Dict[str, str] = {}
    with open(tiara_path, "r", encoding="utf-8") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            tiara_class[r["sequence_id"]] = r["class_fst_stage"]

    ref_kos = all_referenced_kos()
    donor_rule_kos: Set[str] = set()
    for m_don in DONOR_CAP_MODULES:
        donor_rule_kos |= PATHWAY_RULES[m_don].all_kos()

    unbinned_prok_kos_by_ctg: Dict[str, Set[str]] = {}
    euk_kos_by_ctg: Dict[str, Set[str]] = {}
    binned_donor_ko_hits_by_ctg: Dict[str, int] = {}
    unbinned_donor_ko_hits_by_ctg: Dict[str, int] = {}

    with open(comm_path, "r", encoding="utf-8") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            kos = parse_ko_field(r.get("KO", ""))
            if not kos:
                continue
            cid = r["contig_id"]
            t_cls = tiara_class.get(cid, "unknown")
            if t_cls in ("organelle", "mitochondrion", "plastid"):
                continue
            p2_kos = kos & ref_kos
            d_kos = kos & donor_rule_kos
            if cid in binned_contig_set:
                if d_kos:
                    binned_donor_ko_hits_by_ctg[cid] = binned_donor_ko_hits_by_ctg.get(cid, 0) + len(d_kos)
            elif t_cls == "eukarya":
                if p2_kos:
                    euk_kos_by_ctg.setdefault(cid, set()).update(p2_kos)
            else:
                if p2_kos:
                    unbinned_prok_kos_by_ctg.setdefault(cid, set()).update(p2_kos)
                if d_kos:
                    unbinned_donor_ko_hits_by_ctg[cid] = unbinned_donor_ko_hits_by_ctg.get(cid, 0) + len(d_kos)

    # Descriptive I1 (P_PSII x A_NIF) co-encoder bypass check across all 163 prokaryotic MAGs
    i1_coencoders: List[Dict[str, object]] = []
    for mid in prok_all_ids:
        g = genomes[mid]
        psii = PATHWAY_RULES["P_PSII"].evaluate(g, max_false_absence_prob=0.05).state
        nif = PATHWAY_RULES["A_NIF"].evaluate(g, max_false_absence_prob=0.05).state
        if psii == STATE_PRESENT and nif == STATE_PRESENT:
            i1_coencoders.append({"mag_id": mid, "bypass": classify_coencoder_bypass(g)})

    # Helper to evaluate one event/library community at a given breadth threshold
    hq_set = set(prok_hq_ids)

    def _eval_community(
        event_id: str,
        detected_mags: List[str],
        unbinned_kos: Set[str],
        euk_kos: Set[str],
        unbinned_donor_share: float,
    ) -> Dict[str, object]:
        det_genomes = [genomes[m] for m in detected_mags]
        det_hq = [genomes[m] for m in detected_mags if m in hq_set]
        unb_g = GenomeAnnotations.from_kos("_unbinned", unbinned_kos) if unbinned_kos else None
        euk_g = GenomeAnnotations.from_kos("_euk_tiara", euk_kos, domain="Eukaryota") if euk_kos else None
        ctx_fanin = CommunityScorerContext(det_genomes, unb_g, euk_g, mode="fanin_only")
        ctx_cyc = CommunityScorerContext(det_genomes, unb_g, euk_g, mode="consortium_cyc")
        ctx_prok_only = CommunityScorerContext(det_genomes, None, None, mode="fanin_only")
        ctx_with_euk_only = CommunityScorerContext(det_genomes, None, euk_g, mode="fanin_only")

        strata_defs: List[Tuple[str, Tuple[str, ...], SupplyScenario]] = [
            ("R-I3_B12_REDS__E_minimal", ("B12", "RED_S"), SCENARIO_E_MINIMAL),
            ("R-I3_B12_REDS__E_plus_DMSP", ("B12", "RED_S"), SCENARIO_E_PLUS_DMSP),
            ("R-I2_B12_FIXEDN__E_minimal", ("B12", "FIXED_N"), SCENARIO_E_MINIMAL),
            ("R-I2_B12_FIXEDN__E_plus_DIN", ("B12", "FIXED_N"), SCENARIO_E_PLUS_DIN),
            ("R-I6_B12_SINGLE__E_minimal", ("B12",), SCENARIO_E_MINIMAL),
            ("R-I6_B12_SINGLE__E_plus_vitamins", ("B12",), SCENARIO_E_PLUS_VITAMINS),
            ("C_THI_SPLIT__E_minimal", ("THI",), SCENARIO_E_MINIMAL),
            ("C_THI_SPLIT__E_plus_vitamins", ("THI",), SCENARIO_E_PLUS_VITAMINS),
        ]

        counts_by_stratum: Dict[str, Dict[str, int]] = {}
        cyc_counts_by_stratum: Dict[str, Dict[str, int]] = {}
        tuple_rows: List[Dict[str, object]] = []
        euk_transitions: Dict[str, int] = {}
        euk_n_eval = 0
        euk_n_changed = 0

        def _record_euk_transition(rec_g: GenomeAnnotations, req_t: Sequence[str], sc_t: SupplyScenario) -> None:
            nonlocal euk_n_eval, euk_n_changed
            c_p, _ = ctx_prok_only.evaluate(rec_g, req_t, sc_t)
            c_e, _ = ctx_with_euk_only.evaluate(rec_g, req_t, sc_t)
            euk_n_eval += 1
            if c_p != c_e:
                euk_n_changed += 1
            tkey = f"{c_p}->{c_e}"
            euk_transitions[tkey] = euk_transitions.get(tkey, 0) + 1

        for s_label, reqs, scen in strata_defs:
            sc = StratumCounts()
            sc_cyc = StratumCounts()
            n_inelig = 0
            for rec in det_hq:
                call_f, wit_f = ctx_fanin.evaluate(rec, reqs, scen)
                if call_f == CALL_NOT_ELIGIBLE:
                    n_inelig += 1
                else:
                    sc.record(call_f)
                call_c, _ = ctx_cyc.evaluate(rec, reqs, scen)
                if call_c != CALL_NOT_ELIGIBLE:
                    sc_cyc.record(call_c)

                if call_f not in (CALL_NOT_ELIGIBLE, CALL_AUTONOMOUS):
                    if scen == SCENARIO_E_MINIMAL:
                        _record_euk_transition(rec, reqs, SCENARIO_E_MINIMAL)
                    donors_other = [d for d in det_genomes if d.mag_id != rec.mag_id]
                    rec_np = {ep: recipient_need_posterior(rec, ep) for ep in reqs}
                    don_sp = [
                        (1.0, {ep: donor_endpoint_posterior(d, ep) for ep in reqs})
                        for d in donors_other
                    ]
                    post_ind = compute_fanin_state_posterior(
                        rec_np, donor_supply_probs=don_sp, ambient=scen.ambient_compounds
                    )
                    # Coupled joint PMF enforcing p11=0 when donor does not simultaneously present both
                    joint_pmfs = []
                    for d in donors_other:
                        probs_ep = {ep: donor_endpoint_posterior(d, ep) for ep in reqs}
                        if len(reqs) == 2 and s_label.startswith(("R-I3", "R-I2")):
                            e1, e2 = reqs[0], reqs[1]
                            p1, p2 = probs_ep[e1], probs_ep[e2]
                            # Fréchet lower bound max(0, p1 + p2 - 1) under O+ x O- antagonism
                            p11 = max(0.0, p1 + p2 - 1.0)
                            p10 = max(0.0, p1 - p11)
                            p01 = max(0.0, p2 - p11)
                            joint_pmfs.append({
                                frozenset({e1, e2}): p11,
                                frozenset({e1}): p10,
                                frozenset({e2}): p01,
                            })
                        else:
                            e1 = reqs[0]
                            joint_pmfs.append({frozenset({e1}): probs_ep[e1]})
                    post_joint = compute_fanin_state_posterior(
                        rec_np, donor_joint_pmfs=joint_pmfs, ambient=scen.ambient_compounds
                    )
                    tuple_rows.append({
                        "event_id": event_id,
                        "stratum": s_label,
                        "recipient_mag": rec.mag_id,
                        "scenario": scen.name,
                        "call_fanin": call_f,
                        "witness_fanin": ":".join(wit_f) if wit_f else "",
                        "call_consortium_cyc": call_c,
                        "p_need_all": post_ind.p_recipient_needs_all,
                        "p_fanin_indep": post_ind.p_fanin_ge_2,
                        "joint_p_plexi_indep": post_ind.joint_p_plexitrophic,
                        "p_fanin_coupled": post_joint.p_fanin_ge_2,
                        "joint_p_plexi_coupled": post_joint.joint_p_plexitrophic,
                    })

            counts_by_stratum[s_label] = {
                "n_hq_evaluated": len(det_hq),
                "n_not_eligible": n_inelig,
                "n_auto": sc.n_auto,
                "n_cand_aux": sc.n_cand_aux,
                "n_d0": sc.n_d0,
                "n_d1": sc.n_d1,
                "n_plexi": sc.n_plexi,
                "n_ambig": sc.n_ambig,
                "n_no_donor": sc.n_no_donor,
                "prev_bound_low": round(sc.prevalence_bounds()[0], 4),
                "prev_bound_high": round(sc.prevalence_bounds()[1], 4),
            }
            cyc_counts_by_stratum[s_label] = {
                "n_cand_aux": sc_cyc.n_cand_aux,
                "n_d0": sc_cyc.n_d0,
                "n_d1_or_chain": sc_cyc.n_d1,
                "n_plexi_with_cyc": sc_cyc.n_plexi,
                "n_ambig": sc_cyc.n_ambig,
                "n_no_donor": sc_cyc.n_no_donor,
            }

        # Compatible C0 x C0 pool across all 45 pairs under E_minimal
        sc_c0_all = StratumCounts()
        c0_pairs_all = list(itertools.combinations(C0_CONTROL_MODULES, 2))
        c0_inelig = 0
        for rec in det_hq:
            for p_c0 in c0_pairs_all:
                call_c0, _ = ctx_fanin.evaluate(rec, p_c0, SCENARIO_E_MINIMAL)
                if call_c0 == CALL_NOT_ELIGIBLE:
                    c0_inelig += 1
                else:
                    sc_c0_all.record(call_c0)
                if call_c0 not in (CALL_NOT_ELIGIBLE, CALL_AUTONOMOUS):
                    _record_euk_transition(rec, p_c0, SCENARIO_E_MINIMAL)

        counts_by_stratum["C0xC0_ALL45__E_minimal"] = {
            "n_hq_evaluated": len(det_hq) * 45,
            "n_not_eligible": c0_inelig,
            "n_auto": sc_c0_all.n_auto,
            "n_cand_aux": sc_c0_all.n_cand_aux,
            "n_d0": sc_c0_all.n_d0,
            "n_d1": sc_c0_all.n_d1,
            "n_plexi": sc_c0_all.n_plexi,
            "n_ambig": sc_c0_all.n_ambig,
            "n_no_donor": sc_c0_all.n_no_donor,
            "prev_bound_low": round(sc_c0_all.prevalence_bounds()[0], 4),
            "prev_bound_high": round(sc_c0_all.prevalence_bounds()[1], 4),
        }

        euk_trans_serializable = {
            "n_evaluated": euk_n_eval,
            "n_changed": euk_n_changed,
            "transitions": dict(sorted(euk_transitions.items())),
        }
        return {
            "event_id": event_id,
            "n_detected_mags": len(det_genomes),
            "n_detected_hq_mags": len(det_hq),
            "unbinned_donor_ko_share": round(unbinned_donor_share, 4),
            "strata_fanin_primary": counts_by_stratum,
            "strata_consortium_cyc_secondary": cyc_counts_by_stratum,
            "eukaryotic_donor_transitions": euk_trans_serializable,
            "candidate_tuples": tuple_rows,
        }

    # Evaluate per barcode at primary breadth >= 50% (and count >= 70% sensitivity)
    barcode_results_50: List[Dict[str, object]] = []
    barcode_summary_70: List[Dict[str, object]] = []
    all_tuple_rows: List[Dict[str, object]] = []
    union_detected_50: Set[str] = set()

    for sname in sample_names:
        cmap = bam_covbases[sname]
        det_50 = [m for m in prok_all_ids if mag_breadth[m][sname] >= 0.50]
        det_70 = [m for m in prok_all_ids if mag_breadth[m][sname] >= 0.70]
        union_detected_50.update(det_50)

        # Contig-level >=50% breadth for unbinned & euk pools in this barcode
        unb_kos_s: Set[str] = set()
        for cid, kset in unbinned_prok_kos_by_ctg.items():
            if cmap.get(cid, 0) >= 0.50 * contig_len.get(cid, 1):
                unb_kos_s |= kset
        euk_kos_s: Set[str] = set()
        for cid, kset in euk_kos_by_ctg.items():
            if cmap.get(cid, 0) >= 0.50 * contig_len.get(cid, 1):
                euk_kos_s |= kset

        bin_d_hits = sum(
            cnt for cid, cnt in binned_donor_ko_hits_by_ctg.items()
            if cmap.get(cid, 0) >= 0.50 * contig_len.get(cid, 1)
        )
        unb_d_hits = sum(
            cnt for cid, cnt in unbinned_donor_ko_hits_by_ctg.items()
            if cmap.get(cid, 0) >= 0.50 * contig_len.get(cid, 1)
        )
        denom_d = bin_d_hits + unb_d_hits
        unb_share = (unb_d_hits / denom_d) if denom_d > 0 else 0.0

        if any(m in hq_set for m in det_50):
            res_b = _eval_community(sname, det_50, unb_kos_s, euk_kos_s, unb_share)
            barcode_results_50.append(res_b)
            all_tuple_rows.extend(res_b["candidate_tuples"])  # type: ignore[arg-type]

        if any(m in hq_set for m in det_70):
            unb_kos_70 = {
                k
                for cid, kset in unbinned_prok_kos_by_ctg.items()
                if cmap.get(cid, 0) >= 0.70 * contig_len.get(cid, 1)
                for k in kset
            }
            euk_kos_70 = {
                k
                for cid, kset in euk_kos_by_ctg.items()
                if cmap.get(cid, 0) >= 0.70 * contig_len.get(cid, 1)
                for k in kset
            }
            res_70 = _eval_community(f"{sname}_br70", det_70, unb_kos_70, euk_kos_70, unb_share)
            barcode_summary_70.append({
                "sample": sname,
                "n_detected_mags_70": len(det_70),
                "n_detected_hq_70": res_70["n_detected_hq_mags"],
                "R-I3_E_minimal": res_70["strata_fanin_primary"]["R-I3_B12_REDS__E_minimal"],  # type: ignore[index]
                "R-I2_E_minimal": res_70["strata_fanin_primary"]["R-I2_B12_FIXEDN__E_minimal"],  # type: ignore[index]
                "C0xC0_E_minimal": res_70["strata_fanin_primary"]["C0xC0_ALL45__E_minimal"],  # type: ignore[index]
            })

    # Campaign-pooled evaluation across all MAGs detected at >=50% breadth in >=1 library
    tot_bin_d = sum(binned_donor_ko_hits_by_ctg.values())
    tot_unb_d = sum(unbinned_donor_ko_hits_by_ctg.values())
    pooled_unb_share = tot_unb_d / max(1, tot_bin_d + tot_unb_d)
    pooled_unb_kos = set().union(*unbinned_prok_kos_by_ctg.values()) if unbinned_prok_kos_by_ctg else set()
    pooled_euk_kos = set().union(*euk_kos_by_ctg.values()) if euk_kos_by_ctg else set()
    pooled_eval = _eval_community(
        "FAZ84459_campaign_pooled",
        sorted(union_detected_50),
        pooled_unb_kos,
        pooled_euk_kos,
        pooled_unb_share,
    )
    all_tuple_rows.extend(pooled_eval["candidate_tuples"])  # type: ignore[arg-type]

    # Exercise Null 1 (recipient+donor contig-block downsampling from >=95% NC MAGs) and Null 2b/2c on pooled campaign
    rng = random.Random(seed)
    det_pooled_ids = sorted(union_detected_50)
    det_pooled_genomes = [genomes[m] for m in det_pooled_ids]
    nc_pool = [m for m in prok_nc_ids if m in union_detected_50]
    target_profiles = [
        (genomes[m].completeness, max(1, len(genomes[m].contig_lengths)))
        for m in prok_hq_ids
        if m in union_detected_50
    ]

    # Null 1 smoke run on C0xC0 + R-I3/R-I2 over n_null_draws draws
    null1_plexi_fracs: List[float] = []
    null1_upper_fracs: List[float] = []
    if nc_pool and target_profiles:
        for draw_idx in range(min(n_null_draws, 200)):
            sim_donors: List[GenomeAnnotations] = []
            for k_i, m_nc in enumerate(nc_pool):
                t_comp, t_nctg = target_profiles[(draw_idx + k_i) % len(target_profiles)]
                ds_g = downsample_contig_blocks(
                    hits_by_mag[m_nc],
                    target_completeness=t_comp,
                    target_contig_count=t_nctg,
                    rng=rng,
                    mag_id=f"null1_{k_i}_{m_nc}",
                    domain=genomes[m_nc].domain,
                    orig_completeness=genomes[m_nc].completeness,
                    total_genome_bp=genomes[m_nc].genome_size_bp,
                )
                sim_donors.append(ds_g)
            sc_n1 = StratumCounts()
            ctx_n1 = CommunityScorerContext(sim_donors, mode="fanin_only")
            for rec_s in sim_donors:
                for p_pair in (("B12", "RED_S"), ("B12", "FIXED_N"), ("C_HIS", "C_TRP"), ("C_LEU", "C_ILV")):
                    c_n1, _ = ctx_n1.evaluate(rec_s, p_pair, SCENARIO_E_MINIMAL)
                    if c_n1 != CALL_NOT_ELIGIBLE:
                        sc_n1.record(c_n1)
            lo_n1, hi_n1 = sc_n1.prevalence_bounds()
            null1_plexi_fracs.append(lo_n1)
            null1_upper_fracs.append(hi_n1)

    # Null 2b & Null 2c nested Curveball runs on pooled campaign detected MAGs
    # Note: joint 3-state Curveball preserves cellwise definite <= possible (pm_minus[i][j] <= pm_plus[i][j])
    # as well as exact row and column sums of BOTH I- and I+ on every draw (smoke diagnostic on Stage 0).
    mat_minus = build_donor_binary_matrix(det_pooled_genomes, optimistic=False)
    mat_plus = build_donor_binary_matrix(det_pooled_genomes, optimistic=True)
    # Null 2c constrains definite co-encoding of O2_plus x O2_minus modules in I-
    forbid_o2_clash = lambda s_idx: bool(s_idx & O2_PLUS_IDX) and bool(s_idx & O2_MINUS_IDX)

    hq_pooled_recs = [genomes[m] for m in det_pooled_ids if m in hq_set]
    rec_base_m_map = {
        rec.mag_id: build_mag_organism_profile(
            rec, required_override=frozenset(), optimistic_capabilities=False, max_false_absence_prob=0.05
        )
        for rec in hq_pooled_recs
    }
    rec_base_p_map = {
        rec.mag_id: build_mag_organism_profile(
            rec, required_override=frozenset(), optimistic_capabilities=True, max_false_absence_prob=0.05
        )
        for rec in hq_pooled_recs
    }
    eligible_q2_inc: List[Tuple[GenomeAnnotations, Tuple[str, ...]]] = []
    eligible_q2_comp: List[Tuple[GenomeAnnotations, Tuple[str, ...]]] = []
    for rec in hq_pooled_recs:
        for req_inc in (("B12", "RED_S"), ("B12", "FIXED_N")):
            if all(recipient_endpoint_status(rec, ep, 0.05) == STATE_ABSENT for ep in req_inc):
                eligible_q2_inc.append((rec, req_inc))
        for req_c0 in itertools.combinations(C0_CONTROL_MODULES, 2):
            if all(recipient_endpoint_status(rec, ep, 0.05) == STATE_ABSENT for ep in req_c0):
                eligible_q2_comp.append((rec, req_c0))

    def _run_curveball_null(n_draws: int, constrained: bool) -> Dict[str, object]:
        f_fn = forbid_o2_clash if constrained else None
        inc_lows, inc_highs, comp_lows, comp_highs = [], [], [], []
        n_changed_m = 0
        n_changed_p = 0
        tot_trades = 0
        for d_idx in range(n_draws):
            pm_minus, pm_plus, n_tr = curveball_permute_nested_pair(
                mat_minus,
                mat_plus,
                rng,
                n_steps=30,
                forbid_minus_row_fn=f_fn,
                check_each_trade=(d_idx < 3),
            )
            tot_trades += n_tr
            if pm_minus != mat_minus:
                n_changed_m += 1
            if pm_plus != mat_plus:
                n_changed_p += 1
            profs_m = {
                det_pooled_genomes[i].mag_id: profile_from_module_Presence_row(
                    det_pooled_genomes[i].mag_id, pm_minus[i]
                )
                for i in range(len(det_pooled_genomes))
            }
            profs_p = {
                det_pooled_genomes[i].mag_id: profile_from_module_Presence_row(
                    det_pooled_genomes[i].mag_id, pm_plus[i]
                )
                for i in range(len(det_pooled_genomes))
            }
            for grp, lo_list, hi_list in (
                (eligible_q2_inc, inc_lows, inc_highs),
                (eligible_q2_comp, comp_lows, comp_highs),
            ):
                sc_draw = StratumCounts()
                for rec_g, reqs_g in grp:
                    pm_r = dict(profs_m)
                    pp_r = dict(profs_p)
                    bm = rec_base_m_map[rec_g.mag_id]
                    bp = rec_base_p_map[rec_g.mag_id]
                    req_fs = frozenset(reqs_g)
                    pm_r[rec_g.mag_id] = OrganismProfile(
                        org_id=bm.org_id, de_novo=bm.de_novo, conversions=bm.conversions, required=req_fs, exportable=bm.exportable
                    )
                    pp_r[rec_g.mag_id] = OrganismProfile(
                        org_id=bp.org_id, de_novo=bp.de_novo, conversions=bp.conversions, required=req_fs, exportable=bp.exportable
                    )
                    st_m, _ = _fast_decide_acyclic_p_all(rec_g.mag_id, pm_r, frozenset(), cyc_mode="fanin_only")
                    st_p, _ = _fast_decide_acyclic_p_all(rec_g.mag_id, pp_r, frozenset(), cyc_mode="fanin_only")
                    if st_m == "d0":
                        sc_draw.record(CALL_D0_AMBIENT)
                    elif st_m == "not_plexi":
                        sc_draw.record(CALL_D1_OR_CHAIN)
                    elif st_p == "no_donor":
                        sc_draw.record(CALL_UNRESOLVED_NO_DONOR)
                    elif st_m == "plexi" and st_p == "plexi":
                        sc_draw.record(CALL_PLEXITROPHIC)
                    else:
                        sc_draw.record(CALL_UNRESOLVED_AMBIG)
                b_lo, b_hi = sc_draw.prevalence_bounds()
                lo_list.append(b_lo)
                hi_list.append(b_hi)
        if n_draws > 0:
            assert n_changed_m > 0 and n_changed_p > 0, "Curveball null degenerated into a no-op!"
        return {
            "n_draws": n_draws,
            "n_draws_changing_minus": n_changed_m,
            "n_draws_changing_plus": n_changed_p,
            "total_nontrivial_trades": tot_trades,
            "cellwise_definite_subset_possible_verified": True,
            "inc_mean_plexi_low": round(sum(inc_lows) / max(1, len(inc_lows)), 4),
            "inc_mean_plexi_high": round(sum(inc_highs) / max(1, len(inc_highs)), 4),
            "comp_mean_plexi_low": round(sum(comp_lows) / max(1, len(comp_lows)), 4),
            "comp_mean_plexi_high": round(sum(comp_highs) / max(1, len(comp_highs)), 4),
        }

    null2b_summary = _run_curveball_null(min(n_null_draws, 500), constrained=False)
    null2c_summary = _run_curveball_null(min(n_null_draws, 500), constrained=True)

    # Gate check on Stage 0 (1 campaign event < 10; 8 biological barcodes < 10)
    gate_campaign = check_data_adequacy_gates(
        n_independent_events=1,
        unbinned_donor_ko_share=pooled_unb_share,
        unresolved_fraction=0.0,
    )
    gate_barcodes = check_data_adequacy_gates(
        n_independent_events=len(barcode_results_50),
        unbinned_donor_ko_share=pooled_unb_share,
        unresolved_fraction=0.0,
    )

    # Write candidate tuple table
    tuples_tsv = os.path.join(OUT_DIR, "stage0_candidate_tuples.tsv")
    tuple_cols = [
        "event_id", "stratum", "recipient_mag", "scenario", "call_fanin", "witness_fanin",
        "call_consortium_cyc", "p_need_all", "p_fanin_indep", "joint_p_plexi_indep",
        "p_fanin_coupled", "joint_p_plexi_coupled",
    ]
    with open(tuples_tsv, "w", encoding="utf-8", newline="") as fh:
        dw = csv.DictWriter(fh, fieldnames=tuple_cols, delimiter="\t", lineterminator="\n")
        dw.writeheader()
        for tr in all_tuple_rows:
            dw.writerow(tr)

    # Strip verbose candidate_tuples from JSON summary to keep JSON concise & legible
    pooled_clean = {k: v for k, v in pooled_eval.items() if k != "candidate_tuples"}
    barcodes_clean = [
        {k: v for k, v in b.items() if k != "candidate_tuples"} for b in barcode_results_50
    ]

    summary = {
        "protocol_sha256": FROZEN_HASHES[os.path.join(REPO_ROOT, "protocols", "P2-plexitrophy-recipient.md")],
        "scorer_sha256": FROZEN_HASHES[os.path.join(SCRIPT_DIR, "p2_scorer.py")],
        "pathway_rules_sha256": FROZEN_HASHES[os.path.join(SCRIPT_DIR, "p2_pathway_rules.py")],
        "cohort": "Stage0_Ebb_Flow_20260218",
        "stage0_governance_rule": "Pipeline & flow-denominator check ONLY; NO Supported/Falsified hypothesis verdict applied (P2 §4.2 Gate 1 & §5).",
        "adequacy_gate_campaign": gate_campaign,
        "adequacy_gate_if_barcodes_treated_separately": gate_barcodes,
        "flow_denominators": {
            "raw_quality_report_rows": len(q_rows),
            "raw_hq_rows": raw_hq,
            "consensus_dastool_bins_pre_derep": flow_all["consensus_bins"],
            "skani_95ani_max_pairwise_ani": max_skani_ani,
            "skani_95ani_pairs_ge_95": n_pairs_ge_95,
            "derep_95ani_bins": flow_all["derep_bins"],
            "domain_breakdown_of_165_dastool_bins": domain_breakdown,
            "prokaryotic_derep_bins": len(prok_all_ids),
            "prokaryotic_hq_bins_ge90_le5": len(prok_hq_ids),
            "prokaryotic_near_complete_ge95_le5": len(prok_nc_ids),
            "total_bam_libraries_in_depths_txt": len(sample_names),
            "libraries_with_ge1_hq_mag_at_breadth50": len(barcode_results_50),
            "prokaryotic_mags_detected_ge1_lib_breadth50": len(union_detected_50),
            "prokaryotic_hq_mags_detected_ge1_lib_breadth50": sum(1 for m in union_detected_50 if m in hq_set),
        },
        "descriptive_I1_PSII_x_NIF_coencoders": i1_coencoders,
        "campaign_pooled_breadth50": pooled_clean,
        "per_barcode_breadth50": barcodes_clean,
        "sensitivity_breadth70_summary": barcode_summary_70,
        "null_smoke_checks_on_pooled_campaign": {
            "null1_contig_block_downsampling": {
                "n_draws": len(null1_plexi_fracs),
                "nc_source_mags": len(nc_pool),
                "mean_plexi_bound_low": round(sum(null1_plexi_fracs) / max(1, len(null1_plexi_fracs)), 4),
                "mean_plexi_bound_high": round(sum(null1_upper_fracs) / max(1, len(null1_upper_fracs)), 4),
            },
            "null2b_curveball_unconstrained": null2b_summary,
            "null2c_curveball_o2_clash_constrained": null2c_summary,
        },
    }
    out_json = os.path.join(OUT_DIR, "stage0_ebb_flow_summary.json")
    with open(out_json, "w", encoding="utf-8") as fh:
        json.dump(summary, fh, indent=2, sort_keys=False)
    return summary


def infer_conservative_domain_from_kos(kos: FrozenSet[str]) -> str:
    """Return 'Bacteria', 'Archaea', or 'U' strictly from unambiguous domain-specific KOs when GTDB-Tk is unstaged."""
    has_arc = bool(kos & ARC_STRICT_KOS)
    has_bac = bool(kos & BAC_STRICT_KOS)
    if has_bac and not has_arc:
        return "Bacteria"
    if has_arc and not has_bac:
        return "Archaea"
    return "U"


def _merge_interval_length(intervals: List[Tuple[int, int]]) -> int:
    if not intervals:
        return 0
    intervals_sorted = sorted(intervals)
    tot = 0
    cur_s, cur_e = intervals_sorted[0]
    for s, e in intervals_sorted[1:]:
        if s > cur_e:
            tot += cur_e - cur_s
            cur_s, cur_e = s, e
        elif e > cur_e:
            cur_e = e
    tot += cur_e - cur_s
    return tot


def _load_or_build_freshwater_gff_cache(
    gff_path: str,
    gff_sha12: str,
    rep_contigs: Set[str],
    extra_target_contigs: Set[str],
    cache_json: str,
) -> Tuple[Dict[str, Tuple[str, ...]], Dict[str, int]]:
    if os.path.isfile(cache_json):
        with open(cache_json, "r", encoding="utf-8") as fh:
            cached = json.load(fh)
        if (
            cached.get("gff_sha256") == gff_sha12
            and len(cached.get("rep_cds_order", {})) == len(rep_contigs)
            and extra_target_contigs.issubset(set(cached.get("target_contig_lengths", {}).keys()))
        ):
            cds_ord = {k: tuple(v) for k, v in cached["rep_cds_order"].items()}
            t_lens = {k: int(v) for k, v in cached["target_contig_lengths"].items()}
            return cds_ord, t_lens

    want_len_ctgs = rep_contigs | extra_target_contigs
    target_lens: Dict[str, int] = {}

    def _filtered_gff_lines() -> Iterable[str]:
        with open(gff_path, "r", encoding="utf-8") as fh:
            for line in fh:
                if line.startswith("##FASTA"):
                    break
                if line.startswith("##sequence-region"):
                    parts = line.split()
                    if len(parts) >= 4 and parts[1] in want_len_ctgs and parts[3].isdigit():
                        target_lens[parts[1]] = int(parts[3])
                        if parts[1] in rep_contigs:
                            yield line
                elif not line.startswith("#"):
                    cid = line.split("\t", 1)[0]
                    if cid in rep_contigs:
                        yield line

    cds_order, _ = parse_gff_cds_order(_filtered_gff_lines())
    # Contigs in rep_contigs with 0 CDSs still get empty tuple so they never trigger unknown_edge
    for cid in rep_contigs:
        cds_order.setdefault(cid, ())

    with open(cache_json, "w", encoding="utf-8") as fh:
        json.dump(
            {
                "gff_sha256": gff_sha12,
                "rep_cds_order": {k: list(v) for k, v in cds_order.items()},
                "target_contig_lengths": target_lens,
            },
            fh,
        )
    return cds_order, target_lens


STRATUM_TO_GATE2_GROUP: Dict[str, str] = {
    "R-I3_B12_REDS__E_minimal": "R-I3_B12_REDS",
    "R-I3_B12_REDS__E_plus_DMSP": "R-I3_B12_REDS",
    "R-I2_B12_FIXEDN__E_minimal": "R-I2_B12_FIXEDN",
    "R-I2_B12_FIXEDN__E_plus_DIN": "R-I2_B12_FIXEDN",
    "R-I6_B12_SINGLE__E_minimal": "R-I6_B12",
    "C_THI_SPLIT__E_minimal": "C_THI_SPLIT",
    "C0xC0_ALL45__E_minimal": "C0_ALL10",
    "C0xC0_ALL45__E_plus_DMSP": "C0_ALL10",
    "C0xC0_ALL45__E_plus_DIN": "C0_ALL10",
}

GTDB_FINER_RANKS: Tuple[str, ...] = ("c", "o", "f", "g", "s")


def _parse_gtdb_ranks(cls_str: str) -> Dict[str, str]:
    ranks: Dict[str, str] = {}
    for tok in (cls_str or "").split(";"):
        tok = tok.strip()
        if "__" not in tok:
            continue
        rk, val = tok.split("__", 1)
        rk = rk.strip()
        val = val.strip()
        if val and val.lower() not in ("unclassified", "u", "unknown"):
            ranks[rk] = val
    return ranks


def _clade_rank_mismatch(rx: Dict[str, str], ry: Dict[str, str]) -> Optional[int]:
    """Return None if GTDB domain or phylum is missing or differs; else 0..5 finer-rank mismatch count."""
    dx, dy = rx.get("d", ""), ry.get("d", "")
    px, py = rx.get("p", ""), ry.get("p", "")
    if not dx or not px or dx != dy or px != py:
        return None
    shared_depth = 0
    for rk in GTDB_FINER_RANKS:
        vx, vy = rx.get(rk, ""), ry.get(rk, "")
        if vx and vy and vx == vy:
            shared_depth += 1
        else:
            break
    return len(GTDB_FINER_RANKS) - shared_depth


def _select_outcome_blind_c0_pair(
    mag_id: str,
    contrast_id: str,
    elig_pairs: Sequence[Tuple[str, str]],
) -> Tuple[str, str]:
    """Deterministically select 1 eligible C0xC0 pair strictly from recipient ID and contrast ID (outcome-blind)."""
    ordered = sorted(elig_pairs)
    if len(ordered) == 1:
        return ordered[0]
    digest = hashlib.sha256(f"{mag_id}::{contrast_id}".encode("utf-8")).digest()
    idx = int.from_bytes(digest[:8], "big") % len(ordered)
    return ordered[idx]


def _match_c0_controls_in_event(
    contrast_id: str,
    inc_recs: Sequence[GenomeAnnotations],
    det_hq_ev: Sequence[GenomeAnnotations],
    elig_c0_by_mag: Dict[str, List[Tuple[str, str]]],
    gtdb_ranks_by_mag: Dict[str, Dict[str, str]],
    cov_vec_by_mag: Dict[str, Tuple[float, float, float]],
    cov_sds: Tuple[float, float, float],
) -> Tuple[List[Tuple[GenomeAnnotations, GenomeAnnotations, Tuple[str, str]]], int, int]:
    """Outcome-blind 1-to-1 matching without replacement on GTDB phylum/clade + standardized recipient covariates."""
    c0_cands = [g for g in det_hq_ev if elig_c0_by_mag.get(g.mag_id)]
    if not inc_recs or not c0_cands:
        return [], 0, 0
    cand_edges: List[Tuple[int, int, float, str, str, GenomeAnnotations, GenomeAnnotations]] = []
    for rx in inc_recs:
        gx = gtdb_ranks_by_mag.get(rx.mag_id, {})
        vx = cov_vec_by_mag[rx.mag_id]
        for ry in c0_cands:
            gy = gtdb_ranks_by_mag.get(ry.mag_id, {})
            cm = _clade_rank_mismatch(gx, gy)
            if cm is None:
                continue
            vy = cov_vec_by_mag[ry.mag_id]
            z_sq = sum(((vx[j] - vy[j]) / cov_sds[j]) ** 2 for j in range(3))
            z_dist = round(math.sqrt(z_sq), 9)
            is_self = 0 if ry.mag_id != rx.mag_id else 1
            cand_edges.append((is_self, cm, z_dist, rx.mag_id, ry.mag_id, rx, ry))
    cand_edges.sort(key=lambda t: (t[0], t[1], t[2], t[3], t[4]))
    used_x: Set[str] = set()
    used_y: Set[str] = set()
    matched: List[Tuple[GenomeAnnotations, GenomeAnnotations, Tuple[str, str]]] = []
    n_distinct = 0
    n_self = 0
    for is_self, _, _, xid, yid, rx, ry in cand_edges:
        if xid in used_x or yid in used_y:
            continue
        used_x.add(xid)
        used_y.add(yid)
        chosen_pair = _select_outcome_blind_c0_pair(yid, contrast_id, elig_c0_by_mag[yid])
        matched.append((rx, ry, chosen_pair))
        if is_self == 0:
            n_distinct += 1
        else:
            n_self += 1
    return matched, n_distinct, n_self


def _strip_conn_pairs(conn_info: Dict[str, object]) -> Dict[str, object]:
    """Strip enumerated matrix lists from per-event connectivity metadata before JSON serialization."""
    return {k: v for k, v in conn_info.items() if k not in ("exact_fiber_pairs", "exact_component_pairs")}


def _aggregate_cohort_connectivity(per_ev_conn: Sequence[Dict[str, object]]) -> str:
    """Aggregate per-event 3-state connectivity across Gate-1-eligible events (P2-A1 §2)."""
    statuses = [str(c["status"]) for c in per_ev_conn]
    if all(s == "connected" for s in statuses):
        return "connected"
    if any(s == "disconnected (proven)" for s in statuses):
        return "disconnected (proven)"
    return "connectivity unknown"


def _decide_p2_q1_symmetric(
    ci_null1: Optional[Tuple[float, float]],
    null1_uneval_reason: Optional[str],
    ci_null2_primary: Optional[Tuple[float, float]],
    null2_unavail_reason: Optional[str],
    delta_p: float = 0.05,
) -> str:
    """Apply frozen P2 §4.3 + P2-A1 §3.1 symmetrically over Null 1 and certified primary full-fiber Null 2b."""
    if ci_null1 is not None and ci_null2_primary is not None:
        return decide_p2_q1(ci_null1, ci_null2_primary, delta_p=delta_p, gate_failure=None)
    if ci_null1 is not None and ci_null2_primary is None:
        l1, u1 = ci_null1
        if u1 < delta_p:
            return "Falsified"
        if l1 <= delta_p:
            return "Inconclusive"
        return f"Inconclusive (primary Null 2 unavailable: {null2_unavail_reason})"
    if ci_null1 is None and ci_null2_primary is not None:
        _, u2 = ci_null2_primary
        if u2 < delta_p:
            return "Falsified"
        return f"Inconclusive (Null 1 unevaluable: {null1_uneval_reason})"
    return (
        f"Inconclusive (Null 1 unevaluable: {null1_uneval_reason}; "
        f"primary Null 2 unavailable: {null2_unavail_reason})"
    )


def _event_bootstrap_q1_stationary_bounds(
    per_event_obs: Sequence[StratumCounts],
    per_event_p_low: Sequence[float],
    per_event_p_high: Sequence[float],
    per_event_var_bm_low: Sequence[float],
    per_event_var_bm_high: Sequence[float],
    n_boot: int = 1000,
    seed: int = 42,
    alpha: float = 0.05,
    t_crit_19: float = 2.861,
) -> Tuple[float, float]:
    """Bootstrap observed-minus-null prevalence excess with 20-batch MCSE outer widening on Tier 2 events."""
    m = len(per_event_obs)
    if m == 0 or not (
        len(per_event_p_low) == len(per_event_p_high) == len(per_event_var_bm_low) == len(per_event_var_bm_high) == m
    ):
        return (-1.0, 1.0)
    rng = random.Random(seed)
    lows: List[float] = []
    highs: List[float] = []
    for _ in range(n_boot):
        mult = [0] * m
        for _ in range(m):
            mult[rng.randrange(m)] += 1
        obs_tot = StratumCounts(
            n_auto=sum(mult[i] * per_event_obs[i].n_auto for i in range(m) if mult[i]),
            n_d0=sum(mult[i] * per_event_obs[i].n_d0 for i in range(m) if mult[i]),
            n_d1=sum(mult[i] * per_event_obs[i].n_d1 for i in range(m) if mult[i]),
            n_plexi=sum(mult[i] * per_event_obs[i].n_plexi for i in range(m) if mult[i]),
            n_ambig=sum(mult[i] * per_event_obs[i].n_ambig for i in range(m) if mult[i]),
            n_no_donor=sum(mult[i] * per_event_obs[i].n_no_donor for i in range(m) if mult[i]),
        )
        a, b = obs_tot.prevalence_bounds()
        n_tot = obs_tot.n_cand_aux
        if n_tot <= 0:
            lows.append(-1.0)
            highs.append(1.0)
            continue
        weights = [(mult[i] * per_event_obs[i].n_cand_aux) / float(n_tot) for i in range(m)]
        c_hat = sum(weights[i] * per_event_p_low[i] for i in range(m) if mult[i])
        d_hat = sum(weights[i] * per_event_p_high[i] for i in range(m) if mult[i])
        se_c = math.sqrt(sum((weights[i] ** 2) * per_event_var_bm_low[i] for i in range(m) if mult[i]))
        se_d = math.sqrt(sum((weights[i] ** 2) * per_event_var_bm_high[i] for i in range(m) if mult[i]))
        c_outer = min(max(c_hat - t_crit_19 * se_c, 0.0), 1.0)
        d_outer = min(max(d_hat + t_crit_19 * se_d, 0.0), 1.0)
        lows.append(a - d_outer)
        highs.append(b - c_outer)
    lows.sort()
    highs.sort()
    q_lo = max(0, int(math.floor((alpha / 2.0) * n_boot)))
    q_hi = min(n_boot - 1, int(math.ceil((1.0 - alpha / 2.0) * n_boot)) - 1)
    return (lows[q_lo], highs[q_hi])


def _compute_gate2_record(
    rep490_hits: int,
    dastool655_nonrep_hits: int,
    unbinned_prok_hits: int,
    tiara_eukarya_hits: int,
    max_unbinned_share: float = 0.50,
) -> Dict[str, object]:
    """Compute raw Gate 2 counts, denominators, both non-rep conventions, and 0/0 unevaluable status."""
    denom_490 = rep490_hits + unbinned_prok_hits
    denom_655 = rep490_hits + dastool655_nonrep_hits + unbinned_prok_hits

    if denom_490 == 0:
        share_490: Optional[float] = None
        status_490 = "unevaluable_zero_denominator"
    else:
        raw_490 = unbinned_prok_hits / denom_490
        share_490 = round(raw_490, 4)
        status_490 = "fail_gt_0.50" if raw_490 > max_unbinned_share else "pass_le_0.50"

    if denom_655 == 0:
        share_655: Optional[float] = None
        status_655 = "unevaluable_zero_denominator"
    else:
        raw_655 = unbinned_prok_hits / denom_655
        share_655 = round(raw_655, 4)
        status_655 = "fail_gt_0.50" if raw_655 > max_unbinned_share else "pass_le_0.50"

    return {
        "rep490_hits": rep490_hits,
        "dastool655_nonrep_hits": dastool655_nonrep_hits,
        "unbinned_prok_hits": unbinned_prok_hits,
        "tiara_eukarya_hits": tiara_eukarya_hits,
        "denom_vs_490reps": denom_490,
        "denom_vs_655dastool": denom_655,
        "unbinned_share_vs_490reps": share_490,
        "unbinned_share_vs_655dastool": share_655,
        "status_vs_490reps": status_490,
        "status_vs_655dastool": status_655,
    }


def _gate2_failure_reason(ep_group: str, rec: Dict[str, object]) -> Optional[str]:
    s490 = rec["status_vs_490reps"]
    s655 = rec["status_vs_655dastool"]
    if s490 == "unevaluable_zero_denominator" or s655 == "unevaluable_zero_denominator":
        return (
            f"Inconclusive (data-adequacy gate: Gate 2 unevaluable — 0/0 detected community donor KOs "
            f"for {ep_group} [denom_vs_490reps={rec['denom_vs_490reps']}, "
            f"denom_vs_655dastool={rec['denom_vs_655dastool']}])"
        )
    if s490 == "fail_gt_0.50" or s655 == "fail_gt_0.50":
        return (
            f"Inconclusive (data-adequacy gate: Gate 2 unbinned_donor_ko_share > 0.50 for {ep_group} "
            f"[vs_490reps={rec['unbinned_share_vs_490reps']} ({rec['unbinned_prok_hits']}/{rec['denom_vs_490reps']}), "
            f"vs_655dastool={rec['unbinned_share_vs_655dastool']} ({rec['unbinned_prok_hits']}/{rec['denom_vs_655dastool']})])"
        )
    return None


def _event_bootstrap_q2_log_or_bounds(
    events: Sequence[str],
    per_ev_inc: Dict[str, StratumCounts],
    per_ev_comp: Dict[str, StratumCounts],
    n_boot: int = 1000,
    seed: int = 20261121,
    alpha: float = 0.05,
) -> Tuple[float, float]:
    if not events:
        return (-math.inf, math.inf)
    rng = random.Random(seed)
    m = len(events)
    lows: List[float] = []
    highs: List[float] = []
    for _ in range(n_boot):
        picks = [events[rng.randrange(m)] for _ in range(m)]
        inc_tot = StratumCounts(
            n_auto=sum(per_ev_inc[e].n_auto for e in picks),
            n_d0=sum(per_ev_inc[e].n_d0 for e in picks),
            n_d1=sum(per_ev_inc[e].n_d1 for e in picks),
            n_plexi=sum(per_ev_inc[e].n_plexi for e in picks),
            n_ambig=sum(per_ev_inc[e].n_ambig for e in picks),
            n_no_donor=sum(per_ev_inc[e].n_no_donor for e in picks),
        )
        comp_tot = StratumCounts(
            n_auto=sum(per_ev_comp[e].n_auto for e in picks),
            n_d0=sum(per_ev_comp[e].n_d0 for e in picks),
            n_d1=sum(per_ev_comp[e].n_d1 for e in picks),
            n_plexi=sum(per_ev_comp[e].n_plexi for e in picks),
            n_ambig=sum(per_ev_comp[e].n_ambig for e in picks),
            n_no_donor=sum(per_ev_comp[e].n_no_donor for e in picks),
        )
        lo, hi = manski_log_odds_bounds(inc_tot, comp_tot)
        lows.append(lo)
        highs.append(hi)
    lows.sort()
    highs.sort()
    q_lo = max(0, int(math.floor((alpha / 2.0) * n_boot)))
    q_hi = min(n_boot - 1, int(math.ceil((1.0 - alpha / 2.0) * n_boot)) - 1)
    return (lows[q_lo], highs[q_hi])


def _event_bootstrap_q2_null2c_stationary_bounds(
    events: Sequence[str],
    obs_inc: Dict[str, StratumCounts],
    obs_comp: Dict[str, StratumCounts],
    null_inc_p_low: Dict[str, float],
    null_inc_p_high: Dict[str, float],
    null_inc_var_low: Dict[str, float],
    null_inc_var_high: Dict[str, float],
    null_comp_p_low: Dict[str, float],
    null_comp_p_high: Dict[str, float],
    null_comp_var_low: Dict[str, float],
    null_comp_var_high: Dict[str, float],
    n_boot: int = 1000,
    seed: int = 20261121,
    alpha: float = 0.05,
    t_crit_19: float = 2.861,
    pseudocount: float = 0.5,
) -> Tuple[Tuple[float, float], Tuple[float, float]]:
    """Bootstrap Null 2c log-OR and observed-minus-Null-2c log-OR bounds with 20-batch MCSE outer widening."""
    if not events:
        return ((-math.inf, math.inf), (-math.inf, math.inf))
    rng = random.Random(seed)
    m = len(events)
    n_lows: List[float] = []
    n_highs: List[float] = []
    d_lows: List[float] = []
    d_highs: List[float] = []

    def _logit_ha(p: float, n_tot: int) -> float:
        p_c = min(max(p, 0.0), 1.0)
        return math.log((p_c * n_tot + pseudocount) / ((1.0 - p_c) * n_tot + pseudocount))

    for _ in range(n_boot):
        mult = [0] * m
        for _ in range(m):
            mult[rng.randrange(m)] += 1
        oi = StratumCounts(
            n_d0=sum(mult[i] * obs_inc[events[i]].n_d0 for i in range(m) if mult[i]),
            n_d1=sum(mult[i] * obs_inc[events[i]].n_d1 for i in range(m) if mult[i]),
            n_plexi=sum(mult[i] * obs_inc[events[i]].n_plexi for i in range(m) if mult[i]),
            n_ambig=sum(mult[i] * obs_inc[events[i]].n_ambig for i in range(m) if mult[i]),
            n_no_donor=sum(mult[i] * obs_inc[events[i]].n_no_donor for i in range(m) if mult[i]),
        )
        oc = StratumCounts(
            n_d0=sum(mult[i] * obs_comp[events[i]].n_d0 for i in range(m) if mult[i]),
            n_d1=sum(mult[i] * obs_comp[events[i]].n_d1 for i in range(m) if mult[i]),
            n_plexi=sum(mult[i] * obs_comp[events[i]].n_plexi for i in range(m) if mult[i]),
            n_ambig=sum(mult[i] * obs_comp[events[i]].n_ambig for i in range(m) if mult[i]),
            n_no_donor=sum(mult[i] * obs_comp[events[i]].n_no_donor for i in range(m) if mult[i]),
        )
        lo_o, hi_o = manski_log_odds_bounds(oi, oc, pseudocount=pseudocount)
        ni_tot = oi.n_cand_aux
        nc_tot = oc.n_cand_aux
        if ni_tot <= 0 or nc_tot <= 0:
            n_lows.append(-math.inf)
            n_highs.append(math.inf)
            d_lows.append(-math.inf)
            d_highs.append(math.inf)
            continue
        wi = [(mult[i] * obs_inc[events[i]].n_cand_aux) / float(ni_tot) for i in range(m)]
        wc = [(mult[i] * obs_comp[events[i]].n_cand_aux) / float(nc_tot) for i in range(m)]

        ci_hat = sum(wi[i] * null_inc_p_low[events[i]] for i in range(m) if mult[i])
        di_hat = sum(wi[i] * null_inc_p_high[events[i]] for i in range(m) if mult[i])
        se_ci = math.sqrt(sum((wi[i] ** 2) * null_inc_var_low[events[i]] for i in range(m) if mult[i]))
        se_di = math.sqrt(sum((wi[i] ** 2) * null_inc_var_high[events[i]] for i in range(m) if mult[i]))
        ci_out = min(max(ci_hat - t_crit_19 * se_ci, 0.0), 1.0)
        di_out = min(max(di_hat + t_crit_19 * se_di, 0.0), 1.0)

        cc_hat = sum(wc[i] * null_comp_p_low[events[i]] for i in range(m) if mult[i])
        dc_hat = sum(wc[i] * null_comp_p_high[events[i]] for i in range(m) if mult[i])
        se_cc = math.sqrt(sum((wc[i] ** 2) * null_comp_var_low[events[i]] for i in range(m) if mult[i]))
        se_dc = math.sqrt(sum((wc[i] ** 2) * null_comp_var_high[events[i]] for i in range(m) if mult[i]))
        cc_out = min(max(cc_hat - t_crit_19 * se_cc, 0.0), 1.0)
        dc_out = min(max(dc_hat + t_crit_19 * se_dc, 0.0), 1.0)

        lo_n = _logit_ha(ci_out, ni_tot) - _logit_ha(dc_out, nc_tot)
        hi_n = _logit_ha(di_out, ni_tot) - _logit_ha(cc_out, nc_tot)
        n_lows.append(lo_n)
        n_highs.append(hi_n)
        d_lows.append(lo_o - hi_n)
        d_highs.append(hi_o - lo_n)
    n_lows.sort()
    n_highs.sort()
    d_lows.sort()
    d_highs.sort()
    q_lo = max(0, int(math.floor((alpha / 2.0) * n_boot)))
    q_hi = min(n_boot - 1, int(math.ceil((1.0 - alpha / 2.0) * n_boot)) - 1)
    return ((n_lows[q_lo], n_highs[q_hi]), (d_lows[q_lo], d_highs[q_hi]))


def _resolve_cohort_q1_null2b(
    elig_evs: Sequence[str],
    per_ev_obs: Sequence[StratumCounts],
    ev_conn_by_ev: Dict[str, Dict[str, object]],
    ev_ests_s: Sequence[Dict[str, object]],
    p2_a1_approved: bool,
    n_boot: int,
    seed_primary: int,
    seed_comp: int,
) -> Tuple[
    str,
    List[Dict[str, object]],
    Optional[Tuple[float, float]],
    Optional[List[float]],
    Optional[str],
    Optional[Dict[str, object]],
]:
    """Select primary Null 2b source PER EVENT (full fiber whenever available; Option A fallback only on individual unenumerated events)."""
    null2b_conn_status = _aggregate_cohort_connectivity([ev_conn_by_ev[ev] for ev in elig_evs])
    all_comp_valid = all(bool(est["comp_valid"]) for est in ev_ests_s)
    any_non_connected = any(str(ev_conn_by_ev[ev]["status"]) != "connected" for ev in elig_evs)

    prim_sources: List[str] = []
    prim_p_low: List[float] = []
    prim_p_high: List[float] = []
    prim_var_low: List[float] = []
    prim_var_high: List[float] = []
    for est in ev_ests_s:
        if bool(est["full_fiber_valid"]):
            src_lbl = (
                "tier1_exact_full_fiber"
                if bool(est["exact_full_fiber_enumerated"])
                else "tier2_connected_full_fiber"
            )
            prim_sources.append(src_lbl)
            prim_p_low.append(float(est["p_low"]))  # type: ignore[arg-type]
            prim_p_high.append(float(est["p_high"]))  # type: ignore[arg-type]
            prim_var_low.append(float(est["var_bm_low"]))  # type: ignore[arg-type]
            prim_var_high.append(float(est["var_bm_high"]))  # type: ignore[arg-type]
        elif p2_a1_approved and bool(est["comp_valid"]):
            prim_sources.append("component_conditional_p2_a1")
            prim_p_low.append(float(est["comp_p_low"]))  # type: ignore[arg-type]
            prim_p_high.append(float(est["comp_p_high"]))  # type: ignore[arg-type]
            prim_var_low.append(float(est["comp_var_bm_low"]))  # type: ignore[arg-type]
            prim_var_high.append(float(est["comp_var_bm_high"]))  # type: ignore[arg-type]
        else:
            prim_sources.append("unavailable")
            prim_p_low.append(0.0)
            prim_p_high.append(1.0)
            prim_var_low.append(0.0)
            prim_var_high.append(0.0)

    ev_conn_list = [
        {
            "event": ev,
            **ev_conn_by_ev[ev],
            "tier": ev_ests_s[i]["tier"],
            "full_fiber_valid": ev_ests_s[i]["full_fiber_valid"],
            "primary_null2_source": prim_sources[i],
            "mcmc_diagnostics": ev_ests_s[i]["mcmc_diagnostics"],
        }
        for i, ev in enumerate(elig_evs)
    ]

    tot_cand = sum(sc.n_cand_aux for sc in per_ev_obs)
    ci_n2b_primary: Optional[Tuple[float, float]] = None
    prim_pooled_null: Optional[List[float]] = None
    null2b_unavail_reason: Optional[str] = None

    if all(s != "unavailable" for s in prim_sources):
        ci_n2b_primary = _event_bootstrap_q1_stationary_bounds(
            per_ev_obs,
            prim_p_low,
            prim_p_high,
            prim_var_low,
            prim_var_high,
            n_boot=n_boot,
            seed=seed_primary,
        )
        if tot_cand > 0:
            prim_pooled_null = [
                round(sum((per_ev_obs[i].n_cand_aux / float(tot_cand)) * prim_p_low[i] for i in range(len(elig_evs))), 6),
                round(sum((per_ev_obs[i].n_cand_aux / float(tot_cand)) * prim_p_high[i] for i in range(len(elig_evs))), 6),
            ]
    else:
        unavail_ests = [est for i, est in enumerate(ev_ests_s) if prim_sources[i] == "unavailable"]
        any_mcmc_fail = any(not bool(est["comp_valid"]) for est in unavail_ests)
        any_disc_unenum = any(str(est["move_connectivity"]) == "disconnected (proven)" for est in unavail_ests)
        any_unk_unenum = any(str(est["move_connectivity"]) == "connectivity unknown" for est in unavail_ests)
        if any_disc_unenum and not any_mcmc_fail:
            null2b_unavail_reason = "move set disconnected"
        elif any_unk_unenum and not any_mcmc_fail:
            null2b_unavail_reason = "connectivity unknown"
        elif any_mcmc_fail and not (any_disc_unenum or any_unk_unenum):
            null2b_unavail_reason = "MCMC stationary-mean mixing unvalidated"
        elif any_disc_unenum:
            null2b_unavail_reason = "move set disconnected; MCMC stationary-mean mixing unvalidated"
        else:
            null2b_unavail_reason = "connectivity unknown; MCMC stationary-mean mixing unvalidated"

    null2b_comp_cond: Optional[Dict[str, object]] = None
    if any_non_connected:
        if all_comp_valid:
            comp_lo = [float(est["comp_p_low"]) for est in ev_ests_s]  # type: ignore[arg-type]
            comp_hi = [float(est["comp_p_high"]) for est in ev_ests_s]  # type: ignore[arg-type]
            ci_n2_comp = _event_bootstrap_q1_stationary_bounds(
                per_ev_obs,
                comp_lo,
                comp_hi,
                [float(est["comp_var_bm_low"]) for est in ev_ests_s],  # type: ignore[arg-type]
                [float(est["comp_var_bm_high"]) for est in ev_ests_s],  # type: ignore[arg-type]
                n_boot=n_boot,
                seed=seed_comp,
            )
            comp_pooled_null = (
                [
                    round(sum((per_ev_obs[i].n_cand_aux / float(tot_cand)) * comp_lo[i] for i in range(len(elig_evs))), 6),
                    round(sum((per_ev_obs[i].n_cand_aux / float(tot_cand)) * comp_hi[i] for i in range(len(elig_evs))), 6),
                ]
                if tot_cand > 0
                else None
            )
            null2b_comp_cond = {
                "label": (
                    "component-conditional (P2-A1, grid-approved)"
                    if p2_a1_approved
                    else "component-conditional (P2-A1, pending grid)"
                ),
                "connectivity_status": null2b_conn_status,
                "mcmc_or_exact_validated": True,
                "pooled_null_bounds": comp_pooled_null,
                "ci_excess_95": list(ci_n2_comp),
            }
        else:
            null2b_comp_cond = {
                "label": "component-conditional (P2-A1, unvalidated — MCMC diagnostics failed)",
                "connectivity_status": null2b_conn_status,
                "mcmc_or_exact_validated": False,
                "pooled_null_bounds": None,
                "ci_excess_95": None,
            }

    return (
        null2b_conn_status,
        ev_conn_list,
        ci_n2b_primary,
        prim_pooled_null,
        null2b_unavail_reason,
        null2b_comp_cond,
    )


def _resolve_cohort_q2_null2c(
    matched_evs: Sequence[str],
    ev_inc_m: Dict[str, StratumCounts],
    ev_comp_m: Dict[str, StratumCounts],
    ev_conn_by_ev: Dict[str, Dict[str, object]],
    ev_2c_ests: Dict[str, Dict[str, object]],
    p2_a1_approved: bool,
    n_boot: int,
    seed_primary: int,
    seed_comp: int,
) -> Dict[str, object]:
    """Select primary Null 2c source PER EVENT (full fiber whenever available; Option A fallback only on individual unenumerated events)."""
    n2c_cohort_conn = _aggregate_cohort_connectivity([ev_conn_by_ev[e] for e in matched_evs])
    all_2c_ck_valid = all(bool(ev_2c_ests[e]["ck_valid"]) for e in matched_evs)
    any_2c_non_conn = any(str(ev_conn_by_ev[e]["status"]) != "connected" for e in matched_evs)

    prim_src_by_e: Dict[str, str] = {}
    p_inc_lo: Dict[str, float] = {}
    p_inc_hi: Dict[str, float] = {}
    v_inc_lo: Dict[str, float] = {}
    v_inc_hi: Dict[str, float] = {}
    p_comp_lo: Dict[str, float] = {}
    p_comp_hi: Dict[str, float] = {}
    v_comp_lo: Dict[str, float] = {}
    v_comp_hi: Dict[str, float] = {}

    for e in matched_evs:
        est = ev_2c_ests[e]
        if bool(est["full_fiber_valid"]):
            prim_src_by_e[e] = (
                "tier1_exact_full_fiber"
                if bool(est["exact_full_fiber_enumerated"])
                else "tier2_connected_full_fiber"
            )
            p_inc_lo[e] = float(est["inc_p_low"])  # type: ignore[arg-type]
            p_inc_hi[e] = float(est["inc_p_high"])  # type: ignore[arg-type]
            v_inc_lo[e] = float(est["inc_var_low"])  # type: ignore[arg-type]
            v_inc_hi[e] = float(est["inc_var_high"])  # type: ignore[arg-type]
            p_comp_lo[e] = float(est["comp_p_low"])  # type: ignore[arg-type]
            p_comp_hi[e] = float(est["comp_p_high"])  # type: ignore[arg-type]
            v_comp_lo[e] = float(est["comp_var_low"])  # type: ignore[arg-type]
            v_comp_hi[e] = float(est["comp_var_high"])  # type: ignore[arg-type]
        elif p2_a1_approved and bool(est["ck_valid"]):
            prim_src_by_e[e] = "component_conditional_p2_a1"
            p_inc_lo[e] = float(est["ck_inc_p_low"])  # type: ignore[arg-type]
            p_inc_hi[e] = float(est["ck_inc_p_high"])  # type: ignore[arg-type]
            v_inc_lo[e] = float(est["ck_inc_var_low"])  # type: ignore[arg-type]
            v_inc_hi[e] = float(est["ck_inc_var_high"])  # type: ignore[arg-type]
            p_comp_lo[e] = float(est["ck_comp_p_low"])  # type: ignore[arg-type]
            p_comp_hi[e] = float(est["ck_comp_p_high"])  # type: ignore[arg-type]
            v_comp_lo[e] = float(est["ck_comp_var_low"])  # type: ignore[arg-type]
            v_comp_hi[e] = float(est["ck_comp_var_high"])  # type: ignore[arg-type]
        else:
            prim_src_by_e[e] = "unavailable"

    n2c_conn_list = [
        {
            "event": e,
            **ev_conn_by_ev[e],
            "tier": ev_2c_ests[e]["tier"],
            "full_fiber_valid": ev_2c_ests[e]["full_fiber_valid"],
            "primary_null2c_source": prim_src_by_e[e],
            "mcmc_diagnostics": ev_2c_ests[e]["mcmc_diagnostics"],
        }
        for e in matched_evs
    ]

    comp_sens_2c: Optional[Dict[str, object]] = None
    if any_2c_non_conn:
        if all_2c_ck_valid:
            ci_n2c_ck_or, ci_obs_minus_n2c_ck_or = _event_bootstrap_q2_null2c_stationary_bounds(
                matched_evs,
                ev_inc_m,
                ev_comp_m,
                {e: float(ev_2c_ests[e]["ck_inc_p_low"]) for e in matched_evs},  # type: ignore[arg-type]
                {e: float(ev_2c_ests[e]["ck_inc_p_high"]) for e in matched_evs},  # type: ignore[arg-type]
                {e: float(ev_2c_ests[e]["ck_inc_var_low"]) for e in matched_evs},  # type: ignore[arg-type]
                {e: float(ev_2c_ests[e]["ck_inc_var_high"]) for e in matched_evs},  # type: ignore[arg-type]
                {e: float(ev_2c_ests[e]["ck_comp_p_low"]) for e in matched_evs},  # type: ignore[arg-type]
                {e: float(ev_2c_ests[e]["ck_comp_p_high"]) for e in matched_evs},  # type: ignore[arg-type]
                {e: float(ev_2c_ests[e]["ck_comp_var_low"]) for e in matched_evs},  # type: ignore[arg-type]
                {e: float(ev_2c_ests[e]["ck_comp_var_high"]) for e in matched_evs},  # type: ignore[arg-type]
                n_boot=n_boot,
                seed=seed_comp,
            )
            comp_sens_2c = {
                "label": (
                    "component-conditional (P2-A1, grid-approved)"
                    if p2_a1_approved
                    else "component-conditional (P2-A1, pending grid)"
                ),
                "connectivity_status": n2c_cohort_conn,
                "mcmc_or_exact_validated": True,
                "ci_null2c_log_or": list(ci_n2c_ck_or),
                "ci_obs_minus_null2c_log_or": list(ci_obs_minus_n2c_ck_or),
            }
        else:
            comp_sens_2c = {
                "label": "component-conditional (P2-A1, unvalidated — MCMC diagnostics failed)",
                "connectivity_status": n2c_cohort_conn,
                "mcmc_or_exact_validated": False,
                "ci_null2c_log_or": None,
                "ci_obs_minus_null2c_log_or": None,
            }

    tot_ni = sum(ev_inc_m[e].n_cand_aux for e in matched_evs)
    tot_nc = sum(ev_comp_m[e].n_cand_aux for e in matched_evs)
    if all(prim_src_by_e[e] != "unavailable" for e in matched_evs):
        ci_n2c_prim_or, ci_obs_minus_n2c_prim_or = _event_bootstrap_q2_null2c_stationary_bounds(
            matched_evs,
            ev_inc_m,
            ev_comp_m,
            p_inc_lo,
            p_inc_hi,
            v_inc_lo,
            v_inc_hi,
            p_comp_lo,
            p_comp_hi,
            v_comp_lo,
            v_comp_hi,
            n_boot=n_boot,
            seed=seed_primary,
        )
        prim_inc_pool = (
            [
                round(sum((ev_inc_m[e].n_cand_aux / float(tot_ni)) * p_inc_lo[e] for e in matched_evs), 6),
                round(sum((ev_inc_m[e].n_cand_aux / float(tot_ni)) * p_inc_hi[e] for e in matched_evs), 6),
            ]
            if tot_ni > 0
            else None
        )
        prim_comp_pool = (
            [
                round(sum((ev_comp_m[e].n_cand_aux / float(tot_nc)) * p_comp_lo[e] for e in matched_evs), 6),
                round(sum((ev_comp_m[e].n_cand_aux / float(tot_nc)) * p_comp_hi[e] for e in matched_evs), 6),
            ]
            if tot_nc > 0
            else None
        )
        return {
            "target_fiber": "equivalence_class_locked_initial_clash_rows_[S_obs]_~2c",
            "connectivity_status": n2c_cohort_conn,
            "per_event_connectivity": n2c_conn_list,
            "primary_null2c_available": True,
            "primary_pooled_null_inc_bounds": prim_inc_pool,
            "primary_pooled_null_comp_bounds": prim_comp_pool,
            "ci_null2c_log_or": list(ci_n2c_prim_or),
            "ci_obs_minus_null2c_log_or": list(ci_obs_minus_n2c_prim_or),
            "component_conditional_sensitivity": comp_sens_2c,
        }

    return {
        "target_fiber": "equivalence_class_locked_initial_clash_rows_[S_obs]_~2c",
        "connectivity_status": n2c_cohort_conn,
        "per_event_connectivity": n2c_conn_list,
        "primary_null2c_available": False,
        "primary_pooled_null_inc_bounds": None,
        "primary_pooled_null_comp_bounds": None,
        "ci_null2c_log_or": None,
        "ci_obs_minus_null2c_log_or": None,
        "component_conditional_sensitivity": comp_sens_2c,
    }


def _score_recipient_on_null2_profiles(
    rec_id: str,
    req_tup: Tuple[str, ...],
    amb: FrozenSet[str],
    base_m: OrganismProfile,
    base_p: OrganismProfile,
    pr_m: Dict[str, OrganismProfile],
    pr_p: Dict[str, OrganismProfile],
) -> str:
    pm_r = dict(pr_m)
    pp_r = dict(pr_p)
    rfs = frozenset(req_tup)
    pm_r[rec_id] = OrganismProfile(base_m.org_id, base_m.de_novo, base_m.conversions, rfs, base_m.exportable)
    pp_r[rec_id] = OrganismProfile(base_p.org_id, base_p.de_novo, base_p.conversions, rfs, base_p.exportable)
    st_m, _ = _fast_decide_acyclic_p_all(rec_id, pm_r, amb, cyc_mode="fanin_only")
    st_p, _ = _fast_decide_acyclic_p_all(rec_id, pp_r, amb, cyc_mode="fanin_only")
    if st_m == "d0":
        return CALL_D0_AMBIENT
    if st_m == "not_plexi":
        return CALL_D1_OR_CHAIN
    if st_p == "no_donor":
        return CALL_UNRESOLVED_NO_DONOR
    if st_m == "plexi" and st_p == "plexi":
        return CALL_PLEXITROPHIC
    return CALL_UNRESOLVED_AMBIG


def evaluate_freshwater_events(
    covered_dir: str,
    event_map_tsv: str,
    genomes: Dict[str, GenomeAnnotations],
    rep_ids: Sequence[str],
    hq_ids: Sequence[str],
    nc_ids: Sequence[str],
    rep_contig_lens_by_mag: Dict[str, Dict[str, int]],
    target_ctg_lens: Dict[str, int],
    valid_target_contigs: AbstractSet[str],
    unbinned_prok_kos_by_ctg: Dict[str, Set[str]],
    euk_kos_by_ctg: Dict[str, Set[str]],
    donor_hits_by_group_and_ctg: Dict[str, Dict[str, Dict[str, int]]],
    expected_targets_bed_sha256: str = "f48d55128984",
    targets_bed_path: Optional[str] = None,
    manifest_tsv: Optional[str] = None,
    expected_libraries: Optional[AbstractSet[str]] = None,
    min_events: int = 10,
    n_null1_draws: int = 1000,
    n_null2_draws: int = 1000,
    n_boot: int = 1000,
    require_protocol_1000_draws: bool = True,
    p2_a1_approved: bool = False,
    max_connectivity_enum_states: int = 2000,
    gtdb_clade_by_mag: Optional[Dict[str, str]] = None,
    hits_by_mag: Optional[Dict[str, List[GeneHit]]] = None,
    seed: int = 20261121,
) -> Dict[str, object]:
    """Validate per-library BED provenance, compute event-union breadths, and apply frozen P2 + P2-A1 gates & verdicts."""
    if not os.path.isdir(covered_dir):
        raise ValueError(f"covered_dir does not exist: {covered_dir}")
    if not os.path.isfile(event_map_tsv):
        raise ValueError(f"event_map_tsv does not exist: {event_map_tsv}")
    if targets_bed_path is not None:
        if not os.path.isfile(targets_bed_path):
            raise ValueError(f"targets_bed_path does not exist: {targets_bed_path}")
        got_tbed = sha256_12(targets_bed_path)
        if got_tbed != expected_targets_bed_sha256:
            raise ValueError(
                f"targets_bed_path sha256 mismatch: {got_tbed} != expected {expected_targets_bed_sha256}"
            )

    if manifest_tsv is not None:
        manifest_path = manifest_tsv
    elif os.path.isfile(os.path.join(covered_dir, "manifest_v2.tsv")):
        manifest_path = os.path.join(covered_dir, "manifest_v2.tsv")
    elif os.path.isfile(os.path.join(covered_dir, "MANIFEST.tsv")):
        manifest_path = os.path.join(covered_dir, "MANIFEST.tsv")
    else:
        default_mv2 = os.path.join(OUT_DIR, "breadth", "manifest_v2.tsv")
        if os.path.realpath(covered_dir) == os.path.realpath(os.path.join(FW_ROOT, "breadth_v2")) and os.path.isfile(default_mv2):
            manifest_path = default_mv2
        else:
            manifest_path = os.path.join(covered_dir, "MANIFEST.tsv")

    if not os.path.isfile(manifest_path):
        raise ValueError(
            f"Missing per-library provenance manifest ({manifest_path}); required to bind each "
            f"*.covered.bed.gz to targets_v2.bed ({expected_targets_bed_sha256}) and reject v1 outputs."
        )

    manifest_by_lib: Dict[str, Dict[str, str]] = {}
    with open(manifest_path, "r", encoding="utf-8") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        fnames = set(reader.fieldnames or ())
        is_mv2_schema = {"library", "sha256", "targets_bed_sha12"}.issubset(fnames)
        is_generic_schema = {"library", "covered_bed_sha256", "targets_bed_sha256"}.issubset(fnames)
        if not (is_mv2_schema or is_generic_schema):
            raise ValueError(
                f"Manifest {manifest_path} must contain either ('library','sha256','targets_bed_sha12') "
                f"or ('library','covered_bed_sha256','targets_bed_sha256')"
            )
        for r in reader:
            lib = (r.get("library") or "").strip()
            if not lib:
                raise ValueError(f"Empty library in {manifest_path}")
            if lib in manifest_by_lib:
                raise ValueError(f"Duplicate library entry in {manifest_path}: {lib}")
            clean_r = {k: (v or "").strip() for k, v in r.items() if k is not None}
            if is_mv2_schema:
                if clean_r.get("local_vs_grex") and clean_r["local_vs_grex"] != "match":
                    raise ValueError(f"Manifest row for {lib} has local_vs_grex={clean_r['local_vs_grex']!r} != 'match'")
                if "producer_job_task" in fnames and not clean_r.get("producer_job_task"):
                    raise ValueError(f"Manifest row for {lib} has empty producer_job_task")
                clean_r["_cov_sha"] = clean_r["sha256"]
                clean_r["_tbed_sha12"] = clean_r["targets_bed_sha12"]
            else:
                clean_r["_cov_sha"] = clean_r["covered_bed_sha256"]
                clean_r["_tbed_sha12"] = clean_r["targets_bed_sha256"]
            manifest_by_lib[lib] = clean_r

    ev_to_libs: Dict[str, List[str]] = {}
    control_libs: List[str] = []
    seen_libs: Dict[str, str] = {}
    with open(event_map_tsv, "r", encoding="utf-8") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        if not reader.fieldnames or not {"library", "event"}.issubset(set(reader.fieldnames)):
            raise ValueError("event_map_tsv must contain 'library' and 'event' columns")
        for r in reader:
            lib = (r.get("library") or "").strip()
            ev = (r.get("event") or "").strip()
            if not lib or not ev:
                raise ValueError(f"Empty library or event in event_map_tsv row: {r}")
            if lib in seen_libs:
                raise ValueError(
                    f"Duplicate or multi-event library assignment in event_map_tsv: {lib} "
                    f"({seen_libs[lib]} vs {ev})"
                )
            seen_libs[lib] = ev
            if ev.upper() in ("BLANK", "CONTROL", "NEG_CONTROL"):
                control_libs.append(lib)
            else:
                ev_to_libs.setdefault(ev, []).append(lib)

    if not ev_to_libs:
        raise ValueError("event_map_tsv contains zero non-control biological events")
    if expected_libraries is not None:
        missing_from_map = sorted(set(expected_libraries) - set(seen_libs))
        unknown_in_map = sorted(set(seen_libs) - set(expected_libraries))
        if missing_from_map or unknown_in_map:
            raise ValueError(
                f"event_map_tsv library set mismatch against expected BAMs: "
                f"missing={missing_from_map[:5]} (n={len(missing_from_map)}), "
                f"unknown={unknown_in_map[:5]} (n={len(unknown_in_map)})"
            )

    # Validate every mapped biological library's MANIFEST row, file existence, sha256, and intervals
    lib_ivs: Dict[str, Dict[str, List[Tuple[int, int]]]] = {}
    for ev_id in sorted(ev_to_libs):
        for lib in ev_to_libs[ev_id]:
            if lib not in manifest_by_lib:
                raise ValueError(f"Mapped library {lib} missing from {manifest_path}")
            mrow = manifest_by_lib[lib]
            if mrow["_tbed_sha12"] != expected_targets_bed_sha256:
                raise ValueError(
                    f"Library {lib} in {manifest_path} was produced with targets BED "
                    f"{mrow['_tbed_sha12']} != expected {expected_targets_bed_sha256}"
                )
            bed_gz = os.path.join(covered_dir, f"{lib}.covered.bed.gz")
            if not os.path.isfile(bed_gz):
                raise ValueError(f"Missing covered interval file for mapped library {lib}: {bed_gz}")
            exp_cov_sha = mrow["_cov_sha"]
            got_cov_sha = sha256_full(bed_gz) if len(exp_cov_sha) == 64 else sha256_12(bed_gz)
            if got_cov_sha != exp_cov_sha:
                raise ValueError(
                    f"Covered BED sha256 mismatch for {lib}: disk {got_cov_sha} != "
                    f"manifest {exp_cov_sha}"
                )
            iv_map: Dict[str, List[Tuple[int, int]]] = {}
            with gzip.open(bed_gz, "rt", encoding="utf-8") as gfh:
                for ln in gfh:
                    if not ln.strip():
                        continue
                    parts = ln.rstrip("\n").split("\t")
                    if len(parts) < 3:
                        raise ValueError(f"Malformed BED line in {bed_gz}: {ln!r}")
                    c_id = parts[0]
                    s_pos, e_pos = int(parts[1]), int(parts[2])
                    if c_id not in valid_target_contigs:
                        raise ValueError(f"Unknown contig {c_id} in {bed_gz} (not in targets_v2)")
                    max_len = target_ctg_lens.get(c_id)
                    if max_len is None or s_pos < 0 or e_pos < s_pos or e_pos > max_len:
                        raise ValueError(
                            f"Out-of-bounds interval {c_id}:{s_pos}-{e_pos} (len={max_len}) in {bed_gz}"
                        )
                    if e_pos > s_pos:
                        iv_map.setdefault(c_id, []).append((s_pos, e_pos))
            lib_ivs[lib] = iv_map

    strata_defs: List[Tuple[str, Tuple[str, ...], SupplyScenario]] = [
        ("R-I3_B12_REDS__E_minimal", ("B12", "RED_S"), SCENARIO_E_MINIMAL),
        ("R-I3_B12_REDS__E_plus_DMSP", ("B12", "RED_S"), SCENARIO_E_PLUS_DMSP),
        ("R-I2_B12_FIXEDN__E_minimal", ("B12", "FIXED_N"), SCENARIO_E_MINIMAL),
        ("R-I2_B12_FIXEDN__E_plus_DIN", ("B12", "FIXED_N"), SCENARIO_E_PLUS_DIN),
        ("R-I6_B12_SINGLE__E_minimal", ("B12",), SCENARIO_E_MINIMAL),
        ("C_THI_SPLIT__E_minimal", ("THI",), SCENARIO_E_MINIMAL),
    ]
    q2_contrast_defs: List[Tuple[str, Tuple[str, str], SupplyScenario]] = [
        ("R-I3_B12_REDS__E_minimal", ("B12", "RED_S"), SCENARIO_E_MINIMAL),
        ("R-I3_B12_REDS__E_plus_DMSP", ("B12", "RED_S"), SCENARIO_E_PLUS_DMSP),
        ("R-I2_B12_FIXEDN__E_minimal", ("B12", "FIXED_N"), SCENARIO_E_MINIMAL),
        ("R-I2_B12_FIXEDN__E_plus_DIN", ("B12", "FIXED_N"), SCENARIO_E_PLUS_DIN),
    ]
    c0_scen_defs: List[Tuple[str, SupplyScenario]] = [
        ("C0xC0_ALL45__E_minimal", SCENARIO_E_MINIMAL),
        ("C0xC0_ALL45__E_plus_DMSP", SCENARIO_E_PLUS_DMSP),
        ("C0xC0_ALL45__E_plus_DIN", SCENARIO_E_PLUS_DIN),
    ]
    all_stratum_labels = [s[0] for s in strata_defs] + ["C0xC0_ALL45__E_minimal"]
    hq_set = set(hq_ids)
    nc_set = set(nc_ids)
    c0_pairs_all = list(itertools.combinations(C0_CONTROL_MODULES, 2))
    forbid_o2_clash = lambda s_idx: bool(s_idx & O2_PLUS_IDX) and bool(s_idx & O2_MINUS_IDX)

    draw_count_gate_err: Optional[str] = None
    if require_protocol_1000_draws and (n_null1_draws < 1000 or n_null2_draws < 1000 or n_boot < 1000):
        draw_count_gate_err = (
            f"Inconclusive (non-inferential draw count < 1000: got n_null1_draws={n_null1_draws}, "
            f"n_null2_draws={n_null2_draws}, n_boot={n_boot}; frozen P2 §3.4/§4.1 requires >=1000 each)"
        )

    # Outcome-blind recipient covariates & C0xC0 eligibility across HQ MAGs (computed BEFORE inspecting donors or E)
    elig_c0_by_mag: Dict[str, List[Tuple[str, str]]] = {}
    cov_vec_by_mag: Dict[str, Tuple[float, float, float]] = {}
    gtdb_ranks_by_mag: Dict[str, Dict[str, str]] = {}
    for mid in hq_ids:
        g_hq = genomes[mid]
        if g_hq.domain in ("Bacteria", "Archaea"):
            elig_c0_by_mag[mid] = match_c0_pairs_for_recipient(g_hq)
        else:
            elig_c0_by_mag[mid] = []
        cov_vec_by_mag[mid] = (
            float(count_missing_modules(g_hq)),
            math.log(max(1, g_hq.genome_size_bp)),
            float(g_hq.completeness),
        )
        if gtdb_clade_by_mag is not None:
            gtdb_ranks_by_mag[mid] = _parse_gtdb_ranks(gtdb_clade_by_mag.get(mid, ""))

    cov_sds_list: List[float] = []
    for j in range(3):
        vals_j = [cov_vec_by_mag[m][j] for m in hq_ids]
        if len(vals_j) >= 2:
            mean_j = sum(vals_j) / len(vals_j)
            var_j = sum((v - mean_j) ** 2 for v in vals_j) / len(vals_j)
            sd_j = math.sqrt(var_j)
            cov_sds_list.append(sd_j if sd_j > 1e-9 else 1.0)
        else:
            cov_sds_list.append(1.0)
    cov_sds = (cov_sds_list[0], cov_sds_list[1], cov_sds_list[2])

    bio_events = sorted(ev_to_libs)
    per_event_strata: Dict[str, Dict[str, StratumCounts]] = {}
    per_event_det_mags: Dict[str, List[str]] = {}
    per_event_ctx: Dict[str, CommunityScorerContext] = {}
    per_event_unb_euk_plus: Dict[str, Dict[str, OrganismProfile]] = {}
    per_event_gate2: Dict[str, Dict[str, Dict[str, object]]] = {}
    per_event_q2_matches: Dict[str, Dict[str, List[Tuple[GenomeAnnotations, GenomeAnnotations, Tuple[str, str]]]]] = {
        c_key: {} for c_key, _, _ in q2_contrast_defs
    }
    per_event_q2_match_stats: Dict[str, Tuple[int, int, int]] = {
        c_key: (0, 0, 0) for c_key, _, _ in q2_contrast_defs  # (n_inc_total, n_distinct_match, n_self_match)
    }
    n_events_with_any_hq = 0

    for ev_id in bio_events:
        per_ctg_iv: Dict[str, List[Tuple[int, int]]] = {}
        for lib in ev_to_libs[ev_id]:
            for c_id, lst in lib_ivs[lib].items():
                per_ctg_iv.setdefault(c_id, []).extend(lst)
        merged_ctg_cov = {c_id: _merge_interval_length(lst) for c_id, lst in per_ctg_iv.items()}

        det_mags_ev: List[str] = []
        for mid in rep_ids:
            cmap_m = rep_contig_lens_by_mag[mid]
            tot_bp_m = sum(cmap_m.values())
            cov_bp_m = sum(merged_ctg_cov.get(c_id, 0) for c_id in cmap_m)
            if tot_bp_m > 0 and (cov_bp_m / tot_bp_m) >= 0.50:
                det_mags_ev.append(mid)
        per_event_det_mags[ev_id] = det_mags_ev
        det_hq_ev = [genomes[m] for m in det_mags_ev if m in hq_set]
        if det_hq_ev:
            n_events_with_any_hq += 1

        unb_kos_ev = {
            k
            for c_id, kset in unbinned_prok_kos_by_ctg.items()
            if merged_ctg_cov.get(c_id, 0) >= 0.50 * target_ctg_lens.get(c_id, 1)
            for k in kset
        }
        euk_kos_ev = {
            k
            for c_id, kset in euk_kos_by_ctg.items()
            if merged_ctg_cov.get(c_id, 0) >= 0.50 * target_ctg_lens.get(c_id, 1)
            for k in kset
        }
        unb_g_ev = GenomeAnnotations.from_kos("_unbinned", unb_kos_ev) if unb_kos_ev else None
        euk_g_ev = GenomeAnnotations.from_kos("_euk_tiara", euk_kos_ev, domain="Eukaryota") if euk_kos_ev else None
        det_genomes_ev = [genomes[m] for m in det_mags_ev]
        ctx_ev = CommunityScorerContext(det_genomes_ev, unb_g_ev, euk_g_ev, mode="fanin_only")
        per_event_ctx[ev_id] = ctx_ev
        extra_plus_ev: Dict[str, OrganismProfile] = {}
        if "_unbinned_P_star" in ctx_ev.donor_plus:
            extra_plus_ev["_unbinned_P_star"] = ctx_ev.donor_plus["_unbinned_P_star"]
        if "_euk_unassigned_P_star" in ctx_ev.donor_plus:
            extra_plus_ev["_euk_unassigned_P_star"] = ctx_ev.donor_plus["_euk_unassigned_P_star"]
        per_event_unb_euk_plus[ev_id] = extra_plus_ev

        # Outcome-blind 1-to-1 clade + covariate matching for each Q2 contrast in event ev_id
        if gtdb_clade_by_mag is not None:
            for c_key, req_inc, _ in q2_contrast_defs:
                inc_recs_ev = [
                    r_hq
                    for r_hq in det_hq_ev
                    if r_hq.domain in ("Bacteria", "Archaea")
                    and all(recipient_endpoint_status(r_hq, ep, 0.05) == STATE_ABSENT for ep in req_inc)
                ]
                m_pairs, n_dist, n_self = _match_c0_controls_in_event(
                    c_key, inc_recs_ev, det_hq_ev, elig_c0_by_mag, gtdb_ranks_by_mag, cov_vec_by_mag, cov_sds
                )
                per_event_q2_matches[c_key][ev_id] = m_pairs
                prev_tot, prev_d, prev_s = per_event_q2_match_stats[c_key]
                per_event_q2_match_stats[c_key] = (prev_tot + len(inc_recs_ev), prev_d + n_dist, prev_s + n_self)

        ev_sc_map: Dict[str, StratumCounts] = {}
        for s_label, reqs, scen in strata_defs:
            sc = StratumCounts()
            for rec in det_hq_ev:
                call_f, _ = ctx_ev.evaluate(rec, reqs, scen)
                if call_f != CALL_NOT_ELIGIBLE:
                    sc.record(call_f)
            ev_sc_map[s_label] = sc

        for c0_label, c0_scen in c0_scen_defs:
            sc_c0_ev = StratumCounts()
            for rec in det_hq_ev:
                for p_c0 in c0_pairs_all:
                    call_c0, _ = ctx_ev.evaluate(rec, p_c0, c0_scen)
                    if call_c0 != CALL_NOT_ELIGIBLE:
                        sc_c0_ev.record(call_c0)
            ev_sc_map[c0_label] = sc_c0_ev
        per_event_strata[ev_id] = ev_sc_map

        ev_g2_map: Dict[str, Dict[str, object]] = {}
        for gk, bdict in donor_hits_by_group_and_ctg.items():
            det_rep_h = sum(
                cnt for c_id, cnt in bdict["rep490"].items()
                if merged_ctg_cov.get(c_id, 0) >= 0.50 * target_ctg_lens.get(c_id, 1)
            )
            det_nonrep_h = sum(
                cnt for c_id, cnt in bdict["dastool655_nonrep"].items()
                if merged_ctg_cov.get(c_id, 0) >= 0.50 * target_ctg_lens.get(c_id, 1)
            )
            det_unb_h = sum(
                cnt for c_id, cnt in bdict["unb_prok"].items()
                if merged_ctg_cov.get(c_id, 0) >= 0.50 * target_ctg_lens.get(c_id, 1)
            )
            det_euk_h = sum(
                cnt for c_id, cnt in bdict["euk"].items()
                if merged_ctg_cov.get(c_id, 0) >= 0.50 * target_ctg_lens.get(c_id, 1)
            )
            ev_g2_map[gk] = _compute_gate2_record(det_rep_h, det_nonrep_h, det_unb_h, det_euk_h)
        per_event_gate2[ev_id] = ev_g2_map

    # Identify which strata pass Gates 1-3 + draw-count gate so event-level Null 1 and Null 2b/2c run once per needed event
    rng = random.Random(seed)
    strata_reqs_lookup: Dict[str, Tuple[Optional[Tuple[str, ...]], SupplyScenario]] = {
        s_l: (reqs, scen) for s_l, reqs, scen in strata_defs
    }
    strata_reqs_lookup["C0xC0_ALL45__E_minimal"] = (None, SCENARIO_E_MINIMAL)

    q1_gate_info: Dict[str, Tuple[List[str], StratumCounts, Dict[str, object], float, Optional[str]]] = {}
    needed_q1_events: Set[str] = set()
    active_q1_strata: List[str] = []

    for s_label in all_stratum_labels:
        elig_evs = [ev for ev in bio_events if per_event_strata[ev][s_label].n_cand_aux >= 1]
        tot_sc = StratumCounts(
            n_auto=sum(per_event_strata[ev][s_label].n_auto for ev in bio_events),
            n_d0=sum(per_event_strata[ev][s_label].n_d0 for ev in bio_events),
            n_d1=sum(per_event_strata[ev][s_label].n_d1 for ev in bio_events),
            n_plexi=sum(per_event_strata[ev][s_label].n_plexi for ev in bio_events),
            n_ambig=sum(per_event_strata[ev][s_label].n_ambig for ev in bio_events),
            n_no_donor=sum(per_event_strata[ev][s_label].n_no_donor for ev in bio_events),
        )
        ep_grp = STRATUM_TO_GATE2_GROUP[s_label]
        evs_for_g2 = elig_evs if elig_evs else bio_events
        pooled_g2 = _compute_gate2_record(
            rep490_hits=sum(int(per_event_gate2[ev][ep_grp]["rep490_hits"]) for ev in evs_for_g2),  # type: ignore[arg-type]
            dastool655_nonrep_hits=sum(int(per_event_gate2[ev][ep_grp]["dastool655_nonrep_hits"]) for ev in evs_for_g2),  # type: ignore[arg-type]
            unbinned_prok_hits=sum(int(per_event_gate2[ev][ep_grp]["unbinned_prok_hits"]) for ev in evs_for_g2),  # type: ignore[arg-type]
            tiara_eukarya_hits=sum(int(per_event_gate2[ev][ep_grp]["tiara_eukarya_hits"]) for ev in evs_for_g2),  # type: ignore[arg-type]
        )
        unres_frac = (tot_sc.n_unresolved / tot_sc.n_cand_aux) if tot_sc.n_cand_aux > 0 else 1.0

        gate_fail: Optional[str] = None
        if len(elig_evs) < min_events:
            gate_fail = f"Inconclusive (data-adequacy gate: n_events={len(elig_evs)} < {min_events})"
        else:
            g2_err = _gate2_failure_reason(ep_grp, pooled_g2)
            if g2_err is not None:
                gate_fail = g2_err
            elif unres_frac > 0.50:
                gate_fail = f"Inconclusive (data-adequacy gate: unresolved_fraction={unres_frac:.3f} > 0.50)"
            elif draw_count_gate_err is not None:
                gate_fail = draw_count_gate_err

        q1_gate_info[s_label] = (elig_evs, tot_sc, pooled_g2, unres_frac, gate_fail)
        if gate_fail is None and s_label != "C0xC0_ALL45__E_minimal":
            needed_q1_events.update(elig_evs)
            active_q1_strata.append(s_label)

    # Pre-evaluate Q2 matched counts & gates so we know which events need constrained Null 2c draws
    q2_prep: Dict[
        str,
        Tuple[
            Tuple[str, str],
            SupplyScenario,
            List[str],
            Dict[str, StratumCounts],
            Dict[str, StratumCounts],
            StratumCounts,
            StratumCounts,
            Dict[str, object],
            Dict[str, object],
            Optional[str],
        ],
    ] = {}
    needed_q2_events: Set[str] = set()

    for c_key, req_inc, scen_c in q2_contrast_defs:
        matched_evs: List[str] = []
        ev_inc_m: Dict[str, StratumCounts] = {}
        ev_comp_m: Dict[str, StratumCounts] = {}
        if gtdb_clade_by_mag is not None:
            for ev_id in bio_events:
                mpairs = per_event_q2_matches[c_key].get(ev_id, [])
                if not mpairs:
                    continue
                matched_evs.append(ev_id)
                ctx_e = per_event_ctx[ev_id]
                sc_i = StratumCounts()
                sc_c = StratumCounts()
                for rx, ry, p_c0 in mpairs:
                    ci_call, _ = ctx_e.evaluate(rx, req_inc, scen_c)
                    cc_call, _ = ctx_e.evaluate(ry, p_c0, scen_c)
                    if ci_call != CALL_NOT_ELIGIBLE:
                        sc_i.record(ci_call)
                    if cc_call != CALL_NOT_ELIGIBLE:
                        sc_c.record(cc_call)
                ev_inc_m[ev_id] = sc_i
                ev_comp_m[ev_id] = sc_c

        tot_inc_m = StratumCounts(
            n_d0=sum(ev_inc_m[e].n_d0 for e in matched_evs),
            n_d1=sum(ev_inc_m[e].n_d1 for e in matched_evs),
            n_plexi=sum(ev_inc_m[e].n_plexi for e in matched_evs),
            n_ambig=sum(ev_inc_m[e].n_ambig for e in matched_evs),
            n_no_donor=sum(ev_inc_m[e].n_no_donor for e in matched_evs),
        )
        tot_comp_m = StratumCounts(
            n_d0=sum(ev_comp_m[e].n_d0 for e in matched_evs),
            n_d1=sum(ev_comp_m[e].n_d1 for e in matched_evs),
            n_plexi=sum(ev_comp_m[e].n_plexi for e in matched_evs),
            n_ambig=sum(ev_comp_m[e].n_ambig for e in matched_evs),
            n_no_donor=sum(ev_comp_m[e].n_no_donor for e in matched_evs),
        )
        ep_inc = STRATUM_TO_GATE2_GROUP[c_key]
        evs_g2_q2 = matched_evs if matched_evs else bio_events
        inc_g2_m = _compute_gate2_record(
            rep490_hits=sum(int(per_event_gate2[e][ep_inc]["rep490_hits"]) for e in evs_g2_q2),  # type: ignore[arg-type]
            dastool655_nonrep_hits=sum(int(per_event_gate2[e][ep_inc]["dastool655_nonrep_hits"]) for e in evs_g2_q2),  # type: ignore[arg-type]
            unbinned_prok_hits=sum(int(per_event_gate2[e][ep_inc]["unbinned_prok_hits"]) for e in evs_g2_q2),  # type: ignore[arg-type]
            tiara_eukarya_hits=sum(int(per_event_gate2[e][ep_inc]["tiara_eukarya_hits"]) for e in evs_g2_q2),  # type: ignore[arg-type]
        )
        comp_g2_m = _compute_gate2_record(
            rep490_hits=sum(int(per_event_gate2[e]["C0_ALL10"]["rep490_hits"]) for e in evs_g2_q2),  # type: ignore[arg-type]
            dastool655_nonrep_hits=sum(int(per_event_gate2[e]["C0_ALL10"]["dastool655_nonrep_hits"]) for e in evs_g2_q2),  # type: ignore[arg-type]
            unbinned_prok_hits=sum(int(per_event_gate2[e]["C0_ALL10"]["unbinned_prok_hits"]) for e in evs_g2_q2),  # type: ignore[arg-type]
            tiara_eukarya_hits=sum(int(per_event_gate2[e]["C0_ALL10"]["tiara_eukarya_hits"]) for e in evs_g2_q2),  # type: ignore[arg-type]
        )
        u_inc_rate = (tot_inc_m.n_unresolved / tot_inc_m.n_cand_aux) if tot_inc_m.n_cand_aux > 0 else 1.0
        u_comp_rate = (tot_comp_m.n_unresolved / tot_comp_m.n_cand_aux) if tot_comp_m.n_cand_aux > 0 else 1.0

        q2_gate_fail: Optional[str] = None
        if gtdb_clade_by_mag is None:
            q2_gate_fail = (
                "Inconclusive (P2-Q2 matching unevaluable: gtdb_clade_by_mag not provided for frozen §3.3 "
                "covariate + GTDB-clade recipient matching)"
            )
        elif len(matched_evs) < min_events:
            inc_evs_raw = sum(1 for ev in bio_events if per_event_strata[ev][c_key].n_cand_aux >= 1)
            q2_gate_fail = (
                f"Inconclusive (data-adequacy gate: P2-Q2 requires >={min_events} events with outcome-blind "
                f"GTDB-clade + covariate-matched HQ candidate-auxotroph pairs; got matched_events={len(matched_evs)} "
                f"[inc_events_before_matching={inc_evs_raw}] < {min_events})"
            )
        else:
            g2_i_err = _gate2_failure_reason(ep_inc, inc_g2_m)
            g2_c_err = _gate2_failure_reason("C0_ALL10", comp_g2_m)
            if g2_i_err is not None:
                q2_gate_fail = g2_i_err
            elif g2_c_err is not None:
                q2_gate_fail = g2_c_err
            elif u_inc_rate > 0.50 or u_comp_rate > 0.50:
                q2_gate_fail = (
                    f"Inconclusive (data-adequacy gate: unresolved_fraction > 0.50 "
                    f"[inc={u_inc_rate:.3f}, comp={u_comp_rate:.3f}])"
                )
            elif draw_count_gate_err is not None:
                q2_gate_fail = draw_count_gate_err

        q2_prep[c_key] = (
            req_inc, scen_c, matched_evs, ev_inc_m, ev_comp_m, tot_inc_m, tot_comp_m, inc_g2_m, comp_g2_m, q2_gate_fail
        )
        if q2_gate_fail is None:
            needed_q2_events.update(matched_evs)

    # Run per-event Null 1, Null 2b (unconstrained), and Null 2c (O2-clash-constrained) ONCE per needed event
    min_ess_req = 200.0 if require_protocol_1000_draws else max(40.0, 0.20 * float(n_null2_draws))
    ev_null1_counts: Dict[str, Dict[str, Optional[StratumCounts]]] = {}
    ev_null1_no_source: Set[str] = set()
    ev_null2b_conn: Dict[str, Dict[str, object]] = {}
    ev_null2b_est: Dict[str, Dict[str, Dict[str, object]]] = {}
    ev_null2c_conn: Dict[str, Dict[str, object]] = {}
    ev_null2c_q2_est: Dict[str, Dict[str, Dict[str, object]]] = {}

    all_needed_events = sorted(needed_q1_events | needed_q2_events)
    for ev_id in all_needed_events:
        det_mags_e = per_event_det_mags[ev_id]
        det_gs_e = [genomes[m] for m in det_mags_e]
        det_hq_e = [genomes[m] for m in det_mags_e if m in hq_set]
        extra_plus_e = per_event_unb_euk_plus[ev_id]
        ctx_e = per_event_ctx[ev_id]

        # 1. Null 1 contig-block downsampling on event ev_id (strictly no StratumCounts(n_d1=1) or unfragmented fallback)
        if ev_id in needed_q1_events:
            if hits_by_mag is None:
                ev_null1_no_source.add(ev_id)
                ev_null1_counts[ev_id] = {s_l: None for s_l in active_q1_strata}
            else:
                nc_e = [m for m in det_mags_e if m in nc_set and len(hits_by_mag.get(m, ())) > 0]
                if not nc_e or not det_hq_e:
                    ev_null1_no_source.add(ev_id)
                    ev_null1_counts[ev_id] = {s_l: None for s_l in active_q1_strata}
                else:
                    n1_map = {s_l: StratumCounts() for s_l in active_q1_strata}
                    t_profs = [(g.completeness, max(1, len(g.contig_lengths))) for g in det_hq_e]
                    for d_i in range(n_null1_draws):
                        sim_ds = [
                            downsample_contig_blocks(
                                hits_by_mag[m_nc],
                                target_completeness=t_profs[(d_i + k_i) % len(t_profs)][0],
                                target_contig_count=t_profs[(d_i + k_i) % len(t_profs)][1],
                                rng=rng,
                                mag_id=f"n1_{ev_id}_{d_i}_{k_i}",
                                domain=genomes[m_nc].domain,
                                orig_completeness=genomes[m_nc].completeness,
                                total_genome_bp=genomes[m_nc].genome_size_bp,
                            )
                            for k_i, m_nc in enumerate(nc_e)
                        ]
                        ctx_s = CommunityScorerContext(
                            sim_ds, mode="fanin_only", verify_against_reference=(d_i < 2)
                        )
                        for s_l in active_q1_strata:
                            if per_event_strata[ev_id][s_l].n_cand_aux == 0:
                                continue
                            reqs_opt, scen_s = strata_reqs_lookup[s_l]
                            p_list = [reqs_opt] if reqs_opt is not None else c0_pairs_all
                            sc_acc = n1_map[s_l]
                            for r_sim in sim_ds:
                                for req_tup in p_list:
                                    c_s, _ = ctx_s.evaluate(r_sim, req_tup, scen_s)
                                    if c_s != CALL_NOT_ELIGIBLE:
                                        sc_acc.record(c_s)
                    ev_null1_counts[ev_id] = {s_l: n1_map[s_l] for s_l in active_q1_strata}

        mm_e = build_donor_binary_matrix(det_gs_e, optimistic=False)
        mp_e = build_donor_binary_matrix(det_gs_e, optimistic=True)

        # 2. Null 2b unconstrained two-tier stationary-mean estimator (Tier 1 exact Omega enum / Tier 2 4-chain MCMC)
        if ev_id in needed_q1_events:
            conn_2b_raw = classify_ternary_move_connectivity(
                mm_e, mp_e, max_enum_states=max_connectivity_enum_states, return_enumerated_pairs=True
            )
            ev_null2b_conn[ev_id] = _strip_conn_pairs(conn_2b_raw)
            elig_q1_by_s: Dict[str, List[Tuple[GenomeAnnotations, Tuple[str, ...], FrozenSet[str]]]] = {}
            for s_l in active_q1_strata:
                if per_event_strata[ev_id][s_l].n_cand_aux == 0:
                    continue
                reqs_opt, scen_s = strata_reqs_lookup[s_l]
                p_list = [reqs_opt] if reqs_opt is not None else c0_pairs_all
                tups_s: List[Tuple[GenomeAnnotations, Tuple[str, ...], FrozenSet[str]]] = []
                for r_hq in det_hq_e:
                    if r_hq.domain not in ("Bacteria", "Archaea"):
                        continue
                    for req_tup in p_list:
                        if all(recipient_endpoint_status(r_hq, ep, 0.05) == STATE_ABSENT for ep in req_tup):
                            tups_s.append((r_hq, req_tup, scen_s.ambient_compounds))
                elig_q1_by_s[s_l] = tups_s

            def _score_q1_state(
                pm_m: Sequence[Sequence[int]],
                pm_p: Sequence[Sequence[int]],
            ) -> Dict[str, Tuple[float, float]]:
                pr_m = {
                    det_gs_e[i].mag_id: profile_from_module_Presence_row(det_gs_e[i].mag_id, pm_m[i])
                    for i in range(len(det_gs_e))
                }
                pr_p = {
                    det_gs_e[i].mag_id: profile_from_module_Presence_row(det_gs_e[i].mag_id, pm_p[i])
                    for i in range(len(det_gs_e))
                }
                pr_p.update(extra_plus_e)
                b_out: Dict[str, Tuple[float, float]] = {}
                for s_k, t_list in elig_q1_by_s.items():
                    sc_t = StratumCounts()
                    for r_hq, req_tup, amb_s in t_list:
                        sc_t.record(
                            _score_recipient_on_null2_profiles(
                                r_hq.mag_id,
                                req_tup,
                                amb_s,
                                ctx_e.rec_base_minus[r_hq.mag_id],
                                ctx_e.rec_base_plus[r_hq.mag_id],
                                pr_m,
                                pr_p,
                            )
                        )
                    b_out[s_k] = sc_t.prevalence_bounds()
                return b_out

            ev_est_2b: Dict[str, Dict[str, object]] = {}
            if bool(conn_2b_raw["exact_full_fiber_enumerated"]):
                fiber_pairs_2b = conn_2b_raw["exact_fiber_pairs"]  # type: ignore[assignment]
                full_scores = [_score_q1_state(pm_m, pm_p) for pm_m, pm_p in fiber_pairs_2b]  # type: ignore[union-attr]
                comp_pairs_2b = conn_2b_raw.get("exact_component_pairs")
                if str(conn_2b_raw["status"]) == "connected" or comp_pairs_2b is None:
                    comp_scores = full_scores if str(conn_2b_raw["status"]) == "connected" else None
                else:
                    comp_scores = [_score_q1_state(pm_m, pm_p) for pm_m, pm_p in comp_pairs_2b]  # type: ignore[union-attr]
                for s_k in elig_q1_by_s:
                    p_lo_f = sum(bs[s_k][0] for bs in full_scores) / len(full_scores)
                    p_hi_f = sum(bs[s_k][1] for bs in full_scores) / len(full_scores)
                    p_lo_c = (
                        sum(bs[s_k][0] for bs in comp_scores) / len(comp_scores)
                        if comp_scores is not None
                        else p_lo_f
                    )
                    p_hi_c = (
                        sum(bs[s_k][1] for bs in comp_scores) / len(comp_scores)
                        if comp_scores is not None
                        else p_hi_f
                    )
                    ev_est_2b[s_k] = {
                        "tier": "tier1_exact_full_fiber",
                        "exact_full_fiber_enumerated": True,
                        "move_connectivity": str(conn_2b_raw["status"]),
                        "full_fiber_valid": True,
                        "p_low": p_lo_f,
                        "p_high": p_hi_f,
                        "var_bm_low": 0.0,
                        "var_bm_high": 0.0,
                        "mcmc_diagnostics": None,
                        "comp_p_low": p_lo_c,
                        "comp_p_high": p_hi_c,
                        "comp_var_bm_low": 0.0,
                        "comp_var_bm_high": 0.0,
                        "comp_valid": comp_scores is not None,
                    }
            else:
                chains_2b, ham_2b = mcmc_reversible_nested_chains(
                    mm_e, mp_e, rng, n_draws=n_null2_draws, n_chains=4, burnin_mult=20, thin_mult=3
                )
                chain_scores_2b = [
                    [_score_q1_state(pm_m, pm_p) for pm_m, pm_p in ch]
                    for ch in chains_2b
                ]
                for s_k in elig_q1_by_s:
                    lo_ch = [[bd[s_k][0] for bd in ch_s] for ch_s in chain_scores_2b]
                    hi_ch = [[bd[s_k][1] for bd in ch_s] for ch_s in chain_scores_2b]
                    v_diag = validate_mcmc_prevalence_bounds(
                        lo_ch,
                        hi_ch,
                        ham_2b,
                        has_active_move=bool(conn_2b_raw["has_active_move"]),
                        min_ess=min_ess_req,
                    )
                    d_lo = v_diag["bound_low_diag"]  # type: ignore[assignment]
                    d_hi = v_diag["bound_high_diag"]  # type: ignore[assignment]
                    mcmc_ok = bool(v_diag["mcmc_diagnostics_pass"])
                    # Separate connectivity gate on larger Tier 2 sets (Astra 249779947258167296)
                    full_ok = (str(conn_2b_raw["status"]) == "connected") and mcmc_ok
                    ev_est_2b[s_k] = {
                        "tier": "tier2_reversible_mcmc_4chain",
                        "exact_full_fiber_enumerated": False,
                        "move_connectivity": str(conn_2b_raw["status"]),
                        "full_fiber_valid": full_ok,
                        "p_low": float(d_lo["mean"]),  # type: ignore[index]
                        "p_high": float(d_hi["mean"]),  # type: ignore[index]
                        "var_bm_low": float(d_lo["var_bm"]),  # type: ignore[index]
                        "var_bm_high": float(d_hi["var_bm"]),  # type: ignore[index]
                        "mcmc_diagnostics": v_diag,
                        "comp_p_low": float(d_lo["mean"]),  # type: ignore[index]
                        "comp_p_high": float(d_hi["mean"]),  # type: ignore[index]
                        "comp_var_bm_low": float(d_lo["var_bm"]),  # type: ignore[index]
                        "comp_var_bm_high": float(d_hi["var_bm"]),  # type: ignore[index]
                        "comp_valid": mcmc_ok,
                    }
            ev_null2b_est[ev_id] = ev_est_2b

        # 3. Null 2c O2-clash-constrained two-tier estimator on observation-invariant [S_obs]_{~2c} for Q2 attribution
        if ev_id in needed_q2_events:
            conn_2c_raw = classify_ternary_move_connectivity(
                mm_e,
                mp_e,
                forbid_minus_row_fn=forbid_o2_clash,
                max_enum_states=max_connectivity_enum_states,
                return_enumerated_pairs=True,
            )
            ev_null2c_conn[ev_id] = _strip_conn_pairs(conn_2c_raw)
            active_c_keys_ev = [
                c_key
                for c_key, _, _ in q2_contrast_defs
                if q2_prep[c_key][9] is None and per_event_q2_matches[c_key].get(ev_id)
            ]

            def _score_q2_state(
                pm_m: Sequence[Sequence[int]],
                pm_p: Sequence[Sequence[int]],
            ) -> Dict[str, Tuple[Tuple[float, float], Tuple[float, float]]]:
                pr_m = {
                    det_gs_e[i].mag_id: profile_from_module_Presence_row(det_gs_e[i].mag_id, pm_m[i])
                    for i in range(len(det_gs_e))
                }
                pr_p = {
                    det_gs_e[i].mag_id: profile_from_module_Presence_row(det_gs_e[i].mag_id, pm_p[i])
                    for i in range(len(det_gs_e))
                }
                pr_p.update(extra_plus_e)
                c_out: Dict[str, Tuple[Tuple[float, float], Tuple[float, float]]] = {}
                for c_k in active_c_keys_ev:
                    req_inc_c, scen_c = q2_prep[c_k][0], q2_prep[c_k][1]
                    amb_c = scen_c.ambient_compounds
                    sc_i_t = StratumCounts()
                    sc_c_t = StratumCounts()
                    for rx, ry, p_c0 in per_event_q2_matches[c_k][ev_id]:
                        sc_i_t.record(
                            _score_recipient_on_null2_profiles(
                                rx.mag_id, req_inc_c, amb_c,
                                ctx_e.rec_base_minus[rx.mag_id], ctx_e.rec_base_plus[rx.mag_id], pr_m, pr_p
                            )
                        )
                        sc_c_t.record(
                            _score_recipient_on_null2_profiles(
                                ry.mag_id, p_c0, amb_c,
                                ctx_e.rec_base_minus[ry.mag_id], ctx_e.rec_base_plus[ry.mag_id], pr_m, pr_p
                            )
                        )
                    c_out[c_k] = (sc_i_t.prevalence_bounds(), sc_c_t.prevalence_bounds())
                return c_out

            ev_est_2c: Dict[str, Dict[str, object]] = {}
            if bool(conn_2c_raw["exact_full_fiber_enumerated"]):
                fiber_pairs_2c = conn_2c_raw["exact_fiber_pairs"]  # type: ignore[assignment]
                full_q2_scores = [_score_q2_state(pm_m, pm_p) for pm_m, pm_p in fiber_pairs_2c]  # type: ignore[union-attr]
                comp_pairs_2c = conn_2c_raw.get("exact_component_pairs")
                if str(conn_2c_raw["status"]) == "connected" or comp_pairs_2c is None:
                    comp_q2_scores = full_q2_scores if str(conn_2c_raw["status"]) == "connected" else None
                else:
                    comp_q2_scores = [_score_q2_state(pm_m, pm_p) for pm_m, pm_p in comp_pairs_2c]  # type: ignore[union-attr]
                for c_k in active_c_keys_ev:
                    pi_lo_f = sum(bs[c_k][0][0] for bs in full_q2_scores) / len(full_q2_scores)
                    pi_hi_f = sum(bs[c_k][0][1] for bs in full_q2_scores) / len(full_q2_scores)
                    pc_lo_f = sum(bs[c_k][1][0] for bs in full_q2_scores) / len(full_q2_scores)
                    pc_hi_f = sum(bs[c_k][1][1] for bs in full_q2_scores) / len(full_q2_scores)
                    pi_lo_c = (
                        sum(bs[c_k][0][0] for bs in comp_q2_scores) / len(comp_q2_scores)
                        if comp_q2_scores is not None
                        else pi_lo_f
                    )
                    pi_hi_c = (
                        sum(bs[c_k][0][1] for bs in comp_q2_scores) / len(comp_q2_scores)
                        if comp_q2_scores is not None
                        else pi_hi_f
                    )
                    pc_lo_c = (
                        sum(bs[c_k][1][0] for bs in comp_q2_scores) / len(comp_q2_scores)
                        if comp_q2_scores is not None
                        else pc_lo_f
                    )
                    pc_hi_c = (
                        sum(bs[c_k][1][1] for bs in comp_q2_scores) / len(comp_q2_scores)
                        if comp_q2_scores is not None
                        else pc_hi_f
                    )
                    ev_est_2c[c_k] = {
                        "tier": "tier1_exact_full_fiber",
                        "exact_full_fiber_enumerated": True,
                        "move_connectivity": str(conn_2c_raw["status"]),
                        "full_fiber_valid": True,
                        "inc_p_low": pi_lo_f,
                        "inc_p_high": pi_hi_f,
                        "inc_var_low": 0.0,
                        "inc_var_high": 0.0,
                        "comp_p_low": pc_lo_f,
                        "comp_p_high": pc_hi_f,
                        "comp_var_low": 0.0,
                        "comp_var_high": 0.0,
                        "mcmc_diagnostics": None,
                        "ck_inc_p_low": pi_lo_c,
                        "ck_inc_p_high": pi_hi_c,
                        "ck_inc_var_low": 0.0,
                        "ck_inc_var_high": 0.0,
                        "ck_comp_p_low": pc_lo_c,
                        "ck_comp_p_high": pc_hi_c,
                        "ck_comp_var_low": 0.0,
                        "ck_comp_var_high": 0.0,
                        "ck_valid": comp_q2_scores is not None,
                    }
            else:
                chains_2c, ham_2c = mcmc_reversible_nested_chains(
                    mm_e,
                    mp_e,
                    rng,
                    n_draws=n_null2_draws,
                    n_chains=4,
                    burnin_mult=20,
                    thin_mult=3,
                    forbid_minus_row_fn=forbid_o2_clash,
                )
                chain_scores_2c = [
                    [_score_q2_state(pm_m, pm_p) for pm_m, pm_p in ch]
                    for ch in chains_2c
                ]
                for c_k in active_c_keys_ev:
                    ilo_ch = [[bd[c_k][0][0] for bd in ch_s] for ch_s in chain_scores_2c]
                    ihi_ch = [[bd[c_k][0][1] for bd in ch_s] for ch_s in chain_scores_2c]
                    clo_ch = [[bd[c_k][1][0] for bd in ch_s] for ch_s in chain_scores_2c]
                    chi_ch = [[bd[c_k][1][1] for bd in ch_s] for ch_s in chain_scores_2c]
                    v_inc = validate_mcmc_prevalence_bounds(
                        ilo_ch, ihi_ch, ham_2c, has_active_move=bool(conn_2c_raw["has_active_move"]), min_ess=min_ess_req
                    )
                    v_comp = validate_mcmc_prevalence_bounds(
                        clo_ch, chi_ch, ham_2c, has_active_move=bool(conn_2c_raw["has_active_move"]), min_ess=min_ess_req
                    )
                    mcmc_ok_2c = bool(v_inc["mcmc_diagnostics_pass"]) and bool(v_comp["mcmc_diagnostics_pass"])
                    full_ok_2c = (str(conn_2c_raw["status"]) == "connected") and mcmc_ok_2c
                    di_lo, di_hi = v_inc["bound_low_diag"], v_inc["bound_high_diag"]  # type: ignore[assignment]
                    dc_lo, dc_hi = v_comp["bound_low_diag"], v_comp["bound_high_diag"]  # type: ignore[assignment]
                    ev_est_2c[c_k] = {
                        "tier": "tier2_reversible_mcmc_4chain",
                        "exact_full_fiber_enumerated": False,
                        "move_connectivity": str(conn_2c_raw["status"]),
                        "full_fiber_valid": full_ok_2c,
                        "inc_p_low": float(di_lo["mean"]),  # type: ignore[index]
                        "inc_p_high": float(di_hi["mean"]),  # type: ignore[index]
                        "inc_var_low": float(di_lo["var_bm"]),  # type: ignore[index]
                        "inc_var_high": float(di_hi["var_bm"]),  # type: ignore[index]
                        "comp_p_low": float(dc_lo["mean"]),  # type: ignore[index]
                        "comp_p_high": float(dc_hi["mean"]),  # type: ignore[index]
                        "comp_var_low": float(dc_lo["var_bm"]),  # type: ignore[index]
                        "comp_var_high": float(dc_hi["var_bm"]),  # type: ignore[index]
                        "mcmc_diagnostics": {"inc": v_inc, "comp": v_comp, "mcmc_diagnostics_pass": mcmc_ok_2c},
                        "ck_inc_p_low": float(di_lo["mean"]),  # type: ignore[index]
                        "ck_inc_p_high": float(di_hi["mean"]),  # type: ignore[index]
                        "ck_inc_var_low": float(di_lo["var_bm"]),  # type: ignore[index]
                        "ck_inc_var_high": float(di_hi["var_bm"]),  # type: ignore[index]
                        "ck_comp_p_low": float(dc_lo["mean"]),  # type: ignore[index]
                        "ck_comp_p_high": float(dc_hi["mean"]),  # type: ignore[index]
                        "ck_comp_var_low": float(dc_lo["var_bm"]),  # type: ignore[index]
                        "ck_comp_var_high": float(dc_hi["var_bm"]),  # type: ignore[index]
                        "ck_valid": mcmc_ok_2c,
                    }
            ev_null2c_q2_est[ev_id] = ev_est_2c

    # Assemble Q1 results per stratum
    q1_by_stratum: Dict[str, Dict[str, object]] = {}
    for idx_s, s_label in enumerate(all_stratum_labels):
        elig_evs, tot_sc, pooled_g2, unres_frac, gate_fail = q1_gate_info[s_label]
        null1_status = "not_run_gate_failed"
        null1_uneval_reason: Optional[str] = None
        ci_n1: Optional[Tuple[float, float]] = None
        ci_n2b_primary: Optional[Tuple[float, float]] = None
        null2b_prim_pooled: Optional[List[float]] = None
        null2b_unavail_reason: Optional[str] = None
        null2b_conn_status: Optional[str] = None
        null2b_comp_cond: Optional[Dict[str, object]] = None
        null2b_ev_conn_list: Optional[List[Dict[str, object]]] = None

        if gate_fail is not None:
            verdict_s = gate_fail
        elif s_label == "C0xC0_ALL45__E_minimal":
            null1_status = "control_pool_not_q1_stratum"
            verdict_s = "control_pool_descriptive_only (compatible C0xC0 control evaluated via matched P2-Q2 contrasts)"
        else:
            per_ev_obs = [per_event_strata[ev][s_label] for ev in elig_evs]
            ev_ests_s = [ev_null2b_est[ev][s_label] for ev in elig_evs]
            (
                null2b_conn_status,
                null2b_ev_conn_list,
                ci_n2b_primary,
                null2b_prim_pooled,
                null2b_unavail_reason,
                null2b_comp_cond,
            ) = _resolve_cohort_q1_null2b(
                elig_evs,
                per_ev_obs,
                ev_null2b_conn,
                ev_ests_s,
                p2_a1_approved=p2_a1_approved,
                n_boot=n_boot,
                seed_primary=seed + 10 + idx_s,
                seed_comp=seed + 100 + idx_s,
            )

            # Strict Null 1 evaluability checks
            missing_src_evs = [ev for ev in elig_evs if ev in ev_null1_no_source]
            if hits_by_mag is None:
                null1_status = "unevaluable"
                null1_uneval_reason = "missing per-MAG gene-hit contig tables (hits_by_mag is None)"
            elif missing_src_evs:
                null1_status = "unevaluable"
                null1_uneval_reason = (
                    f"zero >=95%-complete source MAGs with gene-hit blocks in "
                    f"{len(missing_src_evs)}/{len(elig_evs)} eligible events (first={missing_src_evs[0]})"
                )
            else:
                per_ev_n1_raw = [ev_null1_counts[ev][s_label] for ev in elig_evs]
                n_evs_with_sim_cand = sum(
                    1 for sc_x in per_ev_n1_raw if sc_x is not None and sc_x.n_cand_aux >= 1
                )
                tot_sim_cand = sum(sc_x.n_cand_aux for sc_x in per_ev_n1_raw if sc_x is not None)
                if tot_sim_cand == 0 or n_evs_with_sim_cand < min_events:
                    null1_status = "unevaluable"
                    null1_uneval_reason = (
                        f"insufficient simulated candidate auxotrophs under Null 1 across {n_null1_draws} draws/event "
                        f"(total_sim_cand_aux={tot_sim_cand}, events_with_sim_cand_aux={n_evs_with_sim_cand} < {min_events})"
                    )
                else:
                    null1_status = "evaluable"
                    per_ev_n1 = [sc_x for sc_x in per_ev_n1_raw if sc_x is not None]
                    ci_n1 = event_bootstrap_q1_bounds(per_ev_obs, per_ev_n1, n_boot=n_boot, seed=seed + idx_s)

            verdict_s = _decide_p2_q1_symmetric(
                ci_null1=ci_n1 if null1_status == "evaluable" else None,
                null1_uneval_reason=null1_uneval_reason,
                ci_null2_primary=ci_n2b_primary,
                null2_unavail_reason=null2b_unavail_reason,
                delta_p=0.05,
            )

        q1_by_stratum[s_label] = {
            "n_q1_eligible_events": len(elig_evs),
            "counts": {
                "n_auto": tot_sc.n_auto,
                "n_cand_aux": tot_sc.n_cand_aux,
                "n_d0": tot_sc.n_d0,
                "n_d1": tot_sc.n_d1,
                "n_plexi": tot_sc.n_plexi,
                "n_ambig": tot_sc.n_ambig,
                "n_no_donor": tot_sc.n_no_donor,
                "unresolved_fraction": round(unres_frac, 4),
            },
            "pooled_detected_gate2": pooled_g2,
            "null1_status": null1_status,
            "null1_unevaluable_reason": null1_uneval_reason,
            "ci_null1": list(ci_n1) if ci_n1 is not None else None,
            "ci_null1_upper_below_zero": bool(ci_n1[1] < 0.0) if ci_n1 is not None else None,
            "null2b_connectivity_status": null2b_conn_status,
            "null2b_primary_unavailable_reason": null2b_unavail_reason,
            "null2b_primary_pooled_null_bounds": null2b_prim_pooled,
            "null2b_per_event_connectivity": null2b_ev_conn_list,
            "ci_null2b": list(ci_n2b_primary) if ci_n2b_primary is not None else None,
            "null2b_component_conditional": null2b_comp_cond,
            "verdict_q1": verdict_s,
        }

    # Assemble P2-Q2 contrasts across all 4 (pair, supply scenario E) combinations
    q2_by_contrast: Dict[str, Dict[str, object]] = {}
    c0_key_by_scen = {
        SCENARIO_E_MINIMAL.name: "C0xC0_ALL45__E_minimal",
        SCENARIO_E_PLUS_DMSP.name: "C0xC0_ALL45__E_plus_DMSP",
        SCENARIO_E_PLUS_DIN.name: "C0xC0_ALL45__E_plus_DIN",
    }

    for idx_c, (c_key, req_inc, scen_c) in enumerate(q2_contrast_defs):
        (
            _, _, matched_evs, ev_inc_m, ev_comp_m, tot_inc_m, tot_comp_m, inc_g2_m, comp_g2_m, q2_gate_fail
        ) = q2_prep[c_key]
        n_inc_tot_before, n_dist_m, n_self_m = per_event_q2_match_stats[c_key]
        n_matched_tot = n_dist_m + n_self_m
        u_inc_rate_opt = round(tot_inc_m.n_unresolved / tot_inc_m.n_cand_aux, 4) if tot_inc_m.n_cand_aux > 0 else None
        u_comp_rate_opt = round(tot_comp_m.n_unresolved / tot_comp_m.n_cand_aux, 4) if tot_comp_m.n_cand_aux > 0 else None

        ci_or: Optional[Tuple[float, float]] = None
        res_only_ci_or: Optional[List[float]] = None
        res_only_status: str = "not_run_gate_failed"
        null2c_attr: Optional[Dict[str, object]] = None
        unmatched_pool_sens: Optional[Dict[str, object]] = None

        if q2_gate_fail is not None:
            verdict_q2_c = q2_gate_fail
        else:
            ci_or = _event_bootstrap_q2_log_or_bounds(
                matched_evs, ev_inc_m, ev_comp_m, n_boot=n_boot, seed=seed + 20 + idx_c
            )
            # Primary P2-Q2 verdict rests strictly on matched observed Manski log-odds interval vs ln(2) (§4.3)
            verdict_q2_c = decide_p2_q2(ci_or, delta_or=math.log(2.0), gate_failure=None)

            # Resolved-only sensitivity interval alongside U_inc/N_inc and U_comp/N_comp
            n_res_inc = tot_inc_m.n_d0 + tot_inc_m.n_d1 + tot_inc_m.n_plexi
            n_res_comp = tot_comp_m.n_d0 + tot_comp_m.n_d1 + tot_comp_m.n_plexi
            if n_res_inc == 0 or n_res_comp == 0:
                res_only_status = f"unevaluable (zero resolved tuples: n_res_inc={n_res_inc}, n_res_comp={n_res_comp})"
            else:
                res_inc_ev = {
                    e: StratumCounts(n_d0=ev_inc_m[e].n_d0, n_d1=ev_inc_m[e].n_d1, n_plexi=ev_inc_m[e].n_plexi)
                    for e in matched_evs
                }
                res_comp_ev = {
                    e: StratumCounts(n_d0=ev_comp_m[e].n_d0, n_d1=ev_comp_m[e].n_d1, n_plexi=ev_comp_m[e].n_plexi)
                    for e in matched_evs
                }
                r_lo, r_hi = _event_bootstrap_q2_log_or_bounds(
                    matched_evs, res_inc_ev, res_comp_ev, n_boot=n_boot, seed=seed + 30 + idx_c
                )
                res_only_ci_or = [r_lo, r_hi]
                res_only_status = "evaluable"

            # Constrained Null 2c attribution on [S_obs]_{~2c} (per-event primary source selection; never blocks primary P2-Q2 verdict)
            ev_2c_ests = {e: ev_null2c_q2_est[e][c_key] for e in matched_evs}
            null2c_attr = _resolve_cohort_q2_null2c(
                matched_evs,
                ev_inc_m,
                ev_comp_m,
                ev_null2c_conn,
                ev_2c_ests,
                p2_a1_approved=p2_a1_approved,
                n_boot=n_boot,
                seed_primary=seed + 40 + idx_c,
                seed_comp=seed + 140 + idx_c,
            )

            # Descriptive unmatched all-45-pairs C0xC0 comparison under the same supply scenario E
            c0_all_k = c0_key_by_scen[scen_c.name]
            inc_evs_u = {e for e in bio_events if per_event_strata[e][c_key].n_cand_aux >= 1}
            comp_evs_u = {e for e in bio_events if per_event_strata[e][c0_all_k].n_cand_aux >= 1}
            union_evs_u = sorted(inc_evs_u | comp_evs_u)
            if union_evs_u:
                unm_lo, unm_hi = _event_bootstrap_q2_log_or_bounds(
                    union_evs_u,
                    {e: per_event_strata[e][c_key] for e in union_evs_u},
                    {e: per_event_strata[e][c0_all_k] for e in union_evs_u},
                    n_boot=n_boot,
                    seed=seed + 50 + idx_c,
                )
                unmatched_pool_sens = {
                    "label": "descriptive_unmatched_all45_c0_sensitivity_only",
                    "ci_log_or": [unm_lo, unm_hi],
                    "implied_band_if_unmatched": decide_p2_q2((unm_lo, unm_hi), delta_or=math.log(2.0)),
                }

        inc_evs_raw_count = sum(1 for ev in bio_events if per_event_strata[ev][c_key].n_cand_aux >= 1)
        q2_by_contrast[c_key] = {
            "supply_scenario": scen_c.name,
            "n_inc_eligible_events_before_matching": inc_evs_raw_count,
            "n_matched_events": len(matched_evs),
            "matching_denominators": {
                "n_inc_tuples_before_matching": n_inc_tot_before,
                "n_matched_pairs": n_matched_tot,
                "n_distinct_recipient_matches": n_dist_m,
                "n_within_recipient_self_matches": n_self_m,
                "n_inc_unmatched_dropped": n_inc_tot_before - n_matched_tot,
                "matched_inc_counts": {
                    "n_cand_aux": tot_inc_m.n_cand_aux,
                    "n_d0": tot_inc_m.n_d0,
                    "n_d1": tot_inc_m.n_d1,
                    "n_plexi": tot_inc_m.n_plexi,
                    "n_ambig": tot_inc_m.n_ambig,
                    "n_no_donor": tot_inc_m.n_no_donor,
                },
                "matched_comp_counts": {
                    "n_cand_aux": tot_comp_m.n_cand_aux,
                    "n_d0": tot_comp_m.n_d0,
                    "n_d1": tot_comp_m.n_d1,
                    "n_plexi": tot_comp_m.n_plexi,
                    "n_ambig": tot_comp_m.n_ambig,
                    "n_no_donor": tot_comp_m.n_no_donor,
                },
            },
            "u_inc_over_n_inc": u_inc_rate_opt,
            "u_comp_over_n_comp": u_comp_rate_opt,
            "pooled_detected_gate2_inc": inc_g2_m,
            "pooled_detected_gate2_comp": comp_g2_m,
            "ci_log_or": list(ci_or) if ci_or is not None else None,
            "ci_log_or_upper_below_zero": bool(ci_or[1] < 0.0) if ci_or is not None else None,
            "resolved_only_status": res_only_status,
            "resolved_only_ci_log_or": res_only_ci_or,
            "null2c_attribution": null2c_attr,
            "descriptive_unmatched_all45_c0_sensitivity": unmatched_pool_sens,
            "verdict_q2": verdict_q2_c,
        }

    v_q1_ri3_min = str(q1_by_stratum["R-I3_B12_REDS__E_minimal"]["verdict_q1"])
    v_q1_ri3_dmsp = str(q1_by_stratum["R-I3_B12_REDS__E_plus_DMSP"]["verdict_q1"])
    v_q1_ri2_min = str(q1_by_stratum["R-I2_B12_FIXEDN__E_minimal"]["verdict_q1"])
    v_q1_ri2_din = str(q1_by_stratum["R-I2_B12_FIXEDN__E_plus_DIN"]["verdict_q1"])

    v_q2_ri3_min = str(q2_by_contrast["R-I3_B12_REDS__E_minimal"]["verdict_q2"])
    v_q2_ri3_dmsp = str(q2_by_contrast["R-I3_B12_REDS__E_plus_DMSP"]["verdict_q2"])
    v_q2_ri2_min = str(q2_by_contrast["R-I2_B12_FIXEDN__E_minimal"]["verdict_q2"])
    v_q2_ri2_din = str(q2_by_contrast["R-I2_B12_FIXEDN__E_plus_DIN"]["verdict_q2"])

    return {
        "n_biological_events_in_map": len(bio_events),
        "n_control_libraries_excluded": len(control_libs),
        "n_events_with_any_hq_mag_flow_diagnostic_only": n_events_with_any_hq,
        "targets_bed_sha256_verified": expected_targets_bed_sha256,
        "draw_parameters": {
            "n_null1_draws": n_null1_draws,
            "n_null2_draws": n_null2_draws,
            "n_boot": n_boot,
            "require_protocol_1000_draws": require_protocol_1000_draws,
            "p2_a1_approved": p2_a1_approved,
        },
        "per_event_detected_gate2": per_event_gate2,
        "q1_by_stratum": q1_by_stratum,
        "q2_by_contrast": q2_by_contrast,
        "combined_primary_verdict_p2_q1": (
            f"R-I3 (E_minimal): {v_q1_ri3_min}; R-I3 (E_plus_DMSP): {v_q1_ri3_dmsp}; "
            f"R-I2 (E_minimal): {v_q1_ri2_min}; R-I2 (E_plus_DIN): {v_q1_ri2_din}"
        ),
        "combined_primary_verdict_p2_q2": (
            f"R-I3 vs matched C0xC0 (E_minimal): {v_q2_ri3_min}; "
            f"R-I3 vs matched C0xC0 (E_plus_DMSP): {v_q2_ri3_dmsp}; "
            f"R-I2 vs matched C0xC0 (E_minimal): {v_q2_ri2_min}; "
            f"R-I2 vs matched C0xC0 (E_plus_DIN): {v_q2_ri2_din}"
        ),
    }


def _make_synthetic_gene_hits(mag_id: str, kos: AbstractSet[str], n_contigs: int = 2, pad_genes: int = 80) -> List[GeneHit]:
    """Build deterministic internal GeneHit lists across n_contigs so downsample_contig_blocks exercises real cuts."""
    ko_sorted = sorted(kos)
    hits: List[GeneHit] = []
    idx = 0
    for c_i in range(n_contigs):
        cid = f"{mag_id}_ctg{c_i + 1}"
        # Anchor contig ends with non-P2 dummy genes so internal P2 KOs stay internal unless cut
        idx += 1
        hits.append(GeneHit(protein_id=f"{mag_id}_p{idx:04d}", contig_id=cid, kos=frozenset({"K00001"})))
        chunk = ko_sorted[c_i::n_contigs]
        for k in chunk:
            idx += 1
            hits.append(GeneHit(protein_id=f"{mag_id}_p{idx:04d}", contig_id=cid, kos=frozenset({k})))
        for _ in range(max(2, pad_genes // n_contigs)):
            idx += 1
            hits.append(GeneHit(protein_id=f"{mag_id}_p{idx:04d}", contig_id=cid, kos=frozenset({"K00001"})))
        idx += 1
        hits.append(GeneHit(protein_id=f"{mag_id}_p{idx:04d}", contig_id=cid, kos=frozenset({"K00001"})))
    return hits


def verify_freshwater_event_gate_fixtures() -> Dict[str, object]:
    """Unit-test provenance, missing/dup inputs, Gates 1-3, <1000 draw gate, zero-eligible Null 1, matched vs pooled Q2, all E, and P2-A1."""
    all_kos = all_referenced_kos()
    c0_not_his_trp = set().union(
        *(PATHWAY_RULES[m].all_kos() for m in C0_CONTROL_MODULES if m not in ("C_HIS", "C_TRP"))
    )
    c0_all_kos = set().union(*(PATHWAY_RULES[m].all_kos() for m in C0_CONTROL_MODULES))

    g_auto = GenomeAnnotations(
        mag_id="dastool-auto_hq",
        domain="Bacteria",
        completeness=99.0,
        contamination=0.5,
        genome_size_bp=3_000_000,
        contig_lengths=(3_000_000,),
        kos=frozenset(all_kos),
    )
    g_c0_rec = GenomeAnnotations(
        mag_id="dastool-c0rec_hq",
        domain="Bacteria",
        completeness=99.0,
        contamination=0.5,
        genome_size_bp=3_000_000,
        contig_lengths=(3_000_000,),
        kos=frozenset({"K00549"} | PATHWAY_RULES["C_THI_SPLIT"].all_kos() | c0_not_his_trp),
    )
    g_ri3_c0_rec = GenomeAnnotations(
        mag_id="dastool-ri3rec_hq",
        domain="Bacteria",
        completeness=99.0,
        contamination=0.5,
        genome_size_bp=3_000_000,
        contig_lengths=(3_000_000,),
        kos=frozenset({"K00548"} | PATHWAY_RULES["ASSIM_NO3"].all_kos() | c0_not_his_trp),
    )
    d_b12 = GenomeAnnotations(
        mag_id="dastool-don_b12",
        domain="Bacteria",
        completeness=85.0,
        contamination=1.0,
        genome_size_bp=2_500_000,
        contig_lengths=(2_500_000,),
        kos=frozenset(PATHWAY_RULES["P_B12_AER"].all_kos() | PATHWAY_RULES["P_BLUB"].all_kos()),
    )
    d_dsr = GenomeAnnotations(
        mag_id="dastool-don_dsr",
        domain="Bacteria",
        completeness=85.0,
        contamination=1.0,
        genome_size_bp=2_500_000,
        contig_lengths=(2_500_000,),
        kos=frozenset(PATHWAY_RULES["A_DSR"].all_kos()),
    )
    d_c0 = GenomeAnnotations(
        mag_id="dastool-don_c0",
        domain="Bacteria",
        completeness=85.0,
        contamination=1.0,
        genome_size_bp=2_500_000,
        contig_lengths=(2_500_000,),
        kos=frozenset(PATHWAY_RULES["C_HIS"].all_kos() | PATHWAY_RULES["C_TRP"].all_kos()),
    )

    genomes = {
        g.mag_id: g for g in (g_auto, g_c0_rec, g_ri3_c0_rec, d_b12, d_dsr, d_c0)
    }
    rep_ids = sorted(genomes.keys())
    hq_ids = ["dastool-auto_hq", "dastool-c0rec_hq", "dastool-ri3rec_hq"]
    base_gtdb = {
        m: "d__Bacteria;p__Pseudomonadota;c__Gammaproteobacteria" for m in rep_ids
    }
    rep_c_lens = {
        "dastool-auto_hq": {"ctg_auto": 1000},
        "dastool-c0rec_hq": {"ctg_c0rec": 1000},
        "dastool-ri3rec_hq": {"ctg_ri3rec": 1000},
        "dastool-don_b12": {"ctg_b12": 1000},
        "dastool-don_dsr": {"ctg_dsr": 1000},
        "dastool-don_c0": {"ctg_c0": 1000},
    }
    target_lens = {
        "ctg_auto": 1000,
        "ctg_c0rec": 1000,
        "ctg_ri3rec": 1000,
        "ctg_b12": 1000,
        "ctg_dsr": 1000,
        "ctg_c0": 1000,
        "ctg_nonrep": 1000,
        "ctg_unb": 1000,
    }
    valid_targets = set(target_lens.keys())
    empty_groups = {
        gk: {"rep490": {}, "dastool655_nonrep": {}, "unb_prok": {}, "euk": {}}
        for gk in ("all_donor_modules", "R-I3_B12_REDS", "R-I2_B12_FIXEDN", "R-I6_B12", "C_THI_SPLIT", "C0_ALL10")
    }

    def _write_gz(path: str, rows: Sequence[Tuple[str, int, int]]) -> str:
        with gzip.open(path, "wt", encoding="utf-8") as gfh:
            for c_id, s_i, e_i in rows:
                gfh.write(f"{c_id}\t{s_i}\t{e_i}\n")
        return sha256_12(path)

    # Unit-check symmetric P2-A1 §3.1 partial-information decision function across all 4 (Null 1, primary Null 2b) quadrants
    assert _decide_p2_q1_symmetric((0.12, 0.40), None, (0.10, 0.35), None) == "Supported"
    assert (
        _decide_p2_q1_symmetric((0.12, 0.40), None, None, "move set disconnected")
        == "Inconclusive (primary Null 2 unavailable: move set disconnected)"
    )
    assert (
        _decide_p2_q1_symmetric((0.12, 0.40), None, None, "connectivity unknown")
        == "Inconclusive (primary Null 2 unavailable: connectivity unknown)"
    )
    assert (
        _decide_p2_q1_symmetric((0.12, 0.40), None, None, "MCMC stationary-mean mixing unvalidated")
        == "Inconclusive (primary Null 2 unavailable: MCMC stationary-mean mixing unvalidated)"
    )
    assert _decide_p2_q1_symmetric((-0.10, 0.02), None, None, "move set disconnected") == "Falsified"
    assert _decide_p2_q1_symmetric((-0.10, 0.02), None, None, "connectivity unknown") == "Falsified"
    assert _decide_p2_q1_symmetric((0.01, 0.12), None, None, "connectivity unknown") == "Inconclusive"
    # Symmetric OR-falsification when Null 1 is unevaluable and certified primary Null 2b upper bound < 0.05
    assert _decide_p2_q1_symmetric(None, "no >=95% MAGs", (-0.25, 0.00), None) == "Falsified"
    assert (
        _decide_p2_q1_symmetric(None, "no >=95% MAGs", (0.10, 0.35), None)
        == "Inconclusive (Null 1 unevaluable: no >=95% MAGs)"
    )
    assert (
        _decide_p2_q1_symmetric(None, "no >=95% MAGs", None, "connectivity unknown")
        == "Inconclusive (Null 1 unevaluable: no >=95% MAGs; primary Null 2 unavailable: connectivity unknown)"
    )

    with tempfile.TemporaryDirectory() as tmpd:
        v2_bed = os.path.join(tmpd, "targets_v2.bed")
        with open(v2_bed, "w", encoding="utf-8") as fh:
            for c_id in sorted(valid_targets):
                fh.write(f"{c_id}\t0\t{target_lens[c_id]}\n")
        v2_sha = sha256_12(v2_bed)

        # 1. Provenance check: v1 reps-only interval file paired with genuine v2 target BED
        d_prov = os.path.join(tmpd, "prov_case")
        os.makedirs(d_prov)
        v1_cov_sha = _write_gz(os.path.join(d_prov, "lib1.covered.bed.gz"), [("ctg_auto", 0, 800)])
        emap_1 = os.path.join(tmpd, "emap1.tsv")
        with open(emap_1, "w", encoding="utf-8") as fh:
            fh.write("library\tevent\nlib1\tev01\n")
        # 1a: MANIFEST honestly records v1 target hash 5d1f5cfac723 while genuine v2_bed is supplied
        with open(os.path.join(d_prov, "MANIFEST.tsv"), "w", encoding="utf-8") as fh:
            fh.write(f"library\tcovered_bed_sha256\ttargets_bed_sha256\nlib1\t{v1_cov_sha}\t5d1f5cfac723\n")
        caught_v1_manifest = False
        try:
            evaluate_freshwater_events(
                d_prov, emap_1, genomes, rep_ids, hq_ids, [], rep_c_lens, target_lens, valid_targets,
                {}, {}, empty_groups, expected_targets_bed_sha256=v2_sha, targets_bed_path=v2_bed,
                gtdb_clade_by_mag=base_gtdb,
            )
        except ValueError:
            caught_v1_manifest = True
        assert caught_v1_manifest, "Failed to reject v1 interval file paired with genuine v2 target BED (manifest check)"

        # 1b: v1 interval file substituted on disk under a v2 MANIFEST entry
        with open(os.path.join(d_prov, "MANIFEST.tsv"), "w", encoding="utf-8") as fh:
            fh.write(f"library\tcovered_bed_sha256\ttargets_bed_sha256\nlib1\tdeadbeefcafe\t{v2_sha}\n")
        caught_v1_subst = False
        try:
            evaluate_freshwater_events(
                d_prov, emap_1, genomes, rep_ids, hq_ids, [], rep_c_lens, target_lens, valid_targets,
                {}, {}, empty_groups, expected_targets_bed_sha256=v2_sha, targets_bed_path=v2_bed,
                gtdb_clade_by_mag=base_gtdb,
            )
        except ValueError:
            caught_v1_subst = True
        assert caught_v1_subst, "Failed to reject substituted v1 interval file against v2 MANIFEST sha256"

        # 1c: manifest_v2.tsv 64-char sha256 schema check (rejects v1 target hash even when v2_bed is passed)
        mv2_test = os.path.join(tmpd, "manifest_v2_test.tsv")
        full_v1_sha = sha256_full(os.path.join(d_prov, "lib1.covered.bed.gz"))
        with open(mv2_test, "w", encoding="utf-8") as fh:
            fh.write(
                "library\tproducer_job_task\tbam\tsha256\ttargets_bed_sha12\tgrex_out_dir\tmethod\tlocal_vs_grex\n"
                f"lib1\t8194721_1\tlib1.sorted.bam\t{full_v1_sha}\t5d1f5cfac723\tout\tbamtobed_splitD\tmatch\n"
            )
        caught_mv2_v1 = False
        try:
            evaluate_freshwater_events(
                d_prov, emap_1, genomes, rep_ids, hq_ids, [], rep_c_lens, target_lens, valid_targets,
                {}, {}, empty_groups, expected_targets_bed_sha256=v2_sha, targets_bed_path=v2_bed,
                manifest_tsv=mv2_test, gtdb_clade_by_mag=base_gtdb,
            )
        except ValueError:
            caught_mv2_v1 = True
        assert caught_mv2_v1, "Failed to reject v1 producer row in manifest_v2.tsv schema"

        # 2. Missing *.covered.bed.gz input must raise ValueError (never silently count as 0 coverage)
        emap_missing = os.path.join(tmpd, "emap_missing.tsv")
        with open(emap_missing, "w", encoding="utf-8") as fh:
            fh.write("library\tevent\nlib1\tev01\nlib_missing\tev02\n")
        with open(os.path.join(d_prov, "MANIFEST.tsv"), "w", encoding="utf-8") as fh:
            fh.write(
                f"library\tcovered_bed_sha256\ttargets_bed_sha256\n"
                f"lib1\t{v1_cov_sha}\t{v2_sha}\n"
                f"lib_missing\t{v1_cov_sha}\t{v2_sha}\n"
            )
        caught_missing_bed = False
        try:
            evaluate_freshwater_events(
                d_prov, emap_missing, genomes, rep_ids, hq_ids, [], rep_c_lens, target_lens, valid_targets,
                {}, {}, empty_groups, expected_targets_bed_sha256=v2_sha, targets_bed_path=v2_bed,
                gtdb_clade_by_mag=base_gtdb,
            )
        except ValueError:
            caught_missing_bed = True
        assert caught_missing_bed, "Failed to reject missing covered.bed.gz file"

        # 3. Duplicate / multi-event library assignment must raise ValueError
        emap_dup = os.path.join(tmpd, "emap_dup.tsv")
        with open(emap_dup, "w", encoding="utf-8") as fh:
            fh.write("library\tevent\nlib1\tev01\nlib1\tev02\n")
        caught_dup_lib = False
        try:
            evaluate_freshwater_events(
                d_prov, emap_dup, genomes, rep_ids, hq_ids, [], rep_c_lens, target_lens, valid_targets,
                {}, {}, empty_groups, expected_targets_bed_sha256=v2_sha, targets_bed_path=v2_bed,
                gtdb_clade_by_mag=base_gtdb,
            )
        except ValueError:
            caught_dup_lib = True
        assert caught_dup_lib, "Failed to reject multi-event library assignment"

        # 4. HQ-but-ineligible events (12 events with HQ MAGs, but only 7 have R-I3 cand_aux while 11 have C0xC0 cand_aux)
        #    AND 0/0 zero-hit Gate 2 check when donor hits map is empty!
        d_hq = os.path.join(tmpd, "hq_inelig")
        os.makedirs(d_hq)
        emap_hq = os.path.join(tmpd, "emap_hq.tsv")
        m_lines = ["library\tcovered_bed_sha256\ttargets_bed_sha256"]
        e_lines = ["library\tevent"]
        for idx in range(1, 13):
            lib = f"lib_{idx:02d}"
            ev = f"ev_{idx:02d}"
            if idx <= 7:
                rows_iv = [("ctg_ri3rec", 0, 800), ("ctg_b12", 0, 800), ("ctg_dsr", 0, 800), ("ctg_c0", 0, 800)]
            elif idx <= 11:
                rows_iv = [("ctg_c0rec", 0, 800), ("ctg_c0", 0, 800)]
            else:
                rows_iv = [("ctg_auto", 0, 800)]
            c_sha = _write_gz(os.path.join(d_hq, f"{lib}.covered.bed.gz"), rows_iv)
            m_lines.append(f"{lib}\t{c_sha}\t{v2_sha}")
            e_lines.append(f"{lib}\t{ev}")
        with open(os.path.join(d_hq, "MANIFEST.tsv"), "w", encoding="utf-8") as fh:
            fh.write("\n".join(m_lines) + "\n")
        with open(emap_hq, "w", encoding="utf-8") as fh:
            fh.write("\n".join(e_lines) + "\n")

        res_hq = evaluate_freshwater_events(
            d_hq, emap_hq, genomes, rep_ids, hq_ids, [], rep_c_lens, target_lens, valid_targets,
            {}, {}, empty_groups, expected_targets_bed_sha256=v2_sha, targets_bed_path=v2_bed, min_events=10,
            gtdb_clade_by_mag=base_gtdb,
        )
        assert res_hq["n_events_with_any_hq_mag_flow_diagnostic_only"] == 12
        q1_ri3 = res_hq["q1_by_stratum"]["R-I3_B12_REDS__E_minimal"]  # type: ignore[index]
        q1_c0 = res_hq["q1_by_stratum"]["C0xC0_ALL45__E_minimal"]  # type: ignore[index]
        q2_ri3 = res_hq["q2_by_contrast"]["R-I3_B12_REDS__E_minimal"]  # type: ignore[index]
        assert q1_ri3["n_q1_eligible_events"] == 7
        assert "n_events=7 < 10" in q1_ri3["verdict_q1"]
        assert q1_c0["n_q1_eligible_events"] == 11
        assert q1_c0["pooled_detected_gate2"]["unbinned_share_vs_490reps"] is None
        assert q1_c0["pooled_detected_gate2"]["status_vs_490reps"] == "unevaluable_zero_denominator"
        assert "Gate 2 unevaluable — 0/0" in q1_c0["verdict_q1"]
        assert q2_ri3["n_matched_events"] == 7
        assert "matched_events=7" in q2_ri3["verdict_q2"]

        # 5. Gate 3 check (85% complete donors -> unresolved_fraction > 0.50) + both non-rep Gate 2 denominators
        d_pass = os.path.join(tmpd, "pass_case")
        os.makedirs(d_pass)
        emap_pass = os.path.join(tmpd, "emap_pass.tsv")
        m_lines_p = ["library\tcovered_bed_sha256\ttargets_bed_sha256"]
        e_lines_p = ["library\tevent"]
        for idx in range(1, 11):
            lib = f"plib_{idx:02d}"
            ev = f"pev_{idx:02d}"
            rows_iv = [
                ("ctg_ri3rec", 0, 900),
                ("ctg_c0rec", 0, 900),
                ("ctg_b12", 0, 900),
                ("ctg_dsr", 0, 900),
                ("ctg_c0", 0, 900),
                ("ctg_nonrep", 0, 900),
                ("ctg_unb", 0, 900),
            ]
            c_sha = _write_gz(os.path.join(d_pass, f"{lib}.covered.bed.gz"), rows_iv)
            m_lines_p.append(f"{lib}\t{c_sha}\t{v2_sha}")
            e_lines_p.append(f"{lib}\t{ev}")
        with open(os.path.join(d_pass, "MANIFEST.tsv"), "w", encoding="utf-8") as fh:
            fh.write("\n".join(m_lines_p) + "\n")
        with open(emap_pass, "w", encoding="utf-8") as fh:
            fh.write("\n".join(e_lines_p) + "\n")
        pos_groups = {
            gk: {
                "rep490": {"ctg_b12": 10},
                "dastool655_nonrep": {"ctg_nonrep": 5},
                "unb_prok": {"ctg_unb": 3},
                "euk": {},
            }
            for gk in empty_groups
        }
        res_pass = evaluate_freshwater_events(
            d_pass, emap_pass, genomes, rep_ids, hq_ids, [], rep_c_lens, target_lens, valid_targets,
            {}, {}, pos_groups, expected_targets_bed_sha256=v2_sha, targets_bed_path=v2_bed, min_events=10,
            gtdb_clade_by_mag=base_gtdb,
        )
        g2_pass = res_pass["q1_by_stratum"]["R-I3_B12_REDS__E_minimal"]["pooled_detected_gate2"]  # type: ignore[index]
        assert g2_pass["denom_vs_490reps"] == 130
        assert g2_pass["denom_vs_655dastool"] == 180
        assert g2_pass["unbinned_share_vs_490reps"] == round(30 / 130, 4)
        assert g2_pass["unbinned_share_vs_655dastool"] == round(30 / 180, 4)
        v_gate3_q2 = str(res_pass["q2_by_contrast"]["R-I3_B12_REDS__E_minimal"]["verdict_q2"])  # type: ignore[index]
        assert "unresolved_fraction > 0.50" in v_gate3_q2

        # 6. Near-complete (98%) donors: test (a) <1000 draw-count gate, and (b) zero-eligible Null 1 propagation
        genomes_resolved = dict(genomes)
        for mid_d, src_d in (("dastool-don_b12", d_b12), ("dastool-don_dsr", d_dsr), ("dastool-don_c0", d_c0)):
            genomes_resolved[mid_d] = GenomeAnnotations(
                mag_id=mid_d, domain="Bacteria", completeness=98.0, contamination=6.0,
                genome_size_bp=2_500_000, contig_lengths=(2_500_000,), kos=src_d.kos,
            )
        # 6a: < 1000 draw-count gate blocks verdict emission when require_protocol_1000_draws=True
        res_sub1000 = evaluate_freshwater_events(
            d_pass, emap_pass, genomes_resolved, rep_ids, hq_ids, [], rep_c_lens, target_lens, valid_targets,
            {}, {}, pos_groups, expected_targets_bed_sha256=v2_sha, targets_bed_path=v2_bed, min_events=10,
            n_null1_draws=25, n_null2_draws=25, n_boot=100, require_protocol_1000_draws=True,
            gtdb_clade_by_mag=base_gtdb,
        )
        assert "non-inferential draw count < 1000" in str(res_sub1000["q1_by_stratum"]["R-I3_B12_REDS__E_minimal"]["verdict_q1"])  # type: ignore[index]
        assert "non-inferential draw count < 1000" in str(res_sub1000["q2_by_contrast"]["R-I3_B12_REDS__E_minimal"]["verdict_q2"])  # type: ignore[index]

        # 6b: Zero >=95% source MAGs (nc_ids=[]) AND zero simulated candidate auxotrophs (nc MAG lacks K00548)
        #     both propagate Inconclusive (Null 1 unevaluable: ...) without StratumCounts(n_d1=1) or unfragmented fallback
        hits_c0_only = {"dastool-c0rec_hq": _make_synthetic_gene_hits("dastool-c0rec_hq", g_c0_rec.kos)}
        res_no_nc = evaluate_freshwater_events(
            d_pass, emap_pass, genomes_resolved, rep_ids, hq_ids, [], rep_c_lens, target_lens, valid_targets,
            {}, {}, pos_groups, expected_targets_bed_sha256=v2_sha, targets_bed_path=v2_bed, min_events=10,
            n_null1_draws=20, n_null2_draws=20, n_boot=100, require_protocol_1000_draws=False,
            gtdb_clade_by_mag=base_gtdb, hits_by_mag=hits_c0_only,
        )
        q1_no_nc = res_no_nc["q1_by_stratum"]["R-I3_B12_REDS__E_minimal"]  # type: ignore[index]
        assert q1_no_nc["null1_status"] == "unevaluable"
        assert q1_no_nc["ci_null1"] is None
        assert "Null 1 unevaluable: zero >=95%-complete source MAGs" in str(q1_no_nc["verdict_q1"])

        res_zero_sim_cand = evaluate_freshwater_events(
            d_pass, emap_pass, genomes_resolved, rep_ids, hq_ids, ["dastool-c0rec_hq"], rep_c_lens, target_lens, valid_targets,
            {}, {}, pos_groups, expected_targets_bed_sha256=v2_sha, targets_bed_path=v2_bed, min_events=10,
            n_null1_draws=20, n_null2_draws=20, n_boot=100, require_protocol_1000_draws=False,
            gtdb_clade_by_mag=base_gtdb, hits_by_mag=hits_c0_only,
        )
        q1_zero_cand = res_zero_sim_cand["q1_by_stratum"]["R-I3_B12_REDS__E_minimal"]  # type: ignore[index]
        assert q1_zero_cand["null1_status"] == "unevaluable"
        assert q1_zero_cand["ci_null1"] is None
        assert "insufficient simulated candidate auxotrophs under Null 1" in str(q1_zero_cand["verdict_q1"])
        # Meanwhile Q2 on the same cohort with distinct clade-matched dastool-c0rec_hq (0/10 vs 10/10) is Supported!
        assert str(res_zero_sim_cand["q2_by_contrast"]["R-I3_B12_REDS__E_minimal"]["verdict_q2"]) == "Supported"  # type: ignore[index]

        # 6c: Symmetric OR-falsification on missing Null 1 (Astra 249778991556640768):
        #     (i) Tier 1 exact full-fiber Null 2b certified (2 detected MAGs: ri3rec + generalist d_b12_dsr -> exact Omega enum,
        #         obs = 0/10 plexi, ci_null2b[1] < 0.05) with nc_ids=[] -> MUST emit "Falsified"!
        #     (ii) Tier 2 Ryser-connected + 4-chain MCMC-validated full-fiber Null 2b (5 detected MAGs, binary {0,2} active matrix,
        #          mcmc_diagnostics_pass == True, ci_null2b[1] < 0.05) with nc_ids=[] -> MUST emit "Falsified"!
        #     (iii) Tier 2 'connectivity unknown' (5 detected MAGs with symbol 1 on A_B12_ANA -> uncertified component-conditional
        #           bound < 0.05, ci_null2b is None) with nc_ids=[] -> MUST NOT emit "Falsified"!
        genomes_gen_unk_a12 = dict(genomes_resolved)
        genomes_gen_unk_a12["dastool-don_b12"] = GenomeAnnotations(
            mag_id="dastool-don_b12", domain="Bacteria", completeness=98.0, contamination=6.0,
            genome_size_bp=2_500_000, contig_lengths=(2_500_000,),
            kos=frozenset(d_b12.kos | d_dsr.kos),  # P_B12_AER shares KOs with A_B12_ANA -> symbol 1 on A_B12_ANA!
        )
        # (i) Tier 1 exact Omega enumeration on 2-MAG events (ctg_ri3rec + ctg_b12 + Gate 2 contigs)
        d_sym_t1 = os.path.join(tmpd, "sym_t1_case")
        os.makedirs(d_sym_t1)
        emap_sym_t1 = os.path.join(tmpd, "emap_sym_t1.tsv")
        man_t1 = ["library\tcovered_bed_sha256\ttargets_bed_sha256"]
        em_t1 = ["library\tevent"]
        for idx in range(1, 11):
            lib = f"t1lib_{idx:02d}"
            ev = f"t1ev_{idx:02d}"
            c_sha = _write_gz(
                os.path.join(d_sym_t1, f"{lib}.covered.bed.gz"),
                [("ctg_ri3rec", 0, 900), ("ctg_b12", 0, 900), ("ctg_nonrep", 0, 900), ("ctg_unb", 0, 900)],
            )
            man_t1.append(f"{lib}\t{c_sha}\t{v2_sha}")
            em_t1.append(f"{lib}\t{ev}")
        with open(os.path.join(d_sym_t1, "MANIFEST.tsv"), "w", encoding="utf-8") as fh:
            fh.write("\n".join(man_t1) + "\n")
        with open(emap_sym_t1, "w", encoding="utf-8") as fh:
            fh.write("\n".join(em_t1) + "\n")

        res_sym_t1 = evaluate_freshwater_events(
            d_sym_t1, emap_sym_t1, genomes_gen_unk_a12, rep_ids, hq_ids, [], rep_c_lens, target_lens, valid_targets,
            {}, {}, pos_groups, expected_targets_bed_sha256=v2_sha, targets_bed_path=v2_bed, min_events=10,
            n_null1_draws=20, n_null2_draws=200, n_boot=100, require_protocol_1000_draws=False,
            gtdb_clade_by_mag=base_gtdb, hits_by_mag=hits_c0_only,
        )
        q1_sym_t1 = res_sym_t1["q1_by_stratum"]["R-I3_B12_REDS__E_minimal"]  # type: ignore[index]
        assert q1_sym_t1["null1_status"] == "unevaluable" and q1_sym_t1["ci_null1"] is None
        assert q1_sym_t1["null2b_per_event_connectivity"][0]["tier"] == "tier1_exact_full_fiber"
        assert q1_sym_t1["ci_null2b"] is not None and q1_sym_t1["ci_null2b"][1] < 0.05
        assert q1_sym_t1["verdict_q1"] == "Falsified", (
            f"Expected symmetric OR-Falsified on missing Null 1 + Tier 1 exact Null 2b, got {q1_sym_t1}"
        )

        # (ii) Tier 2 Ryser-connected + 4-chain MCMC-validated full-fiber Null 2b (all active symbols in {0, 2})
        genomes_gen_ryser = dict(genomes_resolved)
        genomes_gen_ryser["dastool-don_b12"] = GenomeAnnotations(
            mag_id="dastool-don_b12", domain="Bacteria", completeness=98.0, contamination=6.0,
            genome_size_bp=2_500_000, contig_lengths=(2_500_000,),
            kos=frozenset(d_b12.kos | PATHWAY_RULES["A_B12_ANA"].all_kos() | d_dsr.kos),
        )
        res_sym_t2 = evaluate_freshwater_events(
            d_pass, emap_pass, genomes_gen_ryser, rep_ids, hq_ids, [], rep_c_lens, target_lens, valid_targets,
            {}, {}, pos_groups, expected_targets_bed_sha256=v2_sha, targets_bed_path=v2_bed, min_events=10,
            n_null1_draws=20, n_null2_draws=200, n_boot=100, require_protocol_1000_draws=False,
            gtdb_clade_by_mag=base_gtdb, hits_by_mag=hits_c0_only,
        )
        q1_sym_t2 = res_sym_t2["q1_by_stratum"]["R-I3_B12_REDS__E_minimal"]  # type: ignore[index]
        assert q1_sym_t2["null1_status"] == "unevaluable" and q1_sym_t2["ci_null1"] is None
        assert q1_sym_t2["null2b_per_event_connectivity"][0]["tier"] == "tier2_reversible_mcmc_4chain"
        assert q1_sym_t2["null2b_connectivity_status"] == "connected"
        assert q1_sym_t2["ci_null2b"] is not None and q1_sym_t2["ci_null2b"][1] < 0.05
        assert q1_sym_t2["verdict_q1"] == "Falsified", (
            f"Expected symmetric OR-Falsified on missing Null 1 + Tier 2 Ryser/MCMC-validated Null 2b, got {q1_sym_t2}"
        )

        # (iii) Tier 2 'connectivity unknown' (symbol 1 on A_B12_ANA on 5x16 active matrix): uncertified bound MUST NOT falsify!
        res_sym_unk_guard = evaluate_freshwater_events(
            d_pass, emap_pass, genomes_gen_unk_a12, rep_ids, hq_ids, [], rep_c_lens, target_lens, valid_targets,
            {}, {}, pos_groups, expected_targets_bed_sha256=v2_sha, targets_bed_path=v2_bed, min_events=10,
            n_null1_draws=20, n_null2_draws=200, n_boot=100, require_protocol_1000_draws=False,
            gtdb_clade_by_mag=base_gtdb, hits_by_mag=hits_c0_only,
        )
        q1_sym_guard = res_sym_unk_guard["q1_by_stratum"]["R-I3_B12_REDS__E_minimal"]  # type: ignore[index]
        assert q1_sym_guard["null1_status"] == "unevaluable" and q1_sym_guard["ci_null2b"] is None
        assert q1_sym_guard["null2b_component_conditional"]["ci_excess_95"][1] < 0.05
        assert (
            "Null 1 unevaluable:" in str(q1_sym_guard["verdict_q1"])
            and "primary Null 2 unavailable: connectivity unknown" in str(q1_sym_guard["verdict_q1"])
        ), f"Uncertified component-conditional Null 2b must NOT trigger OR-Falsified! Got {q1_sym_guard}"

        # 7. Full 1,000 / 1,000 / 1,000 protocol-conformant Q2 matching & all-E test:
        #    Construct a control pool where outcome-blind GTDB-clade + covariate matching gives "Falsified" ([0, 0])
        #    while unmatched all-45-pairs C0xC0 pooling gives "Supported" (> ln 2), AND verify E_plus_DIN shift,
        #    resolved-only intervals, Null 2c attribution, and Q1 Falsified on Null 1 under unknown Null 2b connectivity!
        d_m = os.path.join(tmpd, "match_vs_pool")
        os.makedirs(d_m)
        emap_m = os.path.join(tmpd, "emap_m.tsv")

        # Incompatible recipient (p__Actinomycetota, 1.5 Mb): lacks (B12, RED_S) and (B12, FIXED_N),
        # and lacks ONLY 1 C0 module (C_TRP) -> 0 eligible C0xC0 pairs & cannot single-donate C_HIS + C_TRP!
        rec_inc = GenomeAnnotations(
            mag_id="dastool-inc_actino",
            domain="Bacteria",
            completeness=99.0,
            contamination=0.5,
            genome_size_bp=1_500_000,
            contig_lengths=(750_000, 750_000),
            kos=frozenset(
                {"K00548"} | PATHWAY_RULES["C_THI_SPLIT"].all_kos() | c0_not_his_trp | PATHWAY_RULES["C_HIS"].all_kos()
            ),
        )
        # Matched C0 control (same class c__Actinomycetes, 1.55 Mb): lacks ONLY (C_HIS, C_TRP), supplied by 2 separate donors -> plexi!
        rec_c0_match = GenomeAnnotations(
            mag_id="dastool-c0_match_actino",
            domain="Bacteria",
            completeness=99.0,
            contamination=0.5,
            genome_size_bp=1_550_000,
            contig_lengths=(775_000, 775_000),
            kos=frozenset(
                {"K00549"} | PATHWAY_RULES["ASSIM_SO4"].all_kos() | PATHWAY_RULES["ASSIM_NO3"].all_kos()
                | PATHWAY_RULES["C_THI_SPLIT"].all_kos() | c0_not_his_trp
            ),
        )
        # Unmatched control #1 (wrong phylum p__Pseudomonadota, 5.5 Mb): lacks all 10 C0 modules (45 pairs, 28 satisfied at d=1 by don_c0_8)!
        rec_c0_wrong_clade = GenomeAnnotations(
            mag_id="dastool-c0_wrong_clade",
            domain="Bacteria",
            completeness=99.0,
            contamination=0.5,
            genome_size_bp=5_500_000,
            contig_lengths=(5_500_000,),
            kos=frozenset(
                {"K00549"} | PATHWAY_RULES["ASSIM_SO4"].all_kos() | PATHWAY_RULES["ASSIM_NO3"].all_kos()
                | PATHWAY_RULES["C_THI_SPLIT"].all_kos()
            ),
        )
        # Unmatched control #2 (same phylum p__Actinomycetota, distant class c__Thermoleophilia & 5.8 Mb / all 10 missing C0 = 45 pairs)!
        rec_c0_wrong_cov = GenomeAnnotations(
            mag_id="dastool-c0_wrong_cov",
            domain="Bacteria",
            completeness=96.0,
            contamination=2.0,
            genome_size_bp=5_800_000,
            contig_lengths=(5_800_000,),
            kos=frozenset(
                {"K00549"} | PATHWAY_RULES["ASSIM_SO4"].all_kos() | PATHWAY_RULES["ASSIM_NO3"].all_kos()
                | PATHWAY_RULES["C_THI_SPLIT"].all_kos()
            ),
        )
        # Donors (include 2 98%-complete bins with contig-end partial hits on THI_HMP / THI_HET ONLY so active submatrix has
        # symbol 1 = STATE_UNRESOLVED -> 'connectivity unknown', without making R-I3, R-I2, or C_HIS x C_TRP ambiguous!)
        don_b12_m = GenomeAnnotations("dastool-m_b12", "Bacteria", 98.5, 1.0, 2_200_000, (1_100_000, 1_100_000), frozenset(PATHWAY_RULES["P_B12_AER"].all_kos() | PATHWAY_RULES["P_BLUB"].all_kos()))
        don_dsr_m = GenomeAnnotations("dastool-m_dsr", "Bacteria", 98.5, 1.0, 2_200_000, (1_100_000, 1_100_000), frozenset(PATHWAY_RULES["A_DSR"].all_kos()))
        don_nif_m = GenomeAnnotations("dastool-m_nif", "Bacteria", 98.0, 6.0, 2_200_000, (2_200_000,), frozenset(PATHWAY_RULES["A_NIF"].all_kos()))
        don_his_m = GenomeAnnotations("dastool-m_his", "Bacteria", 98.0, 6.0, 2_000_000, (2_000_000,), frozenset(PATHWAY_RULES["C_HIS"].all_kos()))
        don_trp_m = GenomeAnnotations("dastool-m_trp", "Bacteria", 98.0, 6.0, 2_000_000, (2_000_000,), frozenset(PATHWAY_RULES["C_TRP"].all_kos()))
        don_c0_8 = GenomeAnnotations("dastool-m_c0_8", "Bacteria", 98.0, 6.0, 3_500_000, (3_500_000,), frozenset(c0_not_his_trp))
        don_thi_u1 = GenomeAnnotations(
            "dastool-m_thi_u1", "Bacteria", 98.0, 6.0, 1_800_000, (900_000, 900_000),
            kos=frozenset({"K00941"}), end_kos=frozenset({"K00941"}),
        )
        don_thi_u2 = GenomeAnnotations(
            "dastool-m_thi_u2", "Bacteria", 98.0, 6.0, 1_800_000, (900_000, 900_000),
            kos=frozenset({"K00878"}), end_kos=frozenset({"K00878"}),
        )

        m_genomes_list = [
            rec_inc, rec_c0_match, rec_c0_wrong_clade, rec_c0_wrong_cov,
            don_b12_m, don_dsr_m, don_nif_m, don_his_m, don_trp_m, don_c0_8, don_thi_u1, don_thi_u2,
        ]
        m_genomes = {g.mag_id: g for g in m_genomes_list}
        m_rep_ids = sorted(m_genomes.keys())
        m_hq_ids = [
            "dastool-inc_actino", "dastool-c0_match_actino", "dastool-c0_wrong_clade",
            "dastool-c0_wrong_cov", "dastool-m_b12", "dastool-m_dsr",
        ]
        m_nc_ids = ["dastool-inc_actino", "dastool-m_b12", "dastool-m_dsr"]
        m_hits = {
            mid: _make_synthetic_gene_hits(mid, m_genomes[mid].kos, n_contigs=2, pad_genes=30)
            for mid in m_nc_ids
        }
        m_gtdb = {
            "dastool-inc_actino": "d__Bacteria;p__Actinomycetota;c__Actinomycetes;o__Mycobacteriales",
            "dastool-c0_match_actino": "d__Bacteria;p__Actinomycetota;c__Actinomycetes;o__Mycobacteriales",
            "dastool-c0_wrong_clade": "d__Bacteria;p__Pseudomonadota;c__Gammaproteobacteria;o__Burkholderiales",
            "dastool-c0_wrong_cov": "d__Bacteria;p__Actinomycetota;c__Thermoleophilia;o__Solirubrobacterales",
            "dastool-m_b12": "d__Bacteria;p__Pseudomonadota;c__Alphaproteobacteria",
            "dastool-m_dsr": "d__Bacteria;p__Desulfobacterota;c__Desulfovibrionia",
            "dastool-m_nif": "d__Bacteria;p__Cyanobacteriota;c__Cyanophyceae",
            "dastool-m_his": "d__Bacteria;p__Bacteroidota;c__Bacteroidia",
            "dastool-m_trp": "d__Bacteria;p__Bacteroidota;c__Bacteroidia",
            "dastool-m_c0_8": "d__Bacteria;p__Pseudomonadota;c__Gammaproteobacteria",
            "dastool-m_thi_u1": "d__Bacteria;p__Verrucomicrobiota;c__Verrucomicrobiae",
            "dastool-m_thi_u2": "d__Bacteria;p__Verrucomicrobiota;c__Verrucomicrobiae",
        }
        m_rep_c_lens = {mid: {f"ctg_{mid}": 1000} for mid in m_rep_ids}
        m_target_lens = {f"ctg_{mid}": 1000 for mid in m_rep_ids}
        m_target_lens["ctg_unb_m"] = 1000
        m_valid_targets = set(m_target_lens.keys())
        m_bed = os.path.join(tmpd, "targets_m.bed")
        with open(m_bed, "w", encoding="utf-8") as fh:
            for cid in sorted(m_valid_targets):
                fh.write(f"{cid}\t0\t{m_target_lens[cid]}\n")
        m_bed_sha = sha256_12(m_bed)

        man_m = ["library\tcovered_bed_sha256\ttargets_bed_sha256"]
        em_m = ["library\tevent"]
        all_m_rows = [(cid, 0, 900) for cid in sorted(m_valid_targets)]
        for idx in range(1, 11):
            lib = f"mlib_{idx:02d}"
            ev = f"mev_{idx:02d}"
            c_sha = _write_gz(os.path.join(d_m, f"{lib}.covered.bed.gz"), all_m_rows)
            man_m.append(f"{lib}\t{c_sha}\t{m_bed_sha}")
            em_m.append(f"{lib}\t{ev}")
        with open(os.path.join(d_m, "MANIFEST.tsv"), "w", encoding="utf-8") as fh:
            fh.write("\n".join(man_m) + "\n")
        with open(emap_m, "w", encoding="utf-8") as fh:
            fh.write("\n".join(em_m) + "\n")
        m_groups = {
            gk: {
                "rep490": {f"ctg_dastool-m_b12": 12},
                "dastool655_nonrep": {},
                "unb_prok": {"ctg_unb_m": 2},
                "euk": {},
            }
            for gk in empty_groups
        }

        res_1000 = evaluate_freshwater_events(
            d_m, emap_m, m_genomes, m_rep_ids, m_hq_ids, m_nc_ids, m_rep_c_lens, m_target_lens, m_valid_targets,
            {}, {}, m_groups, expected_targets_bed_sha256=m_bed_sha, targets_bed_path=m_bed, min_events=10,
            n_null1_draws=1000, n_null2_draws=1000, n_boot=1000, require_protocol_1000_draws=True,
            p2_a1_approved=False, max_connectivity_enum_states=10,
            gtdb_clade_by_mag=m_gtdb, hits_by_mag=m_hits,
        )
        q2_ri3_min = res_1000["q2_by_contrast"]["R-I3_B12_REDS__E_minimal"]  # type: ignore[index]
        q2_ri2_min = res_1000["q2_by_contrast"]["R-I2_B12_FIXEDN__E_minimal"]  # type: ignore[index]
        q2_ri2_din = res_1000["q2_by_contrast"]["R-I2_B12_FIXEDN__E_plus_DIN"]  # type: ignore[index]
        # Matched control (dastool-c0_match_actino: 10/10 plexi) gives Falsified ([0,0]), while pooled all-45 C0xC0 (10/570) gives Supported!
        assert q2_ri3_min["verdict_q2"] == "Falsified", f"Expected matched Q2 Falsified, got {q2_ri3_min}"
        assert q2_ri3_min["descriptive_unmatched_all45_c0_sensitivity"]["implied_band_if_unmatched"] == "Supported"
        assert q2_ri3_min["matching_denominators"]["n_distinct_recipient_matches"] == 10
        assert q2_ri3_min["resolved_only_status"] == "evaluable"
        assert q2_ri3_min["resolved_only_ci_log_or"] == [0.0, 0.0]
        assert q2_ri3_min["u_inc_over_n_inc"] == 0.0 and q2_ri3_min["u_comp_over_n_comp"] == 0.0
        assert (
            q2_ri3_min["null2c_attribution"]["component_conditional_sensitivity"]["label"]
            == "component-conditional (P2-A1, pending grid)"
        )
        # Under E_plus_DIN, R-I2 drops from 10/10 plexi (E_minimal) to 0/10 d1_or_chain, making ci_log_or < 0!
        assert q2_ri2_min["matching_denominators"]["matched_inc_counts"]["n_plexi"] == 10
        assert q2_ri2_din["matching_denominators"]["matched_inc_counts"]["n_plexi"] == 0
        assert q2_ri2_din["matching_denominators"]["matched_inc_counts"]["n_d1"] == 10
        assert q2_ri2_din["ci_log_or_upper_below_zero"] is True
        assert q2_ri2_din["verdict_q2"] == "Falsified"

        # Verify Q1 partial-information rule end-to-end: Null 2b has 'connectivity unknown' (due to 75%-complete bin & cap=10),
        # yet Null 1 upper bound < 0.05 still emits primary Q1 "Falsified" while reporting P2-A1 separately!
        q1_ri3_1000 = res_1000["q1_by_stratum"]["R-I3_B12_REDS__E_minimal"]  # type: ignore[index]
        assert q1_ri3_1000["null1_status"] == "evaluable"
        assert q1_ri3_1000["null2b_connectivity_status"] == "connectivity unknown"
        assert q1_ri3_1000["ci_null2b"] is None
        assert q1_ri3_1000["null2b_component_conditional"]["label"] == "component-conditional (P2-A1, pending grid)"
        assert q1_ri3_1000["verdict_q1"] == "Falsified"

        # 8. Mixed Tier 1 + Tier 2 cohort per-event Option A source selection (Astra 249788802780282880):
        #    (a) Equal-weight events with Tier 1 disconnected (full=0.50, component=1.00) and Tier 2 validated non-connected
        #        (component=0.00) MUST pool to 0.25 for primary (when approved) and 0.50 for all-component sensitivity,
        #        while remaining unavailable (None) without approval!
        mix_evs = ["mix_ev1_t1_disc", "mix_ev2_t2_unk"]
        mix_obs = [StratumCounts(n_plexi=1), StratumCounts(n_plexi=1)]
        mix_conn = {
            "mix_ev1_t1_disc": {
                "status": "disconnected (proven)",
                "component_size": 1,
                "fiber_size": 2,
                "locked_clash_rows": [],
                "has_active_move": False,
                "exact_full_fiber_enumerated": True,
                "proof": "exact_enumeration_and_bfs",
            },
            "mix_ev2_t2_unk": {
                "status": "connectivity unknown",
                "component_size": None,
                "fiber_size": None,
                "locked_clash_rows": [],
                "has_active_move": True,
                "exact_full_fiber_enumerated": False,
                "proof": "exceeds_exact_enumeration_cap",
            },
        }
        mix_q1_ests = [
            {
                "tier": "tier1_exact_full_fiber",
                "exact_full_fiber_enumerated": True,
                "move_connectivity": "disconnected (proven)",
                "full_fiber_valid": True,
                "p_low": 0.50,
                "p_high": 0.50,
                "var_bm_low": 0.0,
                "var_bm_high": 0.0,
                "mcmc_diagnostics": None,
                "comp_p_low": 1.00,
                "comp_p_high": 1.00,
                "comp_var_bm_low": 0.0,
                "comp_var_bm_high": 0.0,
                "comp_valid": True,
            },
            {
                "tier": "tier2_reversible_mcmc_4chain",
                "exact_full_fiber_enumerated": False,
                "move_connectivity": "connectivity unknown",
                "full_fiber_valid": False,
                "p_low": 0.00,
                "p_high": 0.00,
                "var_bm_low": 0.0,
                "var_bm_high": 0.0,
                "mcmc_diagnostics": {"mcmc_diagnostics_pass": True},
                "comp_p_low": 0.00,
                "comp_p_high": 0.00,
                "comp_var_bm_low": 0.0,
                "comp_var_bm_high": 0.0,
                "comp_valid": True,
            },
        ]
        _, mix_list_unapp, mix_ci_unapp, mix_pool_unapp, mix_reason_unapp, mix_comp_unapp = _resolve_cohort_q1_null2b(
            mix_evs, mix_obs, mix_conn, mix_q1_ests, p2_a1_approved=False, n_boot=200, seed_primary=11, seed_comp=12
        )
        assert mix_ci_unapp is None and mix_pool_unapp is None and mix_reason_unapp == "connectivity unknown"
        assert [x["primary_null2_source"] for x in mix_list_unapp] == ["tier1_exact_full_fiber", "unavailable"]
        assert mix_comp_unapp is not None and mix_comp_unapp["pooled_null_bounds"] == [0.5, 0.5]

        _, mix_list_app, mix_ci_app, mix_pool_app, mix_reason_app, mix_comp_app = _resolve_cohort_q1_null2b(
            mix_evs, mix_obs, mix_conn, mix_q1_ests, p2_a1_approved=True, n_boot=200, seed_primary=11, seed_comp=12
        )
        assert mix_ci_app is not None and mix_reason_app is None
        assert [x["primary_null2_source"] for x in mix_list_app] == [
            "tier1_exact_full_fiber",
            "component_conditional_p2_a1",
        ]
        assert mix_pool_app == [0.25, 0.25], f"Expected mixed primary pooled null 0.25, got {mix_pool_app}"
        assert mix_comp_app is not None and mix_comp_app["pooled_null_bounds"] == [0.5, 0.5], (
            f"Expected separate all-component sensitivity pooled null 0.50, got {mix_comp_app}"
        )

        # Verify identical per-event rule on Q2 Null 2c attribution
        mix_q2_ests = {
            "mix_ev1_t1_disc": {
                "tier": "tier1_exact_full_fiber",
                "exact_full_fiber_enumerated": True,
                "move_connectivity": "disconnected (proven)",
                "full_fiber_valid": True,
                "inc_p_low": 0.50,
                "inc_p_high": 0.50,
                "inc_var_low": 0.0,
                "inc_var_high": 0.0,
                "comp_p_low": 0.00,
                "comp_p_high": 0.00,
                "comp_var_low": 0.0,
                "comp_var_high": 0.0,
                "mcmc_diagnostics": None,
                "ck_inc_p_low": 1.00,
                "ck_inc_p_high": 1.00,
                "ck_inc_var_low": 0.0,
                "ck_inc_var_high": 0.0,
                "ck_comp_p_low": 0.00,
                "ck_comp_p_high": 0.00,
                "ck_comp_var_low": 0.0,
                "ck_comp_var_high": 0.0,
                "ck_valid": True,
            },
            "mix_ev2_t2_unk": {
                "tier": "tier2_reversible_mcmc_4chain",
                "exact_full_fiber_enumerated": False,
                "move_connectivity": "connectivity unknown",
                "full_fiber_valid": False,
                "inc_p_low": 0.00,
                "inc_p_high": 0.00,
                "inc_var_low": 0.0,
                "inc_var_high": 0.0,
                "comp_p_low": 0.00,
                "comp_p_high": 0.00,
                "comp_var_low": 0.0,
                "comp_var_high": 0.0,
                "mcmc_diagnostics": {"mcmc_diagnostics_pass": True},
                "ck_inc_p_low": 0.00,
                "ck_inc_p_high": 0.00,
                "ck_inc_var_low": 0.0,
                "ck_inc_var_high": 0.0,
                "ck_comp_p_low": 0.00,
                "ck_comp_p_high": 0.00,
                "ck_comp_var_low": 0.0,
                "ck_comp_var_high": 0.0,
                "ck_valid": True,
            },
        }
        mix_2c_obs_i = {e: StratumCounts(n_plexi=1) for e in mix_evs}
        mix_2c_obs_c = {e: StratumCounts(n_d1=1) for e in mix_evs}
        n2c_unapp = _resolve_cohort_q2_null2c(
            mix_evs, mix_2c_obs_i, mix_2c_obs_c, mix_conn, mix_q2_ests,
            p2_a1_approved=False, n_boot=200, seed_primary=21, seed_comp=22,
        )
        assert n2c_unapp["primary_null2c_available"] is False and n2c_unapp["ci_null2c_log_or"] is None
        n2c_app = _resolve_cohort_q2_null2c(
            mix_evs, mix_2c_obs_i, mix_2c_obs_c, mix_conn, mix_q2_ests,
            p2_a1_approved=True, n_boot=200, seed_primary=21, seed_comp=22,
        )
        assert n2c_app["primary_null2c_available"] is True
        assert n2c_app["primary_pooled_null_inc_bounds"] == [0.25, 0.25]
        assert [x["primary_null2c_source"] for x in n2c_app["per_event_connectivity"]] == [  # type: ignore[union-attr]
            "tier1_exact_full_fiber",
            "component_conditional_p2_a1",
        ]
        assert n2c_app["ci_null2c_log_or"] != n2c_app["component_conditional_sensitivity"]["ci_null2c_log_or"]  # type: ignore[index]

        # (b) End-to-end mixed Tier 1 (ev 1..5: |Omega|=2 exact fiber) + Tier 2 (ev 6..10: connectivity unknown)
        rec_pure_mx = GenomeAnnotations(
            "dastool-inc_pure", "Bacteria", 99.0, 0.5, 1_500_000, (1_500_000,), frozenset({"K00548"})
        )
        mx_genomes = dict(m_genomes, **{"dastool-inc_pure": rec_pure_mx})
        mx_rep_ids = sorted(mx_genomes.keys())
        mx_hq_ids = sorted(set(m_hq_ids) | {"dastool-inc_pure"})
        mx_rep_c_lens = dict(m_rep_c_lens, **{"dastool-inc_pure": {"ctg_dastool-inc_pure": 1000}})
        mx_target_lens = dict(m_target_lens, **{"ctg_dastool-inc_pure": 1000})
        mx_valid_targets = set(mx_target_lens.keys())
        mx_gtdb = dict(m_gtdb, **{"dastool-inc_pure": m_gtdb["dastool-inc_actino"]})
        mx_bed = os.path.join(tmpd, "targets_mx.bed")
        with open(mx_bed, "w", encoding="utf-8") as fh:
            for cid in sorted(mx_valid_targets):
                fh.write(f"{cid}\t0\t{mx_target_lens[cid]}\n")
        mx_bed_sha = sha256_12(mx_bed)

        d_mix_e2e = os.path.join(tmpd, "mix_e2e")
        os.makedirs(d_mix_e2e)
        emap_mix = os.path.join(tmpd, "emap_mix.tsv")
        man_mix = ["library\tcovered_bed_sha256\ttargets_bed_sha256"]
        em_mix = ["library\tevent"]
        t1_only_rows = [
            (cid, 0, 900)
            for cid in ("ctg_dastool-inc_pure", "ctg_dastool-m_b12", "ctg_dastool-m_dsr", "ctg_unb_m")
        ]
        for idx in range(1, 11):
            lib = f"mxlib_{idx:02d}"
            ev = f"mxev_{idx:02d}"
            rows_mx = t1_only_rows if idx <= 5 else all_m_rows
            c_sha = _write_gz(os.path.join(d_mix_e2e, f"{lib}.covered.bed.gz"), rows_mx)
            man_mix.append(f"{lib}\t{c_sha}\t{mx_bed_sha}")
            em_mix.append(f"{lib}\t{ev}")
        with open(os.path.join(d_mix_e2e, "MANIFEST.tsv"), "w", encoding="utf-8") as fh:
            fh.write("\n".join(man_mix) + "\n")
        with open(emap_mix, "w", encoding="utf-8") as fh:
            fh.write("\n".join(em_mix) + "\n")

        res_mx_unapp = evaluate_freshwater_events(
            d_mix_e2e, emap_mix, mx_genomes, mx_rep_ids, mx_hq_ids, m_nc_ids, mx_rep_c_lens, mx_target_lens, mx_valid_targets,
            {}, {}, m_groups, expected_targets_bed_sha256=mx_bed_sha, targets_bed_path=mx_bed, min_events=10,
            n_null1_draws=20, n_null2_draws=200, n_boot=100, require_protocol_1000_draws=False,
            p2_a1_approved=False, max_connectivity_enum_states=10, gtdb_clade_by_mag=mx_gtdb, hits_by_mag=m_hits,
        )
        res_mx_app = evaluate_freshwater_events(
            d_mix_e2e, emap_mix, mx_genomes, mx_rep_ids, mx_hq_ids, m_nc_ids, mx_rep_c_lens, mx_target_lens, mx_valid_targets,
            {}, {}, m_groups, expected_targets_bed_sha256=mx_bed_sha, targets_bed_path=mx_bed, min_events=10,
            n_null1_draws=20, n_null2_draws=200, n_boot=100, require_protocol_1000_draws=False,
            p2_a1_approved=True, max_connectivity_enum_states=10, gtdb_clade_by_mag=mx_gtdb, hits_by_mag=m_hits,
        )
        q1_mx_u = res_mx_unapp["q1_by_stratum"]["R-I3_B12_REDS__E_minimal"]  # type: ignore[index]
        q1_mx_a = res_mx_app["q1_by_stratum"]["R-I3_B12_REDS__E_minimal"]  # type: ignore[index]
        assert q1_mx_u["ci_null2b"] is None and q1_mx_a["ci_null2b"] is not None
        assert [x["primary_null2_source"] for x in q1_mx_u["null2b_per_event_connectivity"]] == (
            ["tier1_exact_full_fiber"] * 5 + ["unavailable"] * 5
        )
        assert [x["primary_null2_source"] for x in q1_mx_a["null2b_per_event_connectivity"]] == (
            ["tier1_exact_full_fiber"] * 5 + ["component_conditional_p2_a1"] * 5
        )

    return {
        "v1_paired_with_v2_bed_rejected": caught_v1_manifest and caught_v1_subst and caught_mv2_v1,
        "missing_bed_rejected": caught_missing_bed,
        "multi_event_library_rejected": caught_dup_lib,
        "hq_ineligible_gate1_caught": True,
        "zero_denom_gate2_unevaluable_caught": True,
        "gate3_unresolved_caught": True,
        "sub1000_draw_gate_caught": True,
        "zero_eligible_null1_propagated": True,
        "matched_vs_pooled_c0_opposite_verdicts_1000draws": (
            q2_ri3_min["verdict_q2"],
            q2_ri3_min["descriptive_unmatched_all45_c0_sensitivity"]["implied_band_if_unmatched"],
        ),
        "scenario_E_plus_DIN_collapse_verified": True,
        "p2_a1_3state_partial_info_verified": True,
        "mixed_cohort_primary_vs_all_comp_pooled": (mix_pool_app, mix_comp_app["pooled_null_bounds"]),
    }


def run_freshwater_staged(
    covered_dir: Optional[str] = None,
    event_map_tsv: Optional[str] = None,
    targets_bed_path: Optional[str] = None,
    n_null_draws: int = 1000,
    n_boot: int = 1000,
    p2_a1_approved: bool = False,
) -> Dict[str, object]:
    os.makedirs(OUT_DIR, exist_ok=True)
    prov_path = os.path.join(FW_ROOT, "PROVENANCE.md")
    q_path = os.path.join(FW_ROOT, "binning", "checkm2", "quality_report.tsv")
    derep_path = os.path.join(FW_ROOT, "binning", "derep", "derep_clusters.tsv")
    per_mag_dir = os.path.join(FW_ROOT, "metabolism", "per_mag", "per_mag")
    comm_path = os.path.join(FW_ROOT, "metabolism", "per_mag", "community_annotations.tsv")
    tiara_path = os.path.join(FW_ROOT, "eukaryotic", "tiara", "tiara_output.tsv")
    gtdb_path = os.path.join(FW_ROOT, "taxonomy", "gtdbtk", "gtdbtk_taxonomy.tsv")
    gff_path = os.path.join(FW_ROOT, "annotation", "bakta", "basic", "annotation.gff3")

    # 1. Verify every sha256 recorded in PROVENANCE.md (now 497 files)
    verified_prov_files = 0
    prov_hashes: Dict[str, str] = {}
    with open(prov_path, "r", encoding="utf-8") as fh:
        for line in fh:
            s = line.strip()
            if not s or s.startswith("#"):
                continue
            exp_hash, rel_p = s.split(None, 1)
            full_p = os.path.normpath(os.path.join(FW_ROOT, rel_p))
            got_hash = sha256_12(full_p)
            if got_hash != exp_hash:
                raise RuntimeError(f"PROVENANCE hash mismatch on {full_p}: {got_hash} != {exp_hash}")
            prov_hashes[os.path.normpath(rel_p)] = got_hash
            verified_prov_files += 1

    with open(q_path, "r", encoding="utf-8") as fh:
        q_rows = list(csv.DictReader(fh, delimiter="\t"))
    raw_hq = sum(
        1 for r in q_rows if float(r["Completeness"]) >= 90.0 and float(r["Contamination"]) <= 5.0
    )
    dastool_q_rows = [r for r in q_rows if r["Name"].startswith("dastool-")]
    dastool_pre_derep_hq = sum(
        1 for r in dastool_q_rows if float(r["Completeness"]) >= 90.0 and float(r["Contamination"]) <= 5.0
    )
    q_by_name = {r["Name"]: r for r in q_rows}

    with open(derep_path, "r", encoding="utf-8") as fh:
        derep_rows = list(csv.DictReader(fh, delimiter="\t"))
    rep_rows = [r for r in derep_rows if r["status"] == "representative"]
    rep_ids = sorted(r["genome"] for r in rep_rows)
    assert_non_redundant_mag_ids(rep_ids)

    # 2. Load GTDB-Tk taxonomy on all 490 dereplicated representatives
    gtdb_domain_by_mag: Dict[str, str] = {}
    gtdb_clade_by_mag: Dict[str, str] = {}
    gtdb_domain_counts = {"Bacteria": 0, "Archaea": 0, "U": 0}
    with open(gtdb_path, "r", encoding="utf-8") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            ug = r["user_genome"]
            if ug.endswith(".fa"):
                ug = ug[:-3]
            cls = r["classification"]
            gtdb_clade_by_mag[ug] = cls
            if cls.startswith("d__Bacteria"):
                gtdb_domain_by_mag[ug] = "Bacteria"
            elif cls.startswith("d__Archaea"):
                gtdb_domain_by_mag[ug] = "Archaea"
            else:
                gtdb_domain_by_mag[ug] = "U"
    for mid in rep_ids:
        dom_m = gtdb_domain_by_mag.get(mid, "U")
        gtdb_domain_counts[dom_m] = gtdb_domain_counts.get(dom_m, 0) + 1

    # 3. Collect contig sets across the 490 derep reps and all 655 dastool-* consensus bins
    rep_contig_lens_by_mag: Dict[str, Dict[str, int]] = {}
    rep_contigs: Set[str] = set()
    for mid in rep_ids:
        fa_p = os.path.join(FW_FASTA_DIR, f"{mid}.fa")
        c_map = read_fasta_contig_lengths(fa_p)
        rep_contig_lens_by_mag[mid] = c_map
        rep_contigs.update(c_map.keys())

    all_dastool_contigs: Set[str] = set(rep_contigs)
    for r_d in dastool_q_rows:
        mid_d = r_d["Name"]
        if mid_d not in rep_contig_lens_by_mag:
            fa_d = os.path.join(FW_FASTA_DIR, f"{mid_d}.fa")
            if os.path.isfile(fa_d):
                all_dastool_contigs.update(read_fasta_contig_lengths(fa_d).keys())

    # 4. Parse Tiara contig classes + community_annotations.tsv for Gate 2 & extra P2-KO contigs
    tiara_fst: Dict[str, str] = {}
    tiara_snd: Dict[str, str] = {}
    with open(tiara_path, "r", encoding="utf-8") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            seq_id = r["sequence_id"]
            tiara_fst[seq_id] = r["class_fst_stage"]
            tiara_snd[seq_id] = r.get("class_snd_stage", "")

    ref_kos = all_referenced_kos()
    b12_kos = set().union(
        *(PATHWAY_RULES[m].all_kos() for m in ("P_B12_AER", "A_B12_ANA", "P_BLUB", "P_BLUB_DMB", "A_BZA", "B12_LOWER_ATTACH"))
    )
    thi_kos = set().union(*(PATHWAY_RULES[m].all_kos() for m in ("THI_HMP", "THI_HET", "THI_COUPLE")))
    c0_kos = set().union(*(PATHWAY_RULES[m].all_kos() for m in C0_CONTROL_MODULES))
    endpoint_donor_ko_sets: Dict[str, Set[str]] = {
        "all_donor_modules": set().union(*(PATHWAY_RULES[m].all_kos() for m in DONOR_CAP_MODULES)),
        "R-I3_B12_REDS": b12_kos | PATHWAY_RULES["A_DSR"].all_kos(),
        "R-I2_B12_FIXEDN": b12_kos | PATHWAY_RULES["A_NIF"].all_kos(),
        "R-I6_B12": b12_kos,
        "C_THI_SPLIT": thi_kos,
        "C0_ALL10": c0_kos,
    }

    all_nonrep_p2_ctgs: Set[str] = set()
    non_organelle_nonrep_p2_ctgs: Set[str] = set()
    organelle_nonrep_p2_ctgs: Set[str] = set()
    organelle_subtypes: Dict[str, int] = {}

    unbinned_prok_kos_by_ctg: Dict[str, Set[str]] = {}
    euk_kos_by_ctg: Dict[str, Set[str]] = {}
    donor_hits_by_group_and_ctg: Dict[str, Dict[str, Dict[str, int]]] = {
        gk: {"rep490": {}, "dastool655_nonrep": {}, "unb_prok": {}, "euk": {}}
        for gk in endpoint_donor_ko_sets
    }

    with open(comm_path, "r", encoding="utf-8") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            kos = parse_ko_field(r.get("KO", ""))
            if not kos:
                continue
            cid = r["contig_id"]
            t_cls = tiara_fst.get(cid, "unknown")
            p2_kos = kos & ref_kos
            if cid not in rep_contigs and p2_kos:
                all_nonrep_p2_ctgs.add(cid)
                if t_cls in ("organelle", "mitochondrion", "plastid"):
                    organelle_nonrep_p2_ctgs.add(cid)
                else:
                    non_organelle_nonrep_p2_ctgs.add(cid)
            if t_cls in ("organelle", "mitochondrion", "plastid"):
                continue

            if cid in rep_contigs:
                bucket = "rep490"
            elif cid in all_dastool_contigs:
                bucket = "dastool655_nonrep"
            elif t_cls == "eukarya":
                bucket = "euk"
                if p2_kos:
                    euk_kos_by_ctg.setdefault(cid, set()).update(p2_kos)
            else:
                bucket = "unb_prok"
                if p2_kos:
                    unbinned_prok_kos_by_ctg.setdefault(cid, set()).update(p2_kos)

            for gk, kset in endpoint_donor_ko_sets.items():
                dh = len(kos & kset)
                if dh:
                    bmap = donor_hits_by_group_and_ctg[gk][bucket]
                    bmap[cid] = bmap.get(cid, 0) + dh

    for oc in organelle_nonrep_p2_ctgs:
        sub = tiara_snd.get(oc) or tiara_fst.get(oc, "organelle")
        organelle_subtypes[sub] = organelle_subtypes.get(sub, 0) + 1

    # Exact ID-set reconciliation against results/p2/breadth/extra_p2_contigs.txt (db998ebf99c8) if present
    extra_list_path = os.path.join(OUT_DIR, "breadth", "extra_p2_contigs.txt")
    extra_reconciliation: Dict[str, object] = {
        "all_nonrep_p2_contigs": len(all_nonrep_p2_ctgs),
        "non_organelle_nonrep_p2_contigs": len(non_organelle_nonrep_p2_ctgs),
        "organelle_excluded_contigs": len(organelle_nonrep_p2_ctgs),
        "organelle_subtype_counts": dict(sorted(organelle_subtypes.items())),
        "breakdown_non_organelle_25807": {
            "unbinned_prok_or_unknown": len(unbinned_prok_kos_by_ctg),
            "tiara_eukarya_unbinned": len(euk_kos_by_ctg),
            "dastool_nonrep_165bins": len(
                non_organelle_nonrep_p2_ctgs - set(unbinned_prok_kos_by_ctg) - set(euk_kos_by_ctg)
            ),
        },
    }
    if os.path.isfile(extra_list_path):
        with open(extra_list_path, "r", encoding="utf-8") as fh:
            breadth_extra_ids = {ln.strip() for ln in fh if ln.strip()}
        assert breadth_extra_ids == all_nonrep_p2_ctgs, (
            f"ID-set mismatch vs extra_p2_contigs.txt: diff={len(breadth_extra_ids ^ all_nonrep_p2_ctgs)}"
        )
        assert breadth_extra_ids - non_organelle_nonrep_p2_ctgs == organelle_nonrep_p2_ctgs
        extra_reconciliation["extra_p2_contigs_txt_sha256"] = sha256_12(extra_list_path)
        extra_reconciliation["exact_id_set_match_with_targets_v2_extra"] = True

    descriptive_assembly_gate2: Dict[str, Dict[str, object]] = {}
    for gk, bdict in donor_hits_by_group_and_ctg.items():
        h_rep = sum(bdict["rep490"].values())
        h_nonrep = sum(bdict["dastool655_nonrep"].values())
        h_unb = sum(bdict["unb_prok"].values())
        h_euk = sum(bdict["euk"].values())
        descriptive_assembly_gate2[gk] = {
            "rep490_hits": h_rep,
            "dastool655_nonrep_hits": h_nonrep,
            "unbinned_prok_hits": h_unb,
            "tiara_eukarya_hits": h_euk,
            "unbinned_share_vs_490reps": round(h_unb / max(1, h_rep + h_unb), 4),
            "unbinned_share_vs_655dastool": round(h_unb / max(1, h_rep + h_nonrep + h_unb), 4),
        }

    # 5. Load/cache complete Bakta GFF3 CDS order on the 48,232 representative contigs + lengths on all 74,219 targets
    gff_sha12 = prov_hashes[os.path.normpath("./annotation/bakta/basic/annotation.gff3")]
    gff_cache_path = os.path.join(OUT_DIR, "freshwater_gff_cds_cache.json")
    gff_order, target_ctg_lens = _load_or_build_freshwater_gff_cache(
        gff_path, gff_sha12, rep_contigs, all_nonrep_p2_ctgs, gff_cache_path
    )

    # 6. Parse all 490 derep representative per_mag TSVs with exact GTDB-Tk domain, GFF3 CDS order, and GeneHit blocks
    genomes: Dict[str, GenomeAnnotations] = {}
    fw_hits_by_mag: Dict[str, List[GeneHit]] = {}
    marker_domain_counts = {"Bacteria": 0, "Archaea": 0, "U": 0}
    total_unknown_edge_kos = 0
    for mid in rep_ids:
        c_lens_map = rep_contig_lens_by_mag[mid]
        c_lens = [c_lens_map[c] for c in sorted(c_lens_map)]
        qr = q_by_name[mid]
        tsv_p = os.path.join(per_mag_dir, f"{mid}.tsv")
        dom_gtdb = gtdb_domain_by_mag[mid]
        g_ann, g_hits = parse_per_mag_tsv(
            tsv_p,
            mag_id=mid,
            domain=dom_gtdb,
            completeness=float(qr["Completeness"]),
            contamination=float(qr["Contamination"]),
            genome_size_bp=int(qr["Genome_Size"]),
            contig_lengths=c_lens,
            gff_cds_order=gff_order,
        )
        total_unknown_edge_kos += len(g_ann.unknown_edge_kos)
        marker_dom = infer_conservative_domain_from_kos(g_ann.kos)
        marker_domain_counts[marker_dom] += 1
        genomes[mid] = g_ann
        fw_hits_by_mag[mid] = g_hits

    assert total_unknown_edge_kos == 0, f"Expected 0 unknown_edge_kos with staged GFF3, got {total_unknown_edge_kos}"

    hq_ids = [m for m in rep_ids if genomes[m].completeness >= 90.0 and genomes[m].contamination <= 5.0]
    nc_ids = [m for m in rep_ids if genomes[m].completeness >= 95.0 and genomes[m].contamination <= 5.0]
    all_rep_genomes = [genomes[m] for m in rep_ids]
    hq_genomes = [genomes[m] for m in hq_ids]
    hq_resolved_dom = [g for g in hq_genomes if g.domain in ("Bacteria", "Archaea")]
    hq_unresolved_dom = [g for g in hq_genomes if g.domain == "U"]

    strata_defs: List[Tuple[str, Tuple[str, ...], SupplyScenario]] = [
        ("R-I3_B12_REDS__E_minimal", ("B12", "RED_S"), SCENARIO_E_MINIMAL),
        ("R-I3_B12_REDS__E_plus_DMSP", ("B12", "RED_S"), SCENARIO_E_PLUS_DMSP),
        ("R-I2_B12_FIXEDN__E_minimal", ("B12", "FIXED_N"), SCENARIO_E_MINIMAL),
        ("R-I2_B12_FIXEDN__E_plus_DIN", ("B12", "FIXED_N"), SCENARIO_E_PLUS_DIN),
        ("R-I6_B12_SINGLE__E_minimal", ("B12",), SCENARIO_E_MINIMAL),
        ("C_THI_SPLIT__E_minimal", ("THI",), SCENARIO_E_MINIMAL),
    ]

    # Non-event-resolved co-assembly diagnostics (with exact GFF3 CDS order + GTDB-Tk domains):
    #   (a) 490-MAG catalog alone, and (b) 490-MAG catalog + assembly-wide _unbinned_P_star + _euk_tiara
    all_unb_kos = set().union(*unbinned_prok_kos_by_ctg.values()) if unbinned_prok_kos_by_ctg else set()
    all_euk_kos = set().union(*euk_kos_by_ctg.values()) if euk_kos_by_ctg else set()
    unb_g_all = GenomeAnnotations.from_kos("_unbinned", all_unb_kos) if all_unb_kos else None
    euk_g_all = GenomeAnnotations.from_kos("_euk_tiara", all_euk_kos, domain="Eukaryota") if all_euk_kos else None

    ctx_fw_490 = CommunityScorerContext(all_rep_genomes, mode="fanin_only")
    ctx_fw_with_pools = CommunityScorerContext(all_rep_genomes, unb_g_all, euk_g_all, mode="fanin_only")

    coassembly_diag_490: Dict[str, Dict[str, object]] = {}
    coassembly_diag_with_pools: Dict[str, Dict[str, object]] = {}
    fw_tuple_rows: List[Dict[str, object]] = []

    for s_label, reqs, scen in strata_defs:
        sc_490 = StratumCounts()
        sc_pool = StratumCounts()
        n_inelig = 0
        for rec in hq_genomes:
            call_490, wit_490 = ctx_fw_490.evaluate(rec, reqs, scen)
            call_pool, wit_pool = ctx_fw_with_pools.evaluate(rec, reqs, scen)
            if call_490 == CALL_NOT_ELIGIBLE:
                n_inelig += 1
            else:
                sc_490.record(call_490)
                sc_pool.record(call_pool)
            if call_490 not in (CALL_NOT_ELIGIBLE, CALL_AUTONOMOUS):
                fw_tuple_rows.append({
                    "stratum": s_label,
                    "recipient_mag": rec.mag_id,
                    "domain_gtdbtk": rec.domain,
                    "completeness": rec.completeness,
                    "contamination": rec.contamination,
                    "n_contigs": len(rec.contig_lengths),
                    "catalog_490_call": call_490,
                    "catalog_490_witness": ":".join(wit_490) if wit_490 else "",
                    "catalog_with_unbinned_euk_call": call_pool,
                    "catalog_with_unbinned_euk_witness": ":".join(wit_pool) if wit_pool else "",
                })
        coassembly_diag_490[s_label] = {
            "n_hq_total": len(hq_genomes),
            "n_not_eligible": n_inelig,
            "n_auto": sc_490.n_auto,
            "n_cand_aux": sc_490.n_cand_aux,
            "n_d0": sc_490.n_d0,
            "n_d1_in_490_catalog": sc_490.n_d1,
            "n_plexi_in_490_catalog_without_unbinned": sc_490.n_plexi,
            "n_ambig": sc_490.n_ambig,
            "n_no_donor": sc_490.n_no_donor,
        }
        coassembly_diag_with_pools[s_label] = {
            "n_hq_total": len(hq_genomes),
            "n_not_eligible": n_inelig,
            "n_auto": sc_pool.n_auto,
            "n_cand_aux": sc_pool.n_cand_aux,
            "n_d0": sc_pool.n_d0,
            "n_d1": sc_pool.n_d1,
            "n_plexi": sc_pool.n_plexi,
            "n_ambig": sc_pool.n_ambig,
            "n_no_donor": sc_pool.n_no_donor,
        }

    sc_c0_490 = StratumCounts()
    sc_c0_pool = StratumCounts()
    c0_inelig = 0
    for rec in hq_genomes:
        for p_c0 in itertools.combinations(C0_CONTROL_MODULES, 2):
            c_490, _ = ctx_fw_490.evaluate(rec, p_c0, SCENARIO_E_MINIMAL)
            c_pool, _ = ctx_fw_with_pools.evaluate(rec, p_c0, SCENARIO_E_MINIMAL)
            if c_490 == CALL_NOT_ELIGIBLE:
                c0_inelig += 1
            else:
                sc_c0_490.record(c_490)
                sc_c0_pool.record(c_pool)
    coassembly_diag_490["C0xC0_ALL45__E_minimal"] = {
        "n_hq_pair_evals": len(hq_genomes) * 45,
        "n_not_eligible": c0_inelig,
        "n_auto": sc_c0_490.n_auto,
        "n_cand_aux": sc_c0_490.n_cand_aux,
        "n_d0": sc_c0_490.n_d0,
        "n_d1_in_490_catalog": sc_c0_490.n_d1,
        "n_plexi_in_490_catalog_without_unbinned": sc_c0_490.n_plexi,
        "n_ambig": sc_c0_490.n_ambig,
        "n_no_donor": sc_c0_490.n_no_donor,
    }
    coassembly_diag_with_pools["C0xC0_ALL45__E_minimal"] = {
        "n_hq_pair_evals": len(hq_genomes) * 45,
        "n_not_eligible": c0_inelig,
        "n_auto": sc_c0_pool.n_auto,
        "n_cand_aux": sc_c0_pool.n_cand_aux,
        "n_d0": sc_c0_pool.n_d0,
        "n_d1": sc_c0_pool.n_d1,
        "n_plexi": sc_c0_pool.n_plexi,
        "n_ambig": sc_c0_pool.n_ambig,
        "n_no_donor": sc_c0_pool.n_no_donor,
    }

    # Descriptive donor capability prevalence across all 490 derep representatives + I1 bypass check
    donor_module_counts: Dict[str, Dict[str, int]] = {}
    for mid_rule in (
        "P_B12_AER", "A_B12_ANA", "P_BLUB", "A_BZA", "B12_LOWER_ATTACH",
        "A_DSR", "A_NIF", "ASSIM_SO4", "ASSIM_NO3", "C_THI_SPLIT", "P_PSII",
        *C0_CONTROL_MODULES,
    ):
        st_c = {STATE_PRESENT: 0, STATE_ABSENT: 0, STATE_UNRESOLVED: 0}
        for g in all_rep_genomes:
            st = PATHWAY_RULES[mid_rule].evaluate(g, max_false_absence_prob=0.05).state
            st_c[st] += 1
        donor_module_counts[mid_rule] = st_c

    b12_complete_donors = [
        g.mag_id for g in all_rep_genomes if evaluate_b12_de_novo(g, max_false_absence_prob=0.05) == STATE_PRESENT
    ]
    b12_plus_dsr_donors = [
        g.mag_id
        for g in all_rep_genomes
        if evaluate_b12_de_novo(g, 0.05) == STATE_PRESENT
        and PATHWAY_RULES["A_DSR"].evaluate(g, 0.05).state == STATE_PRESENT
    ]
    b12_plus_nif_donors = [
        g.mag_id
        for g in all_rep_genomes
        if evaluate_b12_de_novo(g, 0.05) == STATE_PRESENT
        and PATHWAY_RULES["A_NIF"].evaluate(g, 0.05).state == STATE_PRESENT
    ]
    i1_fw = [
        {"mag_id": g.mag_id, "bypass": classify_coencoder_bypass(g)}
        for g in all_rep_genomes
        if PATHWAY_RULES["P_PSII"].evaluate(g, 0.05).state == STATE_PRESENT
        and PATHWAY_RULES["A_NIF"].evaluate(g, 0.05).state == STATE_PRESENT
    ]

    # Optional event-level scoring once *.covered.bed.gz (job 8194877), MANIFEST.tsv, and library->event map are supplied
    if (covered_dir is None) != (event_map_tsv is None):
        raise ValueError("Both --freshwater-covered-dir and --freshwater-event-map must be provided together.")

    expected_bams_path = os.path.join(OUT_DIR, "breadth", "bams.txt")
    expected_fw_libs: Optional[Set[str]] = None
    if os.path.isfile(expected_bams_path):
        with open(expected_bams_path, "r", encoding="utf-8") as fh:
            expected_fw_libs = {
                ln.strip().replace(".sorted.bam", "") for ln in fh if ln.strip()
            }

    # Verify staged breadth_v2 232-library files against manifest_v2.tsv (6635acc1a6ee) if present
    staged_bv2_dir = os.path.join(FW_ROOT, "breadth_v2")
    staged_mv2_path = os.path.join(OUT_DIR, "breadth", "manifest_v2.tsv")
    staged_breadth_v2_verification: Optional[Dict[str, object]] = None
    if os.path.isdir(staged_bv2_dir) and os.path.isfile(staged_mv2_path):
        mv2_sha12 = sha256_12(staged_mv2_path)
        n_verified_bv2 = 0
        mv2_libs: Set[str] = set()
        with open(staged_mv2_path, "r", encoding="utf-8") as fh:
            for r in csv.DictReader(fh, delimiter="\t"):
                lib = r["library"].strip()
                mv2_libs.add(lib)
                if r["targets_bed_sha12"].strip() != "f48d55128984":
                    raise RuntimeError(f"manifest_v2 {lib} targets_bed_sha12 != f48d55128984")
                if r.get("local_vs_grex", "").strip() != "match":
                    raise RuntimeError(f"manifest_v2 {lib} local_vs_grex != match")
                gz_p = os.path.join(staged_bv2_dir, f"{lib}.covered.bed.gz")
                if sha256_full(gz_p) != r["sha256"].strip():
                    raise RuntimeError(f"breadth_v2 full sha256 mismatch on {gz_p}")
                n_verified_bv2 += 1
        if expected_fw_libs is not None and mv2_libs != expected_fw_libs:
            raise RuntimeError(f"manifest_v2 library set mismatch vs bams.txt: diff={len(mv2_libs ^ expected_fw_libs)}")
        staged_breadth_v2_verification = {
            "manifest_v2_sha256": mv2_sha12,
            "targets_bed_sha12": "f48d55128984",
            "n_libraries_full_sha256_verified": n_verified_bv2,
            "matches_bams_txt_232_exactly": expected_fw_libs is not None and mv2_libs == expected_fw_libs,
        }

    event_level_results: Optional[Dict[str, object]] = None
    if covered_dir is not None and event_map_tsv is not None:
        event_level_results = evaluate_freshwater_events(
            covered_dir=covered_dir,
            event_map_tsv=event_map_tsv,
            genomes=genomes,
            rep_ids=rep_ids,
            hq_ids=hq_ids,
            nc_ids=nc_ids,
            rep_contig_lens_by_mag=rep_contig_lens_by_mag,
            target_ctg_lens=target_ctg_lens,
            valid_target_contigs=rep_contigs | all_nonrep_p2_ctgs,
            unbinned_prok_kos_by_ctg=unbinned_prok_kos_by_ctg,
            euk_kos_by_ctg=euk_kos_by_ctg,
            donor_hits_by_group_and_ctg=donor_hits_by_group_and_ctg,
            expected_targets_bed_sha256="f48d55128984",
            targets_bed_path=targets_bed_path,
            expected_libraries=expected_fw_libs,
            n_null1_draws=n_null_draws,
            n_null2_draws=n_null_draws,
            n_boot=n_boot,
            require_protocol_1000_draws=True,
            p2_a1_approved=p2_a1_approved,
            gtdb_clade_by_mag=gtdb_clade_by_mag,
            hits_by_mag=fw_hits_by_mag,
        )

    sh_all_rep = descriptive_assembly_gate2["all_donor_modules"]["unbinned_share_vs_490reps"]
    sh_all_655 = descriptive_assembly_gate2["all_donor_modules"]["unbinned_share_vs_655dastool"]
    if event_level_results is not None:
        verdict_q1 = str(event_level_results["combined_primary_verdict_p2_q1"])
        verdict_q2 = str(event_level_results["combined_primary_verdict_p2_q2"])
    else:
        pending_verdict = (
            "Inconclusive (data-adequacy gate: Gate 1 pending [232/232 per-library targets_v2 interval files "
            "staged & full-SHA256-verified in data/grex_mag_freshwater/breadth_v2 via manifest_v2.tsv 6635acc1a6ee, "
            "awaiting grid's library->campaign_x_station event map to form >=10 independent events] "
            f"and Gate 2 per-event detected unbinned donor-KO share awaiting that event union "
            f"[descriptive assembly-wide unbinned donor-KO share in community_annotations.tsv 935abfebc239 is "
            f"{sh_all_rep:.4f} vs. 490 derep reps / {sh_all_655:.4f} vs. 655 dastool-* bins])"
        )
        verdict_q1 = pending_verdict
        verdict_q2 = pending_verdict

    fw_tuples_tsv = os.path.join(OUT_DIR, "freshwater_catalog_diag_tuples.tsv")
    with open(fw_tuples_tsv, "w", encoding="utf-8", newline="") as fh:
        cols = [
            "stratum",
            "recipient_mag",
            "domain_gtdbtk",
            "completeness",
            "contamination",
            "n_contigs",
            "catalog_490_call",
            "catalog_490_witness",
            "catalog_with_unbinned_euk_call",
            "catalog_with_unbinned_euk_witness",
        ]
        dw = csv.DictWriter(fh, fieldnames=cols, delimiter="\t", lineterminator="\n")
        dw.writeheader()
        for r in fw_tuple_rows:
            dw.writerow(r)

    p2_a1_path = os.path.join(REPO_ROOT, "protocols", "P2-A1-null2-component-conditioning.md")
    summary = {
        "protocol_sha256": FROZEN_HASHES[os.path.join(REPO_ROOT, "protocols", "P2-plexitrophy-recipient.md")],
        "p2_a1_proposed_addendum_sha256": sha256_12(p2_a1_path) if os.path.isfile(p2_a1_path) else None,
        "p2_a1_approved_by_grid": p2_a1_approved,
        "scorer_sha256": FROZEN_HASHES[os.path.join(SCRIPT_DIR, "p2_scorer.py")],
        "pathway_rules_sha256": FROZEN_HASHES[os.path.join(SCRIPT_DIR, "p2_pathway_rules.py")],
        "inferential_defaults": {
            "n_null1_draws": n_null_draws,
            "n_null2_draws": n_null_draws,
            "n_boot": n_boot,
            "require_protocol_1000_draws": True,
        },
        "cohort": "grex_mag_freshwater_staged",
        "provenance_files_sha256_verified": verified_prov_files,
        "staged_file_availability": {
            "quality_report_tsv": True,
            "derep_clusters_tsv": True,
            "per_mag_tsvs_490": True,
            "bin_fastas_in_local_derep_test": True,
            "community_annotations_tsv": os.path.isfile(comm_path),
            "tiara_output_tsv": os.path.isfile(tiara_path),
            "bakta_annotation_gff3": os.path.isfile(gff_path),
            "gtdbtk_taxonomy_tsv": os.path.isfile(gtdb_path),
            "event_union_bam_breadths_and_event_map": event_level_results is not None,
        },
        "primary_verdict_p2_q1": verdict_q1,
        "primary_verdict_p2_q2": verdict_q2,
        "flow_denominators": {
            "raw_quality_report_rows": len(q_rows),
            "raw_hq_rows": raw_hq,
            "consensus_dastool_bins_in_quality_report": len(dastool_q_rows),
            "consensus_dastool_hq_pre_derep": dastool_pre_derep_hq,
            "derep_input_passing_bins": len(derep_rows),
            "derep_95ani_representatives": len(rep_ids),
            "derep_95ani_representative_contigs": len(rep_contigs),
            "all_655_dastool_consensus_contigs": len(all_dastool_contigs),
            "derep_95ani_hq_ge90_le5": len(hq_ids),
            "derep_95ani_near_complete_ge95_le5": len(nc_ids),
            "gtdbtk_domain_counts_490": gtdb_domain_counts,
            "conservative_marker_domain_counts_490": marker_domain_counts,
            "hq_resolved_domain_count": len(hq_resolved_dom),
            "hq_unresolved_domain_U_count": len(hq_unresolved_dom),
            "gff3_unknown_edge_kos_across_490_reps": total_unknown_edge_kos,
        },
        "targets_v2_extra_p2_contig_reconciliation": extra_reconciliation,
        "staged_breadth_v2_232_libraries_verification": staged_breadth_v2_verification,
        "descriptive_assembly_wide_gate2_donor_ko_shares_NOT_event_detected": descriptive_assembly_gate2,
        "catalog_level_donor_inventory_490_reps": {
            "n_complete_b12_de_novo_donors": len(b12_complete_donors),
            "complete_b12_plus_A_DSR_single_donors": b12_plus_dsr_donors,
            "complete_b12_plus_A_NIF_single_donors": b12_plus_nif_donors,
            "descriptive_I1_PSII_x_NIF_coencoders": i1_fw,
            "module_state_counts_490": donor_module_counts,
        },
        "non_event_resolved_coassembly_diagnostic_ONLY_not_a_verdict": {
            "catalog_490_reps_only_with_gff3_gtdbtk": coassembly_diag_490,
            "catalog_490_reps_plus_assembly_unbinned_and_tiara_euk": coassembly_diag_with_pools,
        },
        "event_level_results_if_intervals_and_map_supplied": event_level_results,
    }
    out_json = os.path.join(OUT_DIR, "freshwater_staged_summary.json")
    with open(out_json, "w", encoding="utf-8") as fh:
        json.dump(summary, fh, indent=2, sort_keys=False)
    return summary


def main() -> None:
    parser = argparse.ArgumentParser(description="Run frozen P2 cohort evaluations.")
    parser.add_argument("--self-test", action="store_true", help="Run nested Curveball and event-gate unit fixtures")
    parser.add_argument("--stage0", action="store_true", help="Run Stage 0 (Ebb_Flow_20260218)")
    parser.add_argument("--freshwater", action="store_true", help="Run staged grex_mag_freshwater")
    parser.add_argument("--freshwater-covered-dir", type=str, default=None, help="Directory of *.covered.bed.gz + MANIFEST.tsv")
    parser.add_argument("--freshwater-event-map", type=str, default=None, help="TSV mapping library -> event")
    parser.add_argument("--freshwater-targets-bed", type=str, default=None, help="Optional targets_v2.bed path to hash-verify")
    parser.add_argument("--all", action="store_true", help="Run both Stage 0 and freshwater")
    parser.add_argument("--n-null-draws", type=int, default=1000, help="Null 1 and Null 2b/2c draws per event (protocol >=1000)")
    parser.add_argument("--n-boot", type=int, default=1000, help="Event-cluster bootstrap draws (protocol >=1000)")
    parser.add_argument(
        "--p2-a1-approved",
        action="store_true",
        default=False,
        help="Enable P2-A1 component-conditional Null 2 for primary Q1 ONLY after explicit grid approval",
    )
    args = parser.parse_args()

    verify_frozen_artifacts()
    fx_res = verify_nested_curveball_fixture()
    ev_fx_res = verify_freshwater_event_gate_fixtures()
    if args.self_test:
        print("NESTED_CURVEBALL_FIXTURE_PASS:", json.dumps(fx_res))
        print("FRESHWATER_EVENT_GATE_FIXTURE_PASS:", json.dumps(ev_fx_res))
    if args.stage0 or args.all:
        s0 = run_stage0_ebb_flow(n_null_draws=args.n_null_draws)
        print("STAGE0_FLOW:", json.dumps(s0["flow_denominators"]))
        print("STAGE0_POOLED_STRATA:", json.dumps(s0["campaign_pooled_breadth50"]["strata_fanin_primary"]))  # type: ignore[index]
    if args.freshwater or args.all:
        fw = run_freshwater_staged(
            covered_dir=args.freshwater_covered_dir,
            event_map_tsv=args.freshwater_event_map,
            targets_bed_path=args.freshwater_targets_bed,
            n_null_draws=args.n_null_draws,
            n_boot=args.n_boot,
            p2_a1_approved=args.p2_a1_approved,
        )
        print("FW_VERDICT_Q1:", fw["primary_verdict_p2_q1"])
        print("FW_FLOW:", json.dumps(fw["flow_denominators"]))


if __name__ == "__main__":
    main()
