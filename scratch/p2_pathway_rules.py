#!/usr/bin/env python3
"""Curated KO subunit/route rules for P2 & P4 (pure stdlib, python3 -I safe).

Fixes Trap 8 in danaSeq's kegg_module_completeness.py and implements:
  - notes/p2-incompatible-pairs.md (sha256: 0c91a742e263)
  - protocols/P2-plexitrophy-recipient.md (sha256: 828c0c3e0dde)
  - Committed/terminal key-step checks, domain-aware Bacteria/Archaea/'U' rules
    (Heal et al. 2017, doi:10.1073/pnas.1608462114), and closed-form random-break
    intactness probability with contig-end truncation -> 'unresolved'.

Every KO ID in this module is verified against /data/scratch/refdbs/kofam_db/ko_list.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Dict, FrozenSet, Iterable, List, Optional, Sequence, Set, Tuple, Union

STATE_PRESENT = "present"
STATE_ABSENT = "absent"
STATE_UNRESOLVED = "unresolved"


@dataclass(frozen=True)
class GeneHit:
    """Single annotated CDS feature from per_mag/<mag>.tsv."""
    protein_id: str
    contig_id: str
    kos: FrozenSet[str]
    ecs: FrozenSet[str] = frozenset()
    edge_status: str = "internal"  # 'internal', 'end', or 'unknown'


@dataclass(frozen=True)
class GenomeAnnotations:
    """Parsed prokaryotic MAG annotation summary.

    - end_kos: KOs occurring ONLY at first/last CDS of a contig.
    - unknown_edge_kos: KOs whose contig CDS order could not be resolved (and no
      provably internal copy exists).
    - confirmed_absent_modules: modules with Nanopore interior syntenic deletion proof.
    - false_absence_probs: optional empirical Null-1 false-absence rate per module_id.
    """
    mag_id: str
    domain: str  # 'Bacteria', 'Archaea', or 'U' (unresolved prokaryote)
    completeness: float = 100.0
    contamination: float = 0.0
    genome_size_bp: int = 0
    contig_lengths: Tuple[int, ...] = ()
    kos: FrozenSet[str] = frozenset()
    end_kos: FrozenSet[str] = frozenset()
    unknown_edge_kos: FrozenSet[str] = frozenset()
    confirmed_absent_modules: FrozenSet[str] = frozenset()
    false_absence_probs: Tuple[Tuple[str, float], ...] = ()

    def calibrated_false_absence(self, module_id: str) -> Optional[float]:
        for mid, val in self.false_absence_probs:
            if mid == module_id:
                return val
        return None

    @staticmethod
    def from_kos(
        mag_id: str,
        kos: Iterable[str],
        domain: str = "Bacteria",
        completeness: float = 100.0,
        contamination: float = 0.0,
        genome_size_bp: int = 3_000_000,
        contig_lengths: Optional[Sequence[int]] = None,
        end_kos: Iterable[str] = (),
        unknown_edge_kos: Iterable[str] = (),
        confirmed_absent_modules: Iterable[str] = (),
        false_absence_probs: Optional[Dict[str, float]] = None,
    ) -> "GenomeAnnotations":
        # Complete synthetic genomes default to one 3 Mb chromosome unless overridden
        if contig_lengths is None:
            c_lens: Tuple[int, ...] = (genome_size_bp,) if completeness >= 99.9 and genome_size_bp > 0 else ()
        else:
            c_lens = tuple(contig_lengths)
        fap = tuple(sorted((false_absence_probs or {}).items()))
        return GenomeAnnotations(
            mag_id=mag_id,
            domain=domain,
            completeness=completeness,
            contamination=contamination,
            genome_size_bp=genome_size_bp,
            contig_lengths=c_lens,
            kos=frozenset(kos),
            end_kos=frozenset(end_kos),
            unknown_edge_kos=frozenset(unknown_edge_kos),
            confirmed_absent_modules=frozenset(confirmed_absent_modules),
            false_absence_probs=fap,
        )


def p_intact(
    span_bp: float,
    contig_lengths: Sequence[int],
    completeness: float = 1.0,
    beta: float = 1.0,
) -> float:
    """Exact closed-form intact-locus detection probability under random contig breaks.

    completeness may be given on [0, 1] or [0, 100].
    """
    c = completeness / 100.0 if completeness > 1.0 else completeness
    c = min(max(float(c), 1e-3), 1.0)
    if not contig_lengths:
        return c
    total_bp = float(sum(contig_lengths))
    if total_bp <= 0.0:
        return 0.0
    g_true = total_bp / c
    s = float(span_bp)
    avail = sum(max(float(length) - s + 1.0, 0.0) for length in contig_lengths)
    f_assembled = avail / total_bp
    return max(0.0, min(1.0, (total_bp / g_true) * (1.0 - beta * (1.0 - f_assembled))))


# ---------------------------------------------------------------------------
# Boolean subunit ('+' AND) and isozyme/route (',' OR) expression nodes
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class StepEval:
    state: str  # 'present', 'absent', 'unresolved'
    partial_hit: bool = False
    end_truncated: bool = False
    unknown_edge: bool = False


class Expr:
    def eval(self, g: GenomeAnnotations) -> StepEval:
        raise NotImplementedError

    def all_kos(self) -> Set[str]:
        raise NotImplementedError


@dataclass(frozen=True)
class KO(Expr):
    ko: str

    def eval(self, g: GenomeAnnotations) -> StepEval:
        if self.ko in g.kos:
            return StepEval(
                STATE_PRESENT,
                end_truncated=(self.ko in g.end_kos),
                unknown_edge=(self.ko in g.unknown_edge_kos),
            )
        return StepEval(STATE_ABSENT)

    def all_kos(self) -> Set[str]:
        return {self.ko}


@dataclass(frozen=True)
class Complex(Expr):
    """Multi-subunit enzyme complex ('+' AND grammar).

    - r <= 2: all r subunits required.
    - r >= 3: >= ceil(0.75 * r) subunits required (unless min_subunits overrides,
      e.g., all 3 catalytic subunits for nifHDK and mcrABG) AND all active_site KOs.
    - Any partial hit (0 < present < required) is NEVER scored present; if any observed
      subunit sits at a contig end or has unknown edge status, it cannot be called absent.
    """
    subunits: Tuple[str, ...]
    active_site: Tuple[str, ...] = ()
    min_subunits: Optional[int] = None

    def __init__(
        self,
        subunits: Sequence[str],
        active_site: Sequence[str] = (),
        min_subunits: Optional[int] = None,
    ):
        object.__setattr__(self, "subunits", tuple(subunits))
        object.__setattr__(self, "active_site", tuple(active_site))
        object.__setattr__(self, "min_subunits", min_subunits)

    def required_count(self) -> int:
        r = len(self.subunits)
        if self.min_subunits is not None:
            return self.min_subunits
        return r if r <= 2 else math.ceil(0.75 * r)

    def eval(self, g: GenomeAnnotations) -> StepEval:
        present = [k for k in self.subunits if k in g.kos]
        n_pres = len(present)
        req = self.required_count()
        active_ok = all(k in g.kos for k in self.active_site)
        end_touch = any(k in g.end_kos for k in present)
        unk_edge = any(k in g.unknown_edge_kos for k in present)
        if n_pres >= req and active_ok:
            return StepEval(STATE_PRESENT, end_truncated=end_touch, unknown_edge=unk_edge)
        if n_pres == 0:
            return StepEval(STATE_ABSENT)
        return StepEval(
            STATE_ABSENT,
            partial_hit=True,
            end_truncated=end_touch,
            unknown_edge=unk_edge,
        )

    def all_kos(self) -> Set[str]:
        return set(self.subunits) | set(self.active_site)


@dataclass(frozen=True)
class AllOf(Expr):
    """Sequential/conjunctive sub-expressions inside one branch."""
    parts: Tuple[Expr, ...]

    def __init__(self, *parts: Expr):
        object.__setattr__(self, "parts", tuple(parts))

    def eval(self, g: GenomeAnnotations) -> StepEval:
        evals = [p.eval(g) for p in self.parts]
        end_touch = any(e.end_truncated for e in evals)
        unk_edge = any(e.unknown_edge for e in evals)
        if all(e.state == STATE_PRESENT for e in evals):
            return StepEval(STATE_PRESENT, end_truncated=end_touch, unknown_edge=unk_edge)
        if any(e.state == STATE_UNRESOLVED for e in evals):
            return StepEval(
                STATE_UNRESOLVED,
                partial_hit=True,
                end_truncated=end_touch,
                unknown_edge=unk_edge,
            )
        any_pres = any(e.state == STATE_PRESENT or e.partial_hit for e in evals)
        if any_pres and (end_touch or unk_edge):
            return StepEval(
                STATE_UNRESOLVED,
                partial_hit=True,
                end_truncated=end_touch,
                unknown_edge=unk_edge,
            )
        return StepEval(
            STATE_ABSENT,
            partial_hit=any_pres,
            end_truncated=end_touch,
            unknown_edge=unk_edge,
        )

    def all_kos(self) -> Set[str]:
        out: Set[str] = set()
        for p in self.parts:
            out |= p.all_kos()
        return out


@dataclass(frozen=True)
class AnyOf(Expr):
    """Alternative isozymes or non-homologous routes (',' OR grammar)."""
    options: Tuple[Expr, ...]

    def __init__(self, *options: Expr):
        object.__setattr__(self, "options", tuple(options))

    def eval(self, g: GenomeAnnotations) -> StepEval:
        evals = [o.eval(g) for o in self.options]
        pres = [e for e in evals if e.state == STATE_PRESENT]
        if pres:
            return StepEval(
                STATE_PRESENT,
                end_truncated=all(e.end_truncated for e in pres),
                unknown_edge=all(e.unknown_edge for e in pres),
            )
        if any(e.state == STATE_UNRESOLVED for e in evals):
            return StepEval(
                STATE_UNRESOLVED,
                partial_hit=any(e.partial_hit for e in evals),
                end_truncated=any(e.end_truncated for e in evals),
                unknown_edge=any(e.unknown_edge for e in evals),
            )
        return StepEval(
            STATE_ABSENT,
            partial_hit=any(e.partial_hit for e in evals),
            end_truncated=any(e.end_truncated for e in evals),
            unknown_edge=any(e.unknown_edge for e in evals),
        )

    def all_kos(self) -> Set[str]:
        out: Set[str] = set()
        for o in self.options:
            out |= o.all_kos()
        return out


@dataclass(frozen=True)
class DomainSwitch(Expr):
    """Domain-aware rule ('Bacteria' / 'Archaea' / 'Eukaryota' / 'U' unresolved)."""
    bacteria: Expr
    archaea: Expr
    eukaryota: Optional[Expr] = None

    def eval(self, g: GenomeAnnotations) -> StepEval:
        if g.domain == "Bacteria":
            return self.bacteria.eval(g)
        if g.domain == "Archaea":
            return self.archaea.eval(g)
        if g.domain == "Eukaryota":
            if self.eukaryota is not None:
                return self.eukaryota.eval(g)
            return StepEval(STATE_ABSENT)
        eb = self.bacteria.eval(g)
        ea = self.archaea.eval(g)
        if (
            eb.state == ea.state
            and eb.end_truncated == ea.end_truncated
            and eb.unknown_edge == ea.unknown_edge
        ):
            return eb
        return StepEval(
            STATE_UNRESOLVED,
            partial_hit=(eb.partial_hit or ea.partial_hit or eb.state != ea.state),
            end_truncated=(eb.end_truncated or ea.end_truncated),
            unknown_edge=(eb.unknown_edge or ea.unknown_edge),
        )

    def all_kos(self) -> Set[str]:
        out = self.bacteria.all_kos() | self.archaea.all_kos()
        if self.eukaryota is not None:
            out |= self.eukaryota.all_kos()
        return out


PROKARYOTIC_DOMAINS: FrozenSet[str] = frozenset({"Bacteria", "Archaea"})
ALL_CELLULAR_DOMAINS: FrozenSet[str] = frozenset({"Bacteria", "Archaea", "Eukaryota"})


def estimate_false_absence_prob(rule: "PathwayRule", g: GenomeAnnotations) -> float:
    """Calibrated upper bound on P(called absent | pathway truly in genome).

    Priority:
      1. Biologically excluded domain (`g.domain` known and not in `rule.allowed_domains`) -> 0.0.
      2. Confirmed Nanopore interior syntenic deletion -> 0.0.
      3. Empirical detectability calibration on `g` (`false_absence_probs`) -> calibrated value.
      4. For `Eukaryota` without explicit eukaryotic calibration -> 1.0 (prokaryotic CheckM2
         completeness and intronless operon `p_intact` do not calibrate spliced MetaEuk calls;
         Astra 248978604545196032).
      5. For prokaryotes (`Bacteria`/`Archaea`) with `contig_lengths`: closed-form bound from
         `(contig_lengths, completeness)`; else 1.0.
    """
    if g.domain in ALL_CELLULAR_DOMAINS and g.domain not in rule.allowed_domains:
        return 0.0
    if rule.module_id in g.confirmed_absent_modules:
        return 0.0
    emp = g.calibrated_false_absence(rule.module_id)
    if emp is not None:
        return float(emp)
    if g.domain == "Eukaryota":
        return 1.0
    if not g.contig_lengths:
        return 1.0
    c = g.completeness / 100.0 if g.completeness > 1.0 else g.completeness
    c = min(max(float(c), 1e-3), 1.0)
    f_locus = p_intact(rule.span_bp_est, g.contig_lengths, completeness=1.0)
    f_gene = p_intact(1200.0, g.contig_lengths, completeness=1.0)
    p_false = (1.0 - c) * f_locus + (1.0 - f_locus) * (1.0 - c * f_gene)
    return max(0.0, min(1.0, p_false))


@dataclass(frozen=True)
class PathwayEval:
    module_id: str
    state: str  # 'present', 'absent', 'unresolved'
    step_fraction: float
    n_present: int
    n_steps: int
    key_ok: bool
    partial_complex: bool
    end_truncated: bool
    unknown_edge: bool
    p_intact_est: float
    p_false_absence: float


@dataclass(frozen=True)
class PathwayRule:
    """Multi-step KO pathway or enzyme-complex definition with committed key steps."""
    module_id: str
    name: str
    regime: str  # 'O2_plus', 'O2_minus', 'C0', 'salvage', 'assimilation'
    steps: Tuple[Expr, ...]
    key_indices: Tuple[int, ...] = ()
    present_min_frac: float = 0.80
    absent_max_frac: float = 0.25
    span_bp_est: int = 4000
    allowed_domains: FrozenSet[str] = ALL_CELLULAR_DOMAINS

    def evaluate(
        self,
        g: GenomeAnnotations,
        max_false_absence_prob: Optional[float] = None,
    ) -> PathwayEval:
        if g.domain in ALL_CELLULAR_DOMAINS and g.domain not in self.allowed_domains:
            return PathwayEval(
                module_id=self.module_id,
                state=STATE_ABSENT,
                step_fraction=0.0,
                n_present=0,
                n_steps=len(self.steps),
                key_ok=False,
                partial_complex=False,
                end_truncated=False,
                unknown_edge=False,
                p_intact_est=0.0,
                p_false_absence=0.0,
            )
        evals = [s.eval(g) for s in self.steps]
        n_steps = len(evals)
        n_pres = sum(1 for e in evals if e.state == STATE_PRESENT)
        n_unres = sum(1 for e in evals if e.state == STATE_UNRESOLVED)
        frac = n_pres / n_steps if n_steps else 0.0
        key_ok = all(evals[i].state == STATE_PRESENT for i in self.key_indices)
        any_partial = any(e.partial_hit for e in evals)
        hit_evals = [e for e in evals if e.state == STATE_PRESENT or e.partial_hit]
        end_trunc = any(e.end_truncated for e in hit_evals)
        unk_edge = any(e.unknown_edge for e in hit_evals)
        pi = p_intact(self.span_bp_est, g.contig_lengths, g.completeness)
        p_fa = estimate_false_absence_prob(self, g)

        if frac >= self.present_min_frac and key_ok:
            state = STATE_PRESENT
        elif n_unres > 0:
            state = STATE_UNRESOLVED
        elif frac == 0.0 and not any_partial:
            state = STATE_ABSENT
            if max_false_absence_prob is not None and p_fa >= max_false_absence_prob:
                state = STATE_UNRESOLVED
        elif frac <= self.absent_max_frac and not end_trunc and not unk_edge:
            if any_partial and n_steps < 4:
                state = STATE_UNRESOLVED
            else:
                state = STATE_ABSENT
                if max_false_absence_prob is not None and p_fa >= max_false_absence_prob:
                    state = STATE_UNRESOLVED
        else:
            state = STATE_UNRESOLVED

        return PathwayEval(
            module_id=self.module_id,
            state=state,
            step_fraction=round(frac, 4),
            n_present=n_pres,
            n_steps=n_steps,
            key_ok=key_ok,
            partial_complex=any_partial,
            end_truncated=end_trunc,
            unknown_edge=unk_edge,
            p_intact_est=round(pi, 4),
            p_false_absence=round(p_fa, 4),
        )

    def all_kos(self) -> Set[str]:
        out: Set[str] = set()
        for s in self.steps:
            out |= s.all_kos()
        return out


# ---------------------------------------------------------------------------
# Curated P2 Pathway & Handoff Dictionary (verified against kofam_db/ko_list)
# ---------------------------------------------------------------------------

PATHWAY_RULES: Dict[str, PathwayRule] = {
    # === A. O2-evolving or obligate O2-cosubstrate modules (O_+) ===
    "P_PSII": PathwayRule(
        module_id="P_PSII",
        name="Oxygenic Photosystem II (psbAD + psbBC)",
        regime="O2_plus",
        steps=(
            Complex(("K02703", "K02706")),  # psbA + psbD reaction center
            Complex(("K02704", "K02705")),  # psbB + psbC core antenna
        ),
        key_indices=(0, 1),
        span_bp_est=6000,
    ),
    "P_BLUB_DMB": PathwayRule(
        module_id="P_BLUB_DMB",
        name="Aerobic DMB synthase BluB (free lower-ligand precursor)",
        regime="O2_plus",
        steps=(KO("K04719"),),
        key_indices=(0,),
        span_bp_est=1000,
        allowed_domains=PROKARYOTIC_DOMAINS,
    ),
    "P_BLUB": PathwayRule(
        module_id="P_BLUB",
        name="Aerobic DMB synthesis + activation (Archaea BluB without bacterial CobT)",
        regime="O2_plus",
        steps=(
            DomainSwitch(
                bacteria=Complex(("K04719", "K00768")),  # bluB + cobT in Bacteria
                archaea=KO("K04719"),  # Thaumarchaeota lack bacterial cobT (Heal 2017)
            ),
        ),
        key_indices=(0,),
        span_bp_est=2200,
        allowed_domains=PROKARYOTIC_DOMAINS,
    ),
    "P_B12_AER": PathwayRule(
        module_id="P_B12_AER",
        name="Aerobic corrinoid ring synthesis (cobG + cobNST + ring steps)",
        regime="O2_plus",
        steps=(
            AnyOf(KO("K03394"), KO("K13540")),  # cobI precorrin-2 methyltransferase
            KO("K02229"),  # cobG O2-dependent precorrin-3B synthase [EC:1.14.13.83]
            AnyOf(KO("K05934"), KO("K13540"), KO("K05936"), KO("K02228")),  # cobJ/M/F
            KO("K06042"),  # cobH precorrin-8X methylmutase
            Complex(("K02230", "K09882", "K09883"), active_site=("K02230",), min_subunits=3),  # cobNST
            AnyOf(KO("K02232"), KO("K02227")),  # cobQ / cobD late amidation/aminopropanol
        ),
        key_indices=(1, 4),  # cobG + complete cobNST aerobic cobaltochelatase required
        span_bp_est=12000,
        allowed_domains=PROKARYOTIC_DOMAINS,
    ),
    "P_CU_MMO": PathwayRule(
        module_id="P_CU_MMO",
        name="Aerobic ammonia / methane monooxygenase (amoABC/pmoABC or mmoXYBZ)",
        regime="O2_plus",
        steps=(
            AnyOf(
                Complex(("K10944", "K10945", "K10946"), active_site=("K10944",), min_subunits=3),
                Complex(("K16157", "K16158", "K16159", "K16161"), active_site=("K16157",)),
            ),
        ),
        key_indices=(0,),
        span_bp_est=3500,
        allowed_domains=PROKARYOTIC_DOMAINS,
    ),
    "P_AROM_OXY": PathwayRule(
        module_id="P_AROM_OXY",
        name="Aerobic ring-hydroxylating / catechol dioxygenases",
        regime="O2_plus",
        steps=(
            AnyOf(
                Complex(("K16243", "K16244", "K16242"), active_site=("K16243",), min_subunits=3),
                KO("K03381"),  # catA
                KO("K00446"),  # catE
            ),
        ),
        key_indices=(0,),
        span_bp_est=3000,
    ),
    "P_THIO_AER": PathwayRule(
        module_id="P_THIO_AER",
        name="Aerobic thiazole synthase branch (thiO + thiG)",
        regime="O2_plus",
        steps=(KO("K03153"), KO("K03149")),
        key_indices=(0, 1),
        span_bp_est=2000,
        allowed_domains=PROKARYOTIC_DOMAINS,
    ),

    # === B. O2-inactivated metalloenzyme & radical-SAM anaerobic modules (O_-) ===
    "A_NIF": PathwayRule(
        module_id="A_NIF",
        name="MoFe nitrogenase (nifH + nifD + nifK, all 3 required)",
        regime="O2_minus",
        steps=(Complex(("K02588", "K02586", "K02591"), active_site=("K02586", "K02591"), min_subunits=3),),
        key_indices=(0,),
        span_bp_est=4500,
        allowed_domains=PROKARYOTIC_DOMAINS,
    ),
    "A_WL": PathwayRule(
        module_id="A_WL",
        name="Wood-Ljungdahl anaerobic CODH/ACS complex",
        regime="O2_minus",
        steps=(
            KO("K00198"),  # cooS/acsA
            Complex(("K14138", "K00197", "K00194"), active_site=("K14138",)),  # acsBCD
        ),
        key_indices=(0, 1),
        span_bp_est=5500,
        allowed_domains=PROKARYOTIC_DOMAINS,
    ),
    "A_DSR": PathwayRule(
        module_id="A_DSR",
        name="Dissimilatory sulfate reduction (aprAB + dsrAB)",
        regime="O2_minus",
        steps=(
            Complex(("K00394", "K00395")),  # aprAB
            Complex(("K11180", "K11181")),  # dsrAB
        ),
        key_indices=(0, 1),
        span_bp_est=5000,
        allowed_domains=PROKARYOTIC_DOMAINS,
    ),
    "A_MCR": PathwayRule(
        module_id="A_MCR",
        name="Methanogenesis methyl-CoM reductase (mcrABG, all 3 required)",
        regime="O2_minus",
        steps=(Complex(("K00399", "K00401", "K00402"), active_site=("K00399",), min_subunits=3),),
        key_indices=(0,),
        span_bp_est=4500,
        allowed_domains=PROKARYOTIC_DOMAINS,
    ),
    "A_FEFE_HYD": PathwayRule(
        module_id="A_FEFE_HYD",
        name="Anaerobic [FeFe] hydrogenase (hydAB)",
        regime="O2_minus",
        steps=(Complex(("K00533", "K00534")),),
        key_indices=(0,),
        span_bp_est=2800,
    ),
    "A_NOS_NOR": PathwayRule(
        module_id="A_NOS_NOR",
        name="Denitrification terminal steps (norBC / nosZ)",
        regime="O2_minus",
        steps=(AnyOf(Complex(("K04561", "K02305")), KO("K00376")),),
        key_indices=(0,),
        span_bp_est=2500,
        allowed_domains=PROKARYOTIC_DOMAINS,
    ),
    "A_BZA": PathwayRule(
        module_id="A_BZA",
        name="Anaerobic benzimidazole synthase (bzaAB/F + bzaCDE)",
        regime="O2_minus",
        steps=(
            AnyOf(KO("K22465"), KO("K22466")),
            AnyOf(KO("K25497"), KO("K25498"), KO("K25499")),
        ),
        key_indices=(0, 1),
        span_bp_est=3500,
        allowed_domains=PROKARYOTIC_DOMAINS,
    ),
    "A_B12_ANA": PathwayRule(
        module_id="A_B12_ANA",
        name="Oxygen-independent corrinoid ring synthesis (cbiK/X + cbiG + cbi steps)",
        regime="O2_minus",
        steps=(
            AnyOf(KO("K02190"), KO("K03795")),  # cbiK / cbiX sirohydrochlorin cobaltochelatase
            AnyOf(KO("K02228"), KO("K03394")),  # cbiL / cbiH early methylation
            AnyOf(KO("K02189"), KO("K13541")),  # cbiG anaerobic ring-contraction hydrolase
            AnyOf(KO("K02188"), KO("K02191"), KO("K00595")),  # cbiD / cbiET
            KO("K06042"),  # cbiC cobalt-precorrin-8 methylmutase
            AnyOf(KO("K02232"), KO("K02227")),  # cbiP / cbiB cobyric/cobinamide-P synthesis
        ),
        key_indices=(0, 2, 5),  # cbiK/X + cbiG + late cobinamide-P/cobyric step
        span_bp_est=11000,
        allowed_domains=PROKARYOTIC_DOMAINS,
    ),

    # === Salvage, Remodeling & Assimilation Modules ===
    "B12_LOWER_ATTACH": PathwayRule(
        module_id="B12_LOWER_ATTACH",
        name="Cobinamide + DMB lower-ligand remodeling/attachment (CobT/CobU|CobY/CobS/CobC)",
        regime="salvage",
        steps=(
            DomainSwitch(
                bacteria=KO("K00768"),  # bacterial cobT phosphoribosyltransferase
                archaea=AnyOf(KO("K19712"), KO("K00768")),  # archaeal cobY/salvage context
                eukaryota=KO("K00768"),  # algal COBT remodeling (Helliwell et al. 2016)
            ),
            DomainSwitch(
                bacteria=KO("K02231"),  # cobU/cobP
                archaea=AnyOf(KO("K19712"), KO("K02231")),  # archaeal cobY [EC:2.7.7.62]
                eukaryota=KO("K02231"),
            ),
            KO("K02233"),  # cobS adenosylcobinamide-GDP ribazoletransferase [EC:2.7.8.26]
            AnyOf(KO("K02226"), KO("K22316")),  # cobC / cobZ alpha-ribazole phosphatase
        ),
        key_indices=(2,),  # cobS ribazoletransferase committed attachment enzyme
        span_bp_est=3800,
    ),
    "ASSIM_SO4": PathwayRule(
        module_id="ASSIM_SO4",
        name="Assimilatory sulfate reduction (PAPS cysDN/sat+cysCH or APS K05907 + cysJI/sir)",
        regime="assimilation",
        steps=(
            AnyOf(KO("K00958"), KO("K00955"), KO("K13811"), Complex(("K00956", "K00957"))),  # sat / cysDN
            AnyOf(
                AllOf(AnyOf(KO("K00860"), KO("K00955"), KO("K13811")), KO("K00390")),  # cysC + cysH
                KO("K05907"),  # direct APS reductase [EC:1.8.4.9]
            ),
            AnyOf(Complex(("K00380", "K00381")), KO("K00392")),  # cysJI heterodimer OR sir
        ),
        key_indices=(1, 2),
        span_bp_est=4500,
    ),
    "ASSIM_NO3": PathwayRule(
        module_id="ASSIM_NO3",
        name="Assimilatory nitrate reduction (narB/nasAB -> nirA/NIT-6/nirBD)",
        regime="assimilation",
        steps=(
            AnyOf(KO("K00367"), KO("K10534"), Complex(("K00372", "K00360"))),  # nitrate -> nitrite
            AnyOf(KO("K00366"), KO("K17877"), Complex(("K00362", "K00363"))),  # nitrite -> NH4+
        ),
        key_indices=(0, 1),
        span_bp_est=4000,
    ),
    "THI_HMP": PathwayRule(
        module_id="THI_HMP",
        name="Thiamin HMP pyrimidine branch (thiC/THIC + thiD/THI6)",
        regime="C0",
        steps=(KO("K03147"), AnyOf(KO("K00941"), KO("K14153"))),
        key_indices=(0,),
        span_bp_est=2500,
    ),
    "THI_HET": PathwayRule(
        module_id="THI_HET",
        name="Thiamin HET thiazole branch (bacterial thiG+thiH/O or euk/arch THI4)",
        regime="C0",
        steps=(
            AnyOf(
                AllOf(KO("K03149"), AnyOf(KO("K03150"), KO("K03153"))),
                KO("K03146"),  # eukaryotic/archaeal cysteine-dependent THI4 [EC:2.4.2.60]
                KO("K22699"),  # archaeal sulfide-dependent THI4 [EC:2.4.2.59]
            ),
        ),
        key_indices=(0,),
        span_bp_est=2500,
    ),
    "THI_COUPLE": PathwayRule(
        module_id="THI_COUPLE",
        name="Thiamin coupling & phosphorylation (thiE/THI6 + thiL/THI80)",
        regime="salvage",
        steps=(
            AnyOf(KO("K00788"), KO("K14153"), KO("K14154")),
            AnyOf(KO("K00946"), KO("K00949")),
        ),
        key_indices=(0,),
        span_bp_est=1800,
    ),
    "C_THI_SPLIT": PathwayRule(
        module_id="C_THI_SPLIT",
        name="Complete de novo thiamin synthesis (HMP + HET/THI4 + thiEL/THI80)",
        regime="C0",
        steps=(
            KO("K03147"),  # thiC / THIC
            AnyOf(KO("K00941"), KO("K14153")),  # thiD / bifunctional THID-THIE
            AnyOf(
                AllOf(KO("K03149"), AnyOf(KO("K03150"), KO("K03153"))),  # bacterial thiG + thiH/O
                KO("K03146"),  # eukaryotic/archaeal THI4 [EC:2.4.2.60]
                KO("K22699"),  # sulfide-dependent THI4 [EC:2.4.2.59]
            ),
            AnyOf(KO("K00788"), KO("K14153"), KO("K14154")),  # thiE / bifunctional THIE
            AnyOf(KO("K00946"), KO("K00949")),  # thiL / eukaryotic TPK1 (THI80)
        ),
        key_indices=(0, 2, 3),
        span_bp_est=5500,
    ),

    # === C. Matched Compatible Redox-Neutral Control Pool (C_0) ===
    "C_HIS": PathwayRule(
        module_id="C_HIS",
        name="Histidine biosynthesis (hisG + hisA + hisB + hisD)",
        regime="C0",
        steps=(KO("K00765"), KO("K01814"), KO("K01693"), KO("K00013")),
        key_indices=(0, 3),
        span_bp_est=4500,
    ),
    "C_TRP": PathwayRule(
        module_id="C_TRP",
        name="Tryptophan biosynthesis (trpEG + trpD + trpC + trpAB)",
        regime="C0",
        steps=(
            Complex(("K01657", "K01658")),
            KO("K00766"),
            KO("K01609"),
            Complex(("K01695", "K01696")),
        ),
        key_indices=(0, 3),
        span_bp_est=5500,
    ),
    "C_LEU": PathwayRule(
        module_id="C_LEU",
        name="Leucine biosynthesis (leuA + leuCD + leuB)",
        regime="C0",
        steps=(KO("K01649"), Complex(("K01703", "K01704")), KO("K00052")),
        key_indices=(0, 2),
        span_bp_est=4200,
    ),
    "C_ILV": PathwayRule(
        module_id="C_ILV",
        name="Valine/isoleucine biosynthesis (ilvBH + ilvC + ilvD + ilvE)",
        regime="C0",
        steps=(AnyOf(KO("K01652"), KO("K01653")), KO("K00053"), KO("K01687"), KO("K00826")),
        key_indices=(2, 3),
        span_bp_est=4500,
    ),
    "C_LYS_DAP": PathwayRule(
        module_id="C_LYS_DAP",
        name="Lysine DAP biosynthesis (dapA + dapB + lysA)",
        regime="C0",
        steps=(KO("K01714"), KO("K00215"), KO("K01586")),
        key_indices=(0, 2),
        span_bp_est=3500,
    ),
    "C_ARG": PathwayRule(
        module_id="C_ARG",
        name="Arginine biosynthesis (argF + argG + argH)",
        regime="C0",
        steps=(KO("K00611"), KO("K01940"), KO("K01755")),
        key_indices=(1, 2),
        span_bp_est=3500,
    ),
    "C_PRO": PathwayRule(
        module_id="C_PRO",
        name="Proline biosynthesis (proB + proA + proC)",
        regime="C0",
        steps=(KO("K00931"), KO("K00147"), KO("K00286")),
        key_indices=(0, 2),
        span_bp_est=3200,
    ),
    "C_RIB": PathwayRule(
        module_id="C_RIB",
        name="Riboflavin biosynthesis (ribA + ribBA + ribH + ribE)",
        regime="C0",
        steps=(KO("K01497"), KO("K14652"), KO("K00794"), KO("K00793")),
        key_indices=(2, 3),
        span_bp_est=3500,
    ),
    "C_PAN": PathwayRule(
        module_id="C_PAN",
        name="Pantothenate biosynthesis (panB + panE + panC)",
        regime="C0",
        steps=(KO("K00606"), KO("K00077"), KO("K01918")),
        key_indices=(0, 2),
        span_bp_est=3000,
    ),
    "C_BIO": PathwayRule(
        module_id="C_BIO",
        name="Biotin ring closure (bioF + bioA + bioB)",
        regime="C0",
        steps=(KO("K00652"), KO("K00833"), KO("K01012")),
        key_indices=(2,),
        span_bp_est=3200,
    ),
}

C0_CONTROL_MODULES: Tuple[str, ...] = (
    "C_HIS",
    "C_TRP",
    "C_LEU",
    "C_ILV",
    "C_LYS_DAP",
    "C_ARG",
    "C_PRO",
    "C_RIB",
    "C_PAN",
    "C_BIO",
)

# Section 4 intracellular compartmentalization & respiratory-protection KOs
SPATIAL_COMPARTMENT_KOS: FrozenSet[str] = frozenset(
    {
        "K25172",  # hgd heterocyst glycolipid synthase
        "K27260",
        "K27263",
        "K27266",
        "K27267",
        "K27268",
        "K27269",
        "K27270",
        "K27271",
        "K27272",
        "K04027",  # eutM BMC shell
    }
)

KAI_ABC_COMPLEX = Complex(("K08480", "K08481", "K08482"), min_subunits=3)
CYD_AB_COMPLEX = Complex(("K00425", "K00426"))
CCO_NO_COMPLEX = Complex(("K00404", "K00405"))


def classify_coencoder_bypass(g: GenomeAnnotations) -> str:
    """Classify a single-genome O_+ x O_- co-encoder per Section 4 of p2-incompatible-pairs.md."""
    if g.contamination > 5.0:
        return "chimeric_candidate"
    if g.kos & SPATIAL_COMPARTMENT_KOS:
        return "S_spatial"
    if (
        KAI_ABC_COMPLEX.eval(g).state == STATE_PRESENT
        or CYD_AB_COMPLEX.eval(g).state == STATE_PRESENT
        or CCO_NO_COMPLEX.eval(g).state == STATE_PRESENT
    ):
        return "S_temp_resp"
    return "S_unexplained"


def evaluate_b12_de_novo(
    g: GenomeAnnotations,
    max_false_absence_prob: Optional[float] = 0.05,
) -> str:
    """Return 'present' | 'absent' | 'unresolved' for complete de novo cobalamin synthesis.

    Note: encoding salvage `cobT` (`K00768`) without `bluB` (`K04719`) or `bza` is
    salvage/remodeling capacity, not partial de novo lower-ligand synthesis.
    """
    aer_ring = PATHWAY_RULES["P_B12_AER"].evaluate(g, max_false_absence_prob).state
    ana_ring = PATHWAY_RULES["A_B12_ANA"].evaluate(g, max_false_absence_prob).state
    blub_full = PATHWAY_RULES["P_BLUB"].evaluate(g, max_false_absence_prob).state
    blub_dmb = PATHWAY_RULES["P_BLUB_DMB"].evaluate(g, max_false_absence_prob).state
    bza = PATHWAY_RULES["A_BZA"].evaluate(g, max_false_absence_prob).state

    ring_pres = (aer_ring == STATE_PRESENT) or (ana_ring == STATE_PRESENT)
    ring_abs = (aer_ring == STATE_ABSENT) and (ana_ring == STATE_ABSENT)
    lig_pres = (blub_full == STATE_PRESENT) or (bza == STATE_PRESENT)
    lig_denovo_abs = (blub_dmb == STATE_ABSENT) and (bza == STATE_ABSENT)

    if ring_pres and lig_pres:
        return STATE_PRESENT
    if ring_abs and lig_denovo_abs:
        return STATE_ABSENT
    return STATE_UNRESOLVED


def evaluate_recipient_b12_auxotrophy(
    g: GenomeAnnotations,
    max_false_absence_prob: Optional[float] = 0.05,
) -> str:
    """Recipient eligibility for cobalamin: metH (K00548) without metE (K00549) and lacking de novo B12."""
    has_meth = "K00548" in g.kos
    has_mete = "K00549" in g.kos
    if has_mete or not has_meth:
        return "not_dependent"
    # Also check detectability of metE (1200 bp single CDS) before calling metH-only auxotrophy
    if max_false_absence_prob is not None and "K00549" not in g.confirmed_absent_modules:
        emp_mete = g.calibrated_false_absence("K00549")
        if emp_mete is not None:
            p_fa_mete = emp_mete
        elif not g.contig_lengths:
            p_fa_mete = 1.0
        else:
            c = g.completeness / 100.0 if g.completeness > 1.0 else g.completeness
            c = min(max(float(c), 1e-3), 1.0)
            f_gene = p_intact(2200.0, g.contig_lengths, completeness=1.0)
            p_fa_mete = (1.0 - c) * f_gene + (1.0 - f_gene) * (1.0 - c * f_gene)
        if p_fa_mete >= max_false_absence_prob:
            return STATE_UNRESOLVED
    denovo = evaluate_b12_de_novo(g, max_false_absence_prob=max_false_absence_prob)
    if denovo == STATE_PRESENT:
        return STATE_PRESENT
    if denovo == STATE_ABSENT:
        return STATE_ABSENT
    return STATE_UNRESOLVED


def all_referenced_kos() -> Set[str]:
    kos: Set[str] = {"K00548", "K00549"} | set(SPATIAL_COMPARTMENT_KOS)
    kos |= KAI_ABC_COMPLEX.all_kos() | CYD_AB_COMPLEX.all_kos() | CCO_NO_COMPLEX.all_kos()
    for rule in PATHWAY_RULES.values():
        kos |= rule.all_kos()
    return kos


def verify_all_kos_in_ko_list(ko_list_path: str = "/data/scratch/refdbs/kofam_db/ko_list") -> List[str]:
    """Return sorted list of any referenced KO IDs missing from kofam_db/ko_list."""
    with open(ko_list_path, "r", encoding="utf-8") as fh:
        next(fh, None)
        valid = {line.split("\t", 1)[0].strip() for line in fh if line.strip()}
    return sorted(all_referenced_kos() - valid)
