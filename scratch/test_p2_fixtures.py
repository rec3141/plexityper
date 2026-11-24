#!/usr/bin/env python3
"""Synthetic unit-test fixture suite for P2 pathway rules, P-all decider, and posteriors.

Runs under `python3 -I scratch/test_p2_fixtures.py` with zero external packages
and WITHOUT opening any staged project `per_mag` outcome files (keeping P2 blind).
"""

from __future__ import annotations

import itertools
import math
import os
import random
import sys
import unittest

sys.path.insert(0, os.path.dirname(os.path.realpath(__file__)))

from p2_pathway_rules import (  # noqa: E402
    PATHWAY_RULES,
    STATE_ABSENT,
    STATE_PRESENT,
    STATE_UNRESOLVED,
    GenomeAnnotations,
    classify_coencoder_bypass,
    p_intact,
    verify_all_kos_in_ko_list,
)
from p2_scorer import (  # noqa: E402
    CALL_AUTONOMOUS,
    CALL_D0_AMBIENT,
    CALL_D1_OR_CHAIN,
    CALL_NOT_ELIGIBLE,
    CALL_PLEXITROPHIC,
    CALL_UNRESOLVED_AMBIG,
    CALL_UNRESOLVED_NO_DONOR,
    SCENARIO_E_MINIMAL,
    SCENARIO_E_PLUS_DIN,
    OrganismProfile,
    StratumCounts,
    SupplyScenario,
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
    enumerate_minimal_wirings,
    evaluate_recipient_requirements_in_event,
    filter_derep_prokaryotic_mags,
    is_metabolite_dep_acyclic,
    manski_log_odds_bounds,
    observed_minus_null_id_bounds,
    on_simple_cycle_ge_3,
    parse_gff_cds_order,
    parse_per_mag_rows,
    tarjan_sccs,
)

AER_RING_KOS = {"K03394", "K02229", "K05934", "K06042", "K02230", "K09882", "K09883", "K02232"}
ANA_RING_KOS = {"K02190", "K02228", "K02189", "K02188", "K06042", "K02232"}
BLUB_BAC_KOS = {"K04719", "K00768"}
BZA_KOS = {"K22465", "K25497"}
B12_ATTACH_BAC_KOS = {"K00768", "K02231", "K02233", "K02226"}


class TestP2Fixtures(unittest.TestCase):
    def test_ko_list_verification(self) -> None:
        ko_list = "/data/scratch/refdbs/kofam_db/ko_list"
        self.assertTrue(os.path.exists(ko_list), "kofam_db/ko_list must exist")
        missing = verify_all_kos_in_ko_list(ko_list)
        self.assertEqual(missing, [], f"Unrecognized KOs in rules: {missing}")

    def test_suite1_positive_pathway_and_domain_fixtures(self) -> None:
        g_aer = GenomeAnnotations.from_kos("dastool-b1", AER_RING_KOS | BLUB_BAC_KOS)
        self.assertEqual(PATHWAY_RULES["P_B12_AER"].evaluate(g_aer).state, STATE_PRESENT)
        self.assertEqual(PATHWAY_RULES["P_BLUB"].evaluate(g_aer).state, STATE_PRESENT)

        g_ana = GenomeAnnotations.from_kos("dastool-b2", ANA_RING_KOS | BZA_KOS)
        self.assertEqual(PATHWAY_RULES["A_B12_ANA"].evaluate(g_ana).state, STATE_PRESENT)
        self.assertEqual(PATHWAY_RULES["A_BZA"].evaluate(g_ana).state, STATE_PRESENT)

        # Archaeal Thaumarchaeota exception (Heal et al. 2017): bluB without bacterial cobT
        g_thaum = GenomeAnnotations.from_kos(
            "dastool-arch1", ANA_RING_KOS | {"K04719"}, domain="Archaea"
        )
        self.assertEqual(PATHWAY_RULES["P_BLUB"].evaluate(g_thaum).state, STATE_PRESENT)
        g_bac_no_cobt = GenomeAnnotations.from_kos(
            "dastool-bac_nocobt", ANA_RING_KOS | {"K04719"}, domain="Bacteria"
        )
        self.assertNotEqual(PATHWAY_RULES["P_BLUB"].evaluate(g_bac_no_cobt).state, STATE_PRESENT)
        g_unk_dom = GenomeAnnotations.from_kos(
            "dastool-unk", ANA_RING_KOS | {"K04719"}, domain="U"
        )
        self.assertEqual(PATHWAY_RULES["P_BLUB"].evaluate(g_unk_dom).state, STATE_UNRESOLVED)

        g_nif = GenomeAnnotations.from_kos("dastool-nif", {"K02588", "K02586", "K02591"})
        self.assertEqual(PATHWAY_RULES["A_NIF"].evaluate(g_nif).state, STATE_PRESENT)
        g_thi = GenomeAnnotations.from_kos(
            "dastool-thi", {"K03147", "K00941", "K03149", "K03153", "K00788", "K00946"}
        )
        self.assertEqual(PATHWAY_RULES["C_THI_SPLIT"].evaluate(g_thi).state, STATE_PRESENT)

        # Trap 8 fixes on ASSIM_NO3 (M00531) and ASSIM_SO4 (M00176)
        g_no3_ok = GenomeAnnotations.from_kos("dastool-no3", {"K00367", "K00366"})
        self.assertEqual(PATHWAY_RULES["ASSIM_NO3"].evaluate(g_no3_ok).state, STATE_PRESENT)
        g_nas_only = GenomeAnnotations.from_kos("dastool-nasonly", {"K00372", "K00360"})
        self.assertEqual(PATHWAY_RULES["ASSIM_NO3"].evaluate(g_nas_only).state, STATE_UNRESOLVED)

        g_so4_a = GenomeAnnotations.from_kos(
            "dastool-so4a", {"K00956", "K00957", "K00860", "K00390", "K00392"}
        )
        self.assertEqual(PATHWAY_RULES["ASSIM_SO4"].evaluate(g_so4_a).state, STATE_PRESENT)
        g_so4_b = GenomeAnnotations.from_kos(
            "dastool-so4b", {"K00958", "K05907", "K00380", "K00381"}
        )
        self.assertEqual(PATHWAY_RULES["ASSIM_SO4"].evaluate(g_so4_b).state, STATE_PRESENT)

    def test_suite2_partial_complex_negative_fixtures(self) -> None:
        for kos in ({"K10944"}, {"K10944", "K16157"}):
            g = GenomeAnnotations.from_kos("dastool-mmo-part", kos)
            ev = PATHWAY_RULES["P_CU_MMO"].evaluate(g)
            self.assertNotEqual(ev.state, STATE_PRESENT)
            self.assertTrue(ev.partial_complex)

        g_nifh = GenomeAnnotations.from_kos("dastool-nifh", {"K02588"})
        self.assertEqual(PATHWAY_RULES["A_NIF"].evaluate(g_nifh).state, STATE_UNRESOLVED)

        g_dsr_half = GenomeAnnotations.from_kos("dastool-dsr", {"K00394", "K11180"})
        self.assertNotEqual(PATHWAY_RULES["A_DSR"].evaluate(g_dsr_half).state, STATE_PRESENT)

        g_so4_half = GenomeAnnotations.from_kos("dastool-cysj", {"K00958", "K05907", "K00380"})
        self.assertNotEqual(PATHWAY_RULES["ASSIM_SO4"].evaluate(g_so4_half).state, STATE_PRESENT)

    def test_suite3_contig_end_shuffle_and_detectability_guard(self) -> None:
        # Astra 248964604793110528: KO-bearing CDS between two blank-KO CDSs across all 6 row orders
        base_rows = [
            {"protein_id": "ctgA_00001", "contig_id": "ctgA", "KO": "", "EC": ""},
            {"protein_id": "ctgA_00002", "contig_id": "ctgA", "KO": "ko:K00765", "EC": ""},
            {"protein_id": "ctgA_00003", "contig_id": "ctgA", "KO": "", "EC": ""},
        ]
        for perm in itertools.permutations(base_rows):
            g_perm, hits = parse_per_mag_rows(
                list(perm), mag_id="dastool-mid", completeness=98.0, contig_lengths=[500_000]
            )
            hit_map = {h.protein_id: h.edge_status for h in hits}
            self.assertEqual(hit_map["ctgA_00002"], "internal")
            ev = PATHWAY_RULES["C_HIS"].evaluate(g_perm, max_false_absence_prob=0.05)
            self.assertEqual(ev.state, STATE_ABSENT)
            self.assertFalse(ev.end_truncated)
            self.assertFalse(ev.unknown_edge)

        # Placed at first CDS (ctgA_00001) -> 'end' across all 6 row orders -> unresolved
        end_rows = [
            {"protein_id": "ctgA_00001", "contig_id": "ctgA", "KO": "K00765", "EC": ""},
            {"protein_id": "ctgA_00002", "contig_id": "ctgA", "KO": "", "EC": ""},
            {"protein_id": "ctgA_00003", "contig_id": "ctgA", "KO": "", "EC": ""},
        ]
        for perm in itertools.permutations(end_rows):
            g_end, _ = parse_per_mag_rows(
                list(perm), mag_id="dastool-end", completeness=98.0, contig_lengths=[500_000]
            )
            ev_end = PATHWAY_RULES["C_HIS"].evaluate(g_end, max_false_absence_prob=0.05)
            self.assertEqual(ev_end.state, STATE_UNRESOLVED)
            self.assertTrue(ev_end.end_truncated)

        # Unavailable feature order (opaque locus tags without GFF) -> retains unknown edge status
        opaque_rows = [
            {"protein_id": "LOCUS_ALPHA", "contig_id": "ctgA", "KO": "", "EC": ""},
            {"protein_id": "LOCUS_BETA", "contig_id": "ctgA", "KO": "K00765", "EC": ""},
            {"protein_id": "LOCUS_GAMMA", "contig_id": "ctgA", "KO": "", "EC": ""},
        ]
        g_unk, hits_unk = parse_per_mag_rows(
            opaque_rows, mag_id="dastool-unk-order", completeness=98.0, contig_lengths=[500_000]
        )
        self.assertTrue(all(h.edge_status == "unknown" for h in hits_unk))
        ev_unk = PATHWAY_RULES["C_HIS"].evaluate(g_unk, max_false_absence_prob=0.05)
        self.assertEqual(ev_unk.state, STATE_UNRESOLVED)
        self.assertTrue(ev_unk.unknown_edge)

        # GFF restores exact order for opaque locus tags
        gff_lines = [
            "##sequence-region ctgA 1 500000",
            "ctgA\tBakta\tCDS\t100\t900\t.\t+\t0\tID=LOCUS_ALPHA",
            "ctgA\tBakta\tCDS\t1100\t2100\t.\t+\t0\tID=LOCUS_BETA",
            "ctgA\tBakta\tCDS\t2300\t3200\t.\t+\t0\tID=LOCUS_GAMMA",
        ]
        gff_order, ctg_lens = parse_gff_cds_order(gff_lines)
        g_gff, _ = parse_per_mag_rows(
            list(reversed(opaque_rows)),
            mag_id="dastool-gff",
            completeness=98.0,
            contig_lengths=list(ctg_lens.values()),
            gff_cds_order=gff_order,
        )
        self.assertEqual(PATHWAY_RULES["C_HIS"].evaluate(g_gff, max_false_absence_prob=0.05).state, STATE_ABSENT)
        self.assertGreater(p_intact(4500, list(ctg_lens.values()), completeness=98.0), 0.95)

        # Blocker 4 (Astra 248968746164994048): zero hits WITH uncalibrated or poor detectability
        # MUST NOT produce a clean recipient gap in the entry point!
        g_uncal, _ = parse_per_mag_rows(
            [{"protein_id": "c1_00001", "contig_id": "c1", "KO": "K00548", "EC": ""}],
            mag_id="dastool-uncal",
            completeness=94.0,
            contig_lengths=(),  # No contig lengths and no empirical calibration!
        )
        d_full = GenomeAnnotations.from_kos("dastool-Dfull", AER_RING_KOS | BLUB_BAC_KOS)
        call_uncal, _ = evaluate_recipient_requirements_in_event(g_uncal, ["C_HIS"], [g_uncal, d_full])
        self.assertEqual(call_uncal, CALL_NOT_ELIGIBLE)

        # Fragmented short contigs (1.5 kb contigs) even at 98% completeness -> poor detectability -> NOT_ELIGIBLE!
        g_frag = GenomeAnnotations.from_kos(
            "dastool-frag", {"K00548"}, completeness=98.0, contig_lengths=[1500] * 2000
        )
        call_frag, _ = evaluate_recipient_requirements_in_event(g_frag, ["C_HIS"], [g_frag, d_full])
        self.assertEqual(call_frag, CALL_NOT_ELIGIBLE)

        # Empirical Null-1 calibration (< 0.05) or syntenic deletion unlocks eligibility
        g_emp = GenomeAnnotations.from_kos(
            "dastool-emp",
            {"K00548"},
            completeness=92.0,
            contig_lengths=[50_000] * 60,
            false_absence_probs={"C_HIS": 0.02},
        )
        d_his = GenomeAnnotations.from_kos("dastool-Dhis", {"K00765", "K01814", "K01693", "K00013"})
        call_emp, _ = evaluate_recipient_requirements_in_event(g_emp, ["C_HIS"], [g_emp, d_his])
        self.assertEqual(call_emp, CALL_D1_OR_CHAIN)

    def test_suite4_remodeling_wiring_3a_vs_3b_and_unbinned_guard(self) -> None:
        x_plain = GenomeAnnotations.from_kos("dastool-X1", {"K00548"})
        x_remod = GenomeAnnotations.from_kos("dastool-X2", {"K00548"} | B12_ATTACH_BAC_KOS)

        d1_ring = GenomeAnnotations.from_kos("dastool-D1", AER_RING_KOS)
        d2_dmb_only = GenomeAnnotations.from_kos("dastool-D2a", {"K04719"})
        d2_remod = GenomeAnnotations.from_kos("dastool-D2b", {"K04719"} | B12_ATTACH_BAC_KOS)
        d_full = GenomeAnnotations.from_kos("dastool-Dfull", AER_RING_KOS | BLUB_BAC_KOS)

        # 4a. Single complete B12 donor wins (d=1) alongside split donors
        call, wit = evaluate_recipient_requirements_in_event(
            x_remod, ["B12"], [x_remod, d1_ring, d2_dmb_only, d_full]
        )
        self.assertEqual(call, CALL_D1_OR_CHAIN)
        self.assertIsNotNone(wit)

        # 4b. Wiring 3a (Donor-side remodeling: D1 -> D2_remod -> X_plain) => NOT plexitrophic!
        call_3a, wit_3a = evaluate_recipient_requirements_in_event(
            x_plain, ["B12"], [x_plain, d1_ring, d2_remod]
        )
        self.assertEqual(call_3a, CALL_D1_OR_CHAIN)
        self.assertEqual(wit_3a, ("sink", "dastool-D2b"))

        # 4c. Wiring 3b (Recipient-side remodeling: D1 -> X_remod <- D2_dmb_only) => PLEXITROPHIC!
        call_3b, _ = evaluate_recipient_requirements_in_event(
            x_remod, ["B12"], [x_remod, d1_ring, d2_dmb_only]
        )
        self.assertEqual(call_3b, CALL_PLEXITROPHIC)

        # 4d. Both D2_remod and X_remod can remodel => 3a witness exists => NOT plexitrophic!
        call_both, _ = evaluate_recipient_requirements_in_event(
            x_remod, ["B12"], [x_remod, d1_ring, d2_remod]
        )
        self.assertEqual(call_both, CALL_D1_OR_CHAIN)

        # 4e. Unbinned KO pool carrying complete B12 downgrades d=2 call to unresolved_ambig_donor
        unbinned_b12 = GenomeAnnotations.from_kos("_unbinned", AER_RING_KOS | BLUB_BAC_KOS)
        call_unb, _ = evaluate_recipient_requirements_in_event(
            x_remod, ["B12"], [x_remod, d1_ring, d2_dmb_only], unbinned_genome=unbinned_b12
        )
        self.assertEqual(call_unb, CALL_UNRESOLVED_AMBIG)

        call_unb_only, _ = evaluate_recipient_requirements_in_event(
            x_plain, ["B12"], [x_plain], unbinned_genome=unbinned_b12
        )
        self.assertEqual(call_unb_only, CALL_UNRESOLVED_AMBIG)

    def test_suite5_entry_point_cycle_long_chains_and_non_laundering(self) -> None:
        # Blocker 1 (Astra 248968746164994048): Test an irreducible 3-cycle M_A -> M_B -> M_C -> M_A
        # DIRECTLY through `evaluate_recipient_requirements_in_event(..., mode='consortium_cyc')`!
        #   M_A: needs B12 (metH\metE); has ASSIM_SO4; has A_NIF (exports FIXED_N)
        #   M_B: needs FIXED_N; has metE (no B12 need, no B12 export); has A_DSR (exports RED_S)
        #   M_C: needs RED_S; has ASSIM_NO3 (no FIXED_N need, no A_NIF export); has aerobic B12 (exports B12)
        so4_kos = {"K00958", "K05907", "K00380", "K00381"}
        no3_kos = {"K00367", "K00366"}
        nif_kos = {"K02588", "K02586", "K02591"}
        dsr_kos = {"K00394", "K00395", "K11180", "K11181"}

        m_a = GenomeAnnotations.from_kos("dastool-MA", {"K00548"} | so4_kos | nif_kos)
        m_b = GenomeAnnotations.from_kos("dastool-MB", {"K00549"} | so4_kos | dsr_kos)
        m_c = GenomeAnnotations.from_kos("dastool-MC", {"K00549"} | no3_kos | AER_RING_KOS | BLUB_BAC_KOS)

        call_cyc, _ = evaluate_recipient_requirements_in_event(
            m_a,
            ["B12"],
            [m_a, m_b, m_c],
            scenario=SCENARIO_E_MINIMAL,
            mode="consortium_cyc",
            donor_need_endpoints=("B12", "FIXED_N", "RED_S"),
        )
        self.assertEqual(call_cyc, CALL_PLEXITROPHIC)  # M_A -> M_B -> M_C -> M_A (d=1, kappa=1)!

        # Adding ambient FIXED_N (SCENARIO_E_PLUS_DIN) breaks M_A -> M_B, turning the 3-cycle into
        # the serial chain M_B -> M_C -> M_A (d=1, kappa=0 => CALL_D1_OR_CHAIN)!
        call_broken, _ = evaluate_recipient_requirements_in_event(
            m_a,
            ["B12"],
            [m_a, m_b, m_c],
            scenario=SCENARIO_E_PLUS_DIN,
            mode="consortium_cyc",
            donor_need_endpoints=("B12", "FIXED_N", "RED_S"),
        )
        self.assertEqual(call_broken, CALL_D1_OR_CHAIN)

        # Blocker 2 (Astra 248968746164994048): 6-member serial chain C1->C2->C3->C4->C5->X (m=6 > 4)
        # coexisting with a size-3 fan-in {F1, F2, X} MUST return 'not_plexi' (never truncated)!
        long_chain_profiles = {
            "X": OrganismProfile(
                "X",
                conversions=((frozenset({"A", "B"}), "TARGET"),),
                required=frozenset({"TARGET"}),
            ),
            "F1": OrganismProfile("F1", de_novo=frozenset({"A"})),
            "F2": OrganismProfile("F2", de_novo=frozenset({"B"})),
            "C1": OrganismProfile("C1", de_novo=frozenset({"S1"})),
            "C2": OrganismProfile("C2", conversions=((frozenset({"S1"}), "S2"),)),
            "C3": OrganismProfile("C3", conversions=((frozenset({"S2"}), "S3"),)),
            "C4": OrganismProfile("C4", conversions=((frozenset({"S3"}), "S4"),)),
            "C5": OrganismProfile("C5", conversions=((frozenset({"S4"}), "TARGET"),)),
        }
        st_long, wit_long = decide_nonplexi_p_all("X", long_chain_profiles, cyc_mode="cyc")
        self.assertEqual(st_long, "not_plexi")
        self.assertEqual(wit_long, ("sink", "C5"))

        # Astra 248969335645061120 blocking disagreement fixture on plexisolve.py:
        # A makes/exports p. D makes s unconditionally, has p->q and q->p, exports p,s. X needs p,s.
        # D CANNOT import p, run p->q->p, and re-export p as a single-supplier witness!
        launder_1 = {
            "A": OrganismProfile("A", de_novo=frozenset({"p"}), exportable=frozenset({"p"})),
            "D": OrganismProfile(
                "D",
                de_novo=frozenset({"s"}),
                conversions=((frozenset({"p"}), "q"), (frozenset({"q"}), "p")),
                exportable=frozenset({"p", "s"}),
            ),
            "X": OrganismProfile("X", required=frozenset({"p", "s"}), exportable=frozenset()),
        }
        for mode in ("fanin_only", "cyc", "scc"):
            st_l1, _ = decide_nonplexi_p_all("X", launder_1, cyc_mode=mode)
            self.assertEqual(st_l1, "plexi", f"Failed 1-metabolite laundering guard in {mode}")

        # 2-metabolite mutual laundering generalization:
        # A1 exports p1, A2 exports p2; D has p1->p2 and p2->p1, exports p1,p2; X needs p1,p2.
        # Individually D can produce p2 (from p1) or p1 (from p2), NEVER both in one wiring!
        launder_2 = {
            "A1": OrganismProfile("A1", de_novo=frozenset({"p1"})),
            "A2": OrganismProfile("A2", de_novo=frozenset({"p2"})),
            "D": OrganismProfile(
                "D",
                conversions=((frozenset({"p1"}), "p2"), (frozenset({"p2"}), "p1")),
                exportable=frozenset({"p1", "p2"}),
            ),
            "X": OrganismProfile("X", required=frozenset({"p1", "p2"}), exportable=frozenset()),
        }
        for mode in ("fanin_only", "cyc", "scc"):
            st_l2, _ = decide_nonplexi_p_all("X", launder_2, cyc_mode=mode)
            self.assertEqual(st_l2, "plexi", f"Failed 2-metabolite laundering guard in {mode}")

        # Astra 248973173642936320 4-organism cross-organism laundering counterexample:
        # A makes p; D makes s and has p->q and r->p; B has q->r; X needs p,s.
        # D CANNOT import p from A to start D:p->D:q->B:r->D:p and claim rule-produced p!
        launder_4org = {
            "A": OrganismProfile("A", de_novo=frozenset({"p"})),
            "D": OrganismProfile(
                "D",
                de_novo=frozenset({"s"}),
                conversions=((frozenset({"p"}), "q"), (frozenset({"r"}), "p")),
            ),
            "B": OrganismProfile("B", conversions=((frozenset({"q"}), "r"),)),
            "X": OrganismProfile("X", required=frozenset({"p", "s"}), exportable=frozenset()),
        }
        for mode in ("fanin_only", "cyc", "scc"):
            st_l4, _ = decide_nonplexi_p_all("X", launder_4org, cyc_mode=mode)
            self.assertEqual(st_l4, "plexi", f"Failed 4-org cross-organism laundering guard in {mode}")

        # 5-organism 'cyc' private-attachment cross-export laundering counterexample:
        # X makes u (in 'cyc' X may export u only to D, NOT to C); D has u->q and r->u, plus w->k;
        # B has q->r; C needs u and makes w; X needs k. Also has fallback fan-in {F1(k1), F2(k2)} -> X(k).
        # D cannot both import u from X AND rule-produce u to export u to C!
        launder_5org_cyc = {
            "X": OrganismProfile(
                "X",
                de_novo=frozenset({"u"}),
                conversions=((frozenset({"k1", "k2"}), "k"),),
                required=frozenset({"k"}),
            ),
            "D": OrganismProfile(
                "D",
                conversions=(
                    (frozenset({"u"}), "q"),
                    (frozenset({"r"}), "u"),
                    (frozenset({"w"}), "k"),
                ),
            ),
            "B": OrganismProfile("B", conversions=((frozenset({"q"}), "r"),)),
            "C": OrganismProfile("C", de_novo=frozenset({"w"}), required=frozenset({"u"})),
            "F1": OrganismProfile("F1", de_novo=frozenset({"k1"})),
            "F2": OrganismProfile("F2", de_novo=frozenset({"k2"})),
        }
        self.assertFalse(is_metabolite_dep_acyclic(launder_5org_cyc))
        self.assertEqual(decide_nonplexi_p_all("X", launder_5org_cyc, cyc_mode="cyc")[0], "plexi")

        # Verify is_metabolite_dep_acyclic precondition on acyclic vs cyclic fixtures
        self.assertTrue(is_metabolite_dep_acyclic(long_chain_profiles))
        self.assertTrue(
            is_metabolite_dep_acyclic(
                {g.mag_id: build_mag_organism_profile(g) for g in (m_a, m_b, m_c)}
            )
        )
        self.assertFalse(is_metabolite_dep_acyclic(launder_1))
        self.assertFalse(is_metabolite_dep_acyclic(launder_2))
        self.assertFalse(is_metabolite_dep_acyclic(launder_4org))

        # F8 with r->p removed is acyclic and still 'plexi' via the pure-stdlib Theorem 3 path
        launder_4org_acyclic = dict(launder_4org)
        launder_4org_acyclic["D"] = OrganismProfile(
            "D", de_novo=frozenset({"s"}), conversions=((frozenset({"p"}), "q"),)
        )
        self.assertTrue(is_metabolite_dep_acyclic(launder_4org_acyclic))
        self.assertEqual(decide_nonplexi_p_all("X", launder_4org_acyclic, cyc_mode="cyc")[0], "plexi")

        # Verify enumerate_minimal_wirings / is_wiring_feasible also blocks F8 laundering
        w_f8 = enumerate_minimal_wirings("X", launder_4org)
        self.assertEqual([w.invariants("X")[0] for w in w_f8], [2])

        # Theorem 3' F6 interconverter SAT/UNSAT gadget (Opus 248974707227279360)
        def _make_f6(clauses: list[tuple[int, ...]], nvar: int) -> dict[str, OrganismProfile]:
            x_convs = tuple(
                (frozenset({f"v{abs(lit)}{'+' if lit > 0 else '-'}"}), f"c{j}")
                for j, cl in enumerate(clauses)
                for lit in cl
            )
            orgs: dict[str, OrganismProfile] = {
                "X": OrganismProfile(
                    "X",
                    conversions=x_convs,
                    required=frozenset({"w"} | {f"c{j}" for j in range(len(clauses))}),
                    exportable=frozenset(),
                )
            }
            d_convs = []
            for i in range(1, nvar + 1):
                d_convs.append((frozenset({f"v{i}-"}), f"v{i}+"))
                d_convs.append((frozenset({f"v{i}+"}), f"v{i}-"))
                orgs[f"A{i}+"] = OrganismProfile(f"A{i}+", de_novo=frozenset({f"v{i}+"}))
                orgs[f"A{i}-"] = OrganismProfile(f"A{i}-", de_novo=frozenset({f"v{i}-"}))
            orgs["D"] = OrganismProfile(
                "D",
                de_novo=frozenset({"w"}),
                conversions=tuple(d_convs),
            )
            return orgs

        f6_sat = _make_f6([(1, 2), (-1, 2), (1, -2)], 2)
        f6_unsat = _make_f6([(1, 2), (-1, 2), (1, -2), (-1, -2)], 2)
        for mode in ("fanin_only", "cyc", "scc"):
            self.assertEqual(decide_nonplexi_p_all("X", f6_sat, cyc_mode=mode), ("not_plexi", ("sink", "D")))
            self.assertEqual(decide_nonplexi_p_all("X", f6_unsat, cyc_mode=mode)[0], "plexi")

        # Adversarial Case A3 (Opus 248969156292427776): X <-> D <-> Y (two dyads sharing D)
        # |SCC(X)| = 3 ('scc' calls plexi), whereas simple-cycle kappa(X) = 0 ('cyc' calls not_plexi)!
        dyad_chain = {
            "X": OrganismProfile("X", de_novo=frozenset({"mx"}), required=frozenset({"md"})),
            "D": OrganismProfile(
                "D", de_novo=frozenset({"md", "mdy"}), required=frozenset({"mx", "my"})
            ),
            "Y": OrganismProfile("Y", de_novo=frozenset({"my"}), required=frozenset({"mdy"})),
        }
        self.assertEqual(decide_nonplexi_p_all("X", dyad_chain, cyc_mode="cyc")[0], "not_plexi")
        self.assertEqual(decide_nonplexi_p_all("X", dyad_chain, cyc_mode="scc")[0], "plexi")

        # Commensal donor D -> A feeding 3-cycle A -> B -> C -> A via enumerate_minimal_wirings
        cyc3_profiles = {
            "D": OrganismProfile("D", de_novo=frozenset({"M0"})),
            "A": OrganismProfile("A", de_novo=frozenset({"M_AB"}), required=frozenset({"M0", "M_CA"})),
            "B": OrganismProfile("B", de_novo=frozenset({"M_BC"}), required=frozenset({"M_AB"})),
            "C": OrganismProfile("C", de_novo=frozenset({"M_CA"}), required=frozenset({"M_BC"})),
        }
        w_a = enumerate_minimal_wirings("A", cyc3_profiles)
        self.assertEqual(len(w_a), 1)
        self.assertEqual(w_a[0].invariants("A"), (2, 4, 3))
        self.assertEqual(w_a[0].invariants("B"), (1, 4, 3))
        self.assertEqual(w_a[0].kappa("B"), 1)
        self.assertEqual(decide_nonplexi_p_all("B", cyc3_profiles, cyc_mode="cyc")[0], "plexi")

    def test_suite6_downsampling_position_and_boundary_invariance(self) -> None:
        # Blocker 3 (Astra 248968746164994048):
        # 1) Retention must span early, middle, and late positions (NOT prefix-only!)
        # 2) Two original contigs must NEVER be merged into one simulated contig!
        # 3) Reported completeness must reflect ACTUAL achieved gene retention!
        rows = []
        # Two original contigs: orig1 (10 genes, K00001..K00010), orig2 (10 genes, K00011..K00020)
        for i in range(1, 11):
            rows.append({"protein_id": f"orig1_{i:05d}", "contig_id": "orig1", "KO": f"K000{i:02d}", "EC": ""})
        for i in range(1, 11):
            rows.append({"protein_id": f"orig2_{i:05d}", "contig_id": "orig2", "KO": f"K000{i + 10:02d}", "EC": ""})

        _, full_hits = parse_per_mag_rows(rows, mag_id="dastool-ref20", complete_cds_table=True)
        orig1_kos = {f"K000{i:02d}" for i in range(1, 11)}
        orig2_kos = {f"K000{i:02d}" for i in range(11, 21)}

        rng = random.Random(20261119)
        first_kept = 0
        mid_kept = 0
        last_kept = 0
        n_trials = 200

        for t in range(n_trials):
            ds = downsample_contig_blocks(
                full_hits,
                target_completeness=50.0,
                target_contig_count=10,
                rng=rng,
                mag_id=f"ds_{t}",
                orig_completeness=96.0,
                total_genome_bp=200_000,
            )
            if "K00001" in ds.kos:
                first_kept += 1
            if "K00010" in ds.kos:
                mid_kept += 1
            if "K00020" in ds.kos:
                last_kept += 1
            # Achieved completeness must match 96.0 * (len(ds.kos) / 20.0)
            expected_comp = round(96.0 * (len(ds.kos) / 20.0), 4)
            self.assertAlmostEqual(ds.completeness, expected_comp, places=3)

        # Each position should be retained roughly ~50% of the time (well inside [0.32, 0.68])
        for count in (first_kept, mid_kept, last_kept):
            rate = count / n_trials
            self.assertGreater(rate, 0.32)
            self.assertLess(rate, 0.68)

        # Verify directly that target_contig_count=1 (fewer than 2 original contigs!) STILL keeps >=2 contigs
        # and NEVER merges orig1 and orig2 onto one simulated contig!
        ds_nomerge = downsample_contig_blocks(
            full_hits,
            target_completeness=100.0,
            target_contig_count=1,
            rng=rng,
            orig_completeness=100.0,
        )
        self.assertEqual(len(ds_nomerge.contig_lengths), 2)
        self.assertTrue(orig1_kos.issubset(ds_nomerge.kos) and orig2_kos.issubset(ds_nomerge.kos))

    def test_suite7_trap9_scenarios_bounds_and_verdicts(self) -> None:
        with self.assertRaises(ValueError):
            assert_non_redundant_mag_ids(["dastool-semibin_001", "semibin_001"])

        q_rows = [
            {"Name": "dastool-b1", "Completeness": "94.0", "Contamination": "1.2"},
            {"Name": "semibin_001", "Completeness": "94.0", "Contamination": "1.2"},
            {"Name": "dastool-euk1", "Completeness": "96.0", "Contamination": "0.5"},
            {"Name": "dastool-arch1", "Completeness": "92.0", "Contamination": "2.0"},
        ]
        tax = {
            "dastool-b1": "d__Bacteria;p__Pseudomonadota",
            "dastool-euk1": "d__Eukaryota;p__Chlorophyta",
            "dastool-arch1": "d__Archaea;p__Thermoproteota",
        }
        kept, counts = filter_derep_prokaryotic_mags(
            q_rows,
            derep_representatives={"dastool-b1", "dastool-euk1", "dastool-arch1"},
            taxonomy_by_mag=tax,
            min_completeness=90.0,
            max_contamination=5.0,
        )
        self.assertEqual(kept, ["dastool-arch1", "dastool-b1"])
        self.assertEqual(counts["prokaryotic_bins"], 2)

        x_ri2 = GenomeAnnotations.from_kos("dastool-Xri2", {"K00548"})
        d_b12 = GenomeAnnotations.from_kos("dastool-Db12", AER_RING_KOS | BLUB_BAC_KOS)
        d_nif = GenomeAnnotations.from_kos("dastool-Dnif", {"K02588", "K02586", "K02591"})

        self.assertEqual(
            evaluate_recipient_requirements_in_event(
                x_ri2, ["B12", "FIXED_N"], [x_ri2, d_b12, d_nif], SCENARIO_E_MINIMAL
            )[0],
            CALL_PLEXITROPHIC,
        )
        self.assertEqual(
            evaluate_recipient_requirements_in_event(
                x_ri2, ["B12", "FIXED_N"], [x_ri2, d_b12, d_nif], SCENARIO_E_PLUS_DIN
            )[0],
            CALL_D1_OR_CHAIN,
        )
        self.assertEqual(
            evaluate_recipient_requirements_in_event(
                x_ri2,
                ["B12", "FIXED_N"],
                [x_ri2, d_b12, d_nif],
                SupplyScenario("E_both", frozenset({"B12", "FIXED_N"})),
            )[0],
            CALL_D0_AMBIENT,
        )

        self.assertEqual(classify_coencoder_bypass(GenomeAnnotations.from_kos("c", {"K25172"})), "S_spatial")
        self.assertEqual(classify_coencoder_bypass(GenomeAnnotations.from_kos("a", {"K00425", "K00426"})), "S_temp_resp")
        self.assertEqual(classify_coencoder_bypass(GenomeAnnotations.from_kos("x", set(), contamination=8.5)), "chimeric_candidate")

        rng = random.Random(123)
        mat = [[1, 1, 0, 0, 1], [0, 1, 1, 0, 0], [1, 0, 1, 1, 0], [0, 0, 0, 1, 1]]
        perm_mat = curveball_permute_binary_matrix(mat, rng, n_steps=40)
        self.assertEqual([sum(r) for r in perm_mat], [sum(r) for r in mat])
        self.assertEqual([sum(perm_mat[i][j] for i in range(4)) for j in range(5)], [sum(mat[i][j] for i in range(4)) for j in range(5)])

        obs = StratumCounts(n_auto=5, n_d0=2, n_d1=10, n_plexi=6, n_ambig=2, n_no_donor=0)
        nul = StratumCounts(n_auto=5, n_d0=2, n_d1=15, n_plexi=2, n_ambig=1, n_no_donor=0)
        lo_id, hi_id = observed_minus_null_id_bounds(obs, nul)
        self.assertAlmostEqual(lo_id, 0.15)
        self.assertAlmostEqual(hi_id, 0.30)
        self.assertLess(*manski_log_odds_bounds(obs, nul))

        gate_msg = check_data_adequacy_gates(4, 0.1, 0.1)
        self.assertIsNotNone(gate_msg)
        self.assertEqual(decide_p2_q1((-0.02, 0.03), (0.01, 0.09), delta_p=0.05), "Falsified")
        self.assertEqual(decide_p2_q1((0.06, 0.14), (0.07, 0.18), delta_p=0.05), "Supported")
        self.assertEqual(decide_p2_q2((0.80, 1.40)), "Supported")

    def test_suite8_graded_biological_state_posteriors(self) -> None:
        # Smooth pathway posterior increases monotonically with step count
        p0 = compute_pathway_posterior(PATHWAY_RULES["C_HIS"], GenomeAnnotations.from_kos("g0", set()))
        p2 = compute_pathway_posterior(
            PATHWAY_RULES["C_HIS"], GenomeAnnotations.from_kos("g2", {"K00765", "K01814"})
        )
        p4 = compute_pathway_posterior(
            PATHWAY_RULES["C_HIS"],
            GenomeAnnotations.from_kos("g4", {"K00765", "K01814", "K01693", "K00013"}),
        )
        self.assertLess(p0, 0.05)
        self.assertLess(p0, p2)
        self.assertLess(p2, p4)
        self.assertGreater(p4, 0.99)

        # Exact Proposition 1 check: N=10 i.i.d. donors, q=0.3, k=2 -> P(fanin>=2) = 0.333719!
        donors_10 = [(1.0, {"R1": 0.3, "R2": 0.3}) for _ in range(10)]
        post = compute_fanin_state_posterior({"R1": 1.0, "R2": 1.0}, donors_10)
        self.assertAlmostEqual(post.p_fanin_ge_2, 0.333719, places=5)
        self.assertAlmostEqual(post.p_d1_single_supplier + post.p_fanin_ge_2 + post.p_uncovered, 1.0, places=5)

        # Invariance to duplicate internal routes (Astra 248970428366110720):
        # One certain donor for R1 + one certain donor for R2 + zero dual-carrier -> P(fanin>=2) == 1.0;
        # adding a second certain R1-only donor creates a duplicate 2-partner route {D1b, D2},
        # yet the biological-state probability P(no single partner suffices) stays 1.0!
        post_base = compute_fanin_state_posterior(
            {"R1": 1.0, "R2": 1.0}, [(1.0, {"R1": 1.0, "R2": 0.0}), (1.0, {"R1": 0.0, "R2": 1.0})]
        )
        post_dup = compute_fanin_state_posterior(
            {"R1": 1.0, "R2": 1.0},
            [
                (1.0, {"R1": 1.0, "R2": 0.0}),
                (1.0, {"R1": 1.0, "R2": 0.0}),
                (1.0, {"R1": 0.0, "R2": 1.0}),
            ],
        )
        self.assertAlmostEqual(post_base.p_fanin_ge_2, 1.0, places=6)
        self.assertAlmostEqual(post_dup.p_fanin_ge_2, 1.0, places=6)

        # Coupled within-donor joint PMFs (Astra 248973173642936320):
        # Mutual exclusion p11=0 with marginals q1=0.3, q2=0.3 across N=10 donors:
        # Exact P(fanin>=2) = 1 - 2*(0.7^10) + 0.4^10 = 0.943610 (matches Table B of fanin_closed_form.out)
        coupled_pmfs = [
            {
                frozenset(): 0.4,
                frozenset({"R1"}): 0.3,
                frozenset({"R2"}): 0.3,
                frozenset({"R1", "R2"}): 0.0,
            }
            for _ in range(10)
        ]
        post_coupled = compute_fanin_state_posterior(
            {"R1": 1.0, "R2": 1.0}, donor_joint_pmfs=coupled_pmfs
        )
        self.assertAlmostEqual(post_coupled.p_fanin_ge_2, 0.943610, places=5)
        self.assertAlmostEqual(post_coupled.p_d1_single_supplier, 0.0, places=6)

    def test_suite9_random_acyclic_and_cyclic_s1_oracle_fuzz(self) -> None:
        """Cross-check decide_nonplexi_p_all against an independent single-node S1 brute oracle."""

        def _brute_s1_oracle(
            focal_id: str,
            profiles: dict[str, OrganismProfile],
            amb: frozenset[str],
            cyc_mode: str,
        ) -> str:
            orgs = sorted(profiles.keys())
            uses: dict[str, set[str]] = {}
            prods: dict[str, set[str]] = {}
            for y in orgs:
                u = set(profiles[y].required)
                for ins, _ in profiles[y].rules():
                    u |= set(ins)
                uses[y] = u - set(amb)
                prods[y] = {p for _, p in profiles[y].rules()} - set(amb)

            cand_edges = sorted(
                {
                    (z, y)
                    for z in orgs
                    for y in orgs
                    if z != y
                    and any(
                        profiles[z].can_export(m)
                        for m in (prods[z] & uses[y])
                    )
                }
            )
            forced_r: set[tuple[str, str]] = set()
            choice_r: list[tuple[str, str]] = []
            for z in orgs:
                for m in prods[z]:
                    importable = m in uses[z] and any(
                        w != z and m in prods[w] and profiles[w].can_export(m) for w in orgs
                    )
                    if importable:
                        choice_r.append((z, m))
                    else:
                        forced_r.add((z, m))

            def _dag_smax(allowed_edges: set[tuple[str, str]], r_nodes: set[tuple[str, str]]) -> set[str]:
                s_cur = set(orgs)
                while True:
                    avail = {y: set(amb) for y in s_cur}
                    changed = True
                    while changed:
                        changed = False
                        for y in s_cur:
                            for ins, p in profiles[y].rules():
                                if (y, p) in r_nodes and p not in avail[y] and ins.issubset(avail[y]):
                                    avail[y].add(p)
                                    changed = True
                            for m in uses[y] - avail[y]:
                                if (y, m) in r_nodes:
                                    continue
                                for z in s_cur:
                                    if (
                                        z != y
                                        and (z, y) in allowed_edges
                                        and (z, m) in r_nodes
                                        and m in avail[z]
                                        and profiles[z].can_export(m)
                                    ):
                                        avail[y].add(m)
                                        changed = True
                                        break
                    dead = {y for y in s_cur if not profiles[y].required.issubset(avail[y])}
                    if not dead:
                        return s_cur
                    s_cur -= dead

            all_edges = set(cand_edges)
            viable_any = any(
                focal_id in _dag_smax(all_edges, forced_r | set(rc))
                for k_r in range(len(choice_r) + 1)
                for rc in itertools.combinations(choice_r, k_r)
            )
            if not viable_any:
                return "no_donor"
            if focal_id in _dag_smax(set(), forced_r | set(choice_r)):
                return "d0"

            for k_e in range(len(cand_edges) + 1):
                for e_sub in itertools.combinations(cand_edges, k_e):
                    e_set = set(e_sub)
                    if len({z for (z, y) in e_set if y == focal_id}) > 1:
                        continue
                    if cyc_mode == "scc":
                        sccs = tarjan_sccs(orgs, e_set)
                        if next((len(c) for c in sccs if focal_id in c), 1) > 2:
                            continue
                    elif cyc_mode == "cyc":
                        if on_simple_cycle_ge_3(orgs, e_set, focal_id):
                            continue
                    for k_r in range(len(choice_r) + 1):
                        for rc in itertools.combinations(choice_r, k_r):
                            if focal_id in _dag_smax(e_set, forced_r | set(rc)):
                                return "not_plexi"
            return "plexi"

        rng = random.Random(20261119)
        mets = ["m0", "m1", "m2", "m3", "m4"]
        n_checked = 0
        n_cyclic = 0
        while n_checked < 80:
            amb = frozenset(m for m in mets if rng.random() < 0.1)
            profs: dict[str, OrganismProfile] = {}
            for idx_o in range(4):
                oid = "X" if idx_o == 0 else f"O{idx_o}"
                de_novo = set()
                convs = []
                for m in mets:
                    r_draw = rng.random()
                    if r_draw < 0.22:
                        de_novo.add(m)
                    elif r_draw < 0.52:
                        ins = frozenset(rng.sample([x for x in mets if x != m], rng.choice([1, 2])))
                        convs.append((ins, m))
                prods = de_novo | {p for _, p in convs}
                exp = frozenset(p for p in prods if rng.random() < 0.85)
                reqs = frozenset(
                    rng.sample(mets, rng.choice([1, 2]) if idx_o == 0 else rng.choice([0, 1, 2]))
                )
                profs[oid] = OrganismProfile(oid, frozenset(de_novo), tuple(convs), reqs, exp)

            # Bound candidate edge count so exact brute-force oracle stays sub-second
            cand_cnt = sum(
                1
                for z in profs
                for y in profs
                if z != y
                and any(
                    profs[z].can_export(m)
                    for m in (
                        (profs[z].de_novo | {p for _, p in profs[z].conversions})
                        & (profs[y].required | {a for ins, _ in profs[y].conversions for a in ins})
                        - amb
                    )
                )
            )
            if cand_cnt > 9:
                continue
            n_checked += 1
            if not is_metabolite_dep_acyclic(profs, amb):
                n_cyclic += 1
            for mode in ("fanin_only", "cyc", "scc"):
                got, _ = decide_nonplexi_p_all("X", profs, amb, cyc_mode=mode)
                want = _brute_s1_oracle("X", profs, amb, mode)
                self.assertEqual(got, want, f"Mismatch on trial {n_checked} mode={mode}")
        self.assertGreater(n_cyclic, 15, "Fuzz suite must exercise cyclic rule graphs too")

    def test_suite10_eukaryotic_donor_class(self) -> None:
        """Verify organism-resolved, unknown-detectability, and biologically excluded eukaryotic donor rules."""
        # Prokaryotic recipient needing THI + C_HIS; bacterial donors split (D_thi, D_his -> d=2 alone)
        x_prok = GenomeAnnotations.from_kos("dastool-X", set())
        d_thi = GenomeAnnotations.from_kos(
            "dastool-Dthi", {"K03147", "K00941", "K03149", "K03153", "K00788", "K00946"}
        )
        d_his = GenomeAnnotations.from_kos("dastool-Dhis", {"K00765", "K01814", "K01693", "K00013"})
        d_nif = GenomeAnnotations.from_kos("dastool-Dnif", {"K02588", "K02586", "K02591"})

        # Astra 248978604545196032 reproducer: empty Eukaryota KO record retains eukaryote-eligible C_HIS
        # in I^+ (STATE_UNRESOLVED) while excluding prokaryote-only A_NIF and P_B12_AER (STATE_ABSENT)!
        euk_empty = GenomeAnnotations.from_kos("euk-empty", set(), domain="Eukaryota")
        self.assertEqual(
            PATHWAY_RULES["C_HIS"].evaluate(euk_empty, max_false_absence_prob=0.05).state,
            STATE_UNRESOLVED,
        )
        self.assertEqual(
            PATHWAY_RULES["A_NIF"].evaluate(euk_empty, max_false_absence_prob=0.05).state,
            STATE_ABSENT,
        )
        self.assertEqual(
            PATHWAY_RULES["P_B12_AER"].evaluate(euk_empty, max_false_absence_prob=0.05).state,
            STATE_ABSENT,
        )
        prof_empty_plus = build_mag_organism_profile(
            euk_empty, optimistic_capabilities=True, max_false_absence_prob=0.05
        )
        self.assertIn("C_HIS", prof_empty_plus.de_novo)
        self.assertNotIn("FIXED_N", prof_empty_plus.de_novo)
        self.assertNotIn("B12", prof_empty_plus.de_novo)

        # Eukaryotic genomes are rejected as focal recipients in Pass 1
        euk_rec = GenomeAnnotations.from_kos("euk-alga1", set(), domain="Eukaryota")
        self.assertEqual(
            evaluate_recipient_requirements_in_event(euk_rec, ["THI"], [d_thi])[0],
            CALL_NOT_ELIGIBLE,
        )

        # 1. Organism-resolved alga encoding eukaryotic THI (THIC + THI6 + THI4 + TPK1) + C_HIS -> d=1!
        alga_full = GenomeAnnotations.from_kos(
            "euk-alga-resolved",
            {"K03147", "K14153", "K03146", "K00949", "K00765", "K01814", "K01693", "K00013"},
            domain="Eukaryota",
            completeness=45.0,
            contig_lengths=(),
        )
        call_res, wit_res = evaluate_recipient_requirements_in_event(
            x_prok,
            ["THI", "C_HIS"],
            [x_prok, d_thi, d_his],
            eukaryotic_donors=[alga_full],
        )
        self.assertEqual(call_res, CALL_D1_OR_CHAIN)
        self.assertEqual(wit_res, ("sink", "euk-alga-resolved"))

        # 2. Astra 248978604545196032 mandatory fixture: one missing eukaryote-eligible route (C_HIS)
        # in an uncalibrated eukaryotic bin changes robust bacterial fan-in (CALL_PLEXITROPHIC)
        # to CALL_UNRESOLVED_AMBIG (retained in I^+), UNLESS calibrated or annotated-only mode!
        alga_missing_his = GenomeAnnotations.from_kos(
            "euk-alga-frag",
            {"K03147", "K14153", "K03146", "K00949"},  # THI present, 0 hits on C_HIS
            domain="Eukaryota",
            completeness=45.0,
            contig_lengths=(),
        )
        self.assertEqual(
            evaluate_recipient_requirements_in_event(x_prok, ["THI", "C_HIS"], [x_prok, d_thi, d_his])[0],
            CALL_PLEXITROPHIC,
        )
        call_frag_iplus, wit_frag = evaluate_recipient_requirements_in_event(
            x_prok,
            ["THI", "C_HIS"],
            [x_prok, d_thi, d_his],
            eukaryotic_donors=[alga_missing_his],
        )
        self.assertEqual(call_frag_iplus, CALL_UNRESOLVED_AMBIG)
        self.assertEqual(wit_frag, ("sink", "euk-alga-frag"))

        # Explicitly calibrated eukaryotic absence (false_absence_probs < 0.05) OR annotated-only bound keeps plexi
        alga_cal_no_his = GenomeAnnotations.from_kos(
            "euk-alga-cal",
            {"K03147", "K14153", "K03146", "K00949"},
            domain="Eukaryota",
            false_absence_probs={"C_HIS": 0.01},
        )
        self.assertEqual(
            evaluate_recipient_requirements_in_event(
                x_prok, ["THI", "C_HIS"], [x_prok, d_thi, d_his], eukaryotic_donors=[alga_cal_no_his]
            )[0],
            CALL_PLEXITROPHIC,
        )
        self.assertEqual(
            evaluate_recipient_requirements_in_event(
                x_prok,
                ["THI", "C_HIS"],
                [x_prok, d_thi, d_his],
                eukaryotic_donors=[alga_missing_his],
                euk_annotated_routes_only=True,
            )[0],
            CALL_PLEXITROPHIC,
        )

        # Biologically excluded prokaryote-only route (A_NIF -> FIXED_N) on the SAME uncalibrated alga
        # is NOT retained in I^+, so ["FIXED_N", "THI"] remains robust CALL_PLEXITROPHIC!
        self.assertEqual(
            evaluate_recipient_requirements_in_event(
                x_prok,
                ["FIXED_N", "THI"],
                [x_prok, d_nif, d_thi],
                eukaryotic_donors=[alga_missing_his],
            )[0],
            CALL_PLEXITROPHIC,
        )

        # 3. Release-gated euk_exportable_endpoints: if C_HIS is marked non-released by euk, stays plexi
        call_gated, _ = evaluate_recipient_requirements_in_event(
            x_prok,
            ["THI", "C_HIS"],
            [x_prok, d_thi, d_his],
            eukaryotic_donors=[alga_full],
            euk_exportable_endpoints=frozenset({"THI"}),
        )
        self.assertEqual(call_gated, CALL_PLEXITROPHIC)

        # 4. Unassigned multi-taxon pooled MetaEuk hits (with annotated-only KOs) go STRICTLY to I^+ -> unresolved_ambig
        call_pooled, _ = evaluate_recipient_requirements_in_event(
            x_prok,
            ["THI", "C_HIS"],
            [x_prok, d_thi, d_his],
            unassigned_euk_genome=alga_full,
        )
        self.assertEqual(call_pooled, CALL_UNRESOLVED_AMBIG)

        # 5. Matched-event transition audit
        audit = compare_calls_with_eukaryotic_donors(
            [
                (x_prok, ["THI", "C_HIS"], [x_prok, d_thi, d_his], [alga_full], None),
                (x_prok, ["THI", "C_HIS"], [x_prok, d_thi, d_his], [alga_missing_his], None),
            ]
        )
        self.assertEqual(audit["n_evaluated"], 2)
        self.assertEqual(audit["n_changed"], 2)
        self.assertEqual(
            audit["transitions"],
            {
                (CALL_PLEXITROPHIC, CALL_D1_OR_CHAIN): 1,
                (CALL_PLEXITROPHIC, CALL_UNRESOLVED_AMBIG): 1,
            },
        )


if __name__ == "__main__":
    unittest.main(verbosity=2)
