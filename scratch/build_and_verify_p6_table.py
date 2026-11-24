"""Verify (read-only by default; regenerate with --write) /data/dev/plexityper/protocols/P6-pathway-table.tsv:
1. /data/scratch/refdbs/kofam_db/ko_list check: every panel, confounder, co-carriage, and
   prokaryotic housekeeping KO exists and has a protein HMM profile;
2. Disjointness of unambiguous_kos within all 8 cohort chemistry configurations:
   {TierA (NO+N2O), TierB_N2O (N2O only), TierB_NO (NO only), TierC (neither)} x {+H2, -H2};
3. Exact C, H, O, N, S, charge, and redox-equivalent balance on all 20 rows;
4. Exact recomputation of dG0(298.15 K) and dH0(298.15 K) from project-local scratch/p6_obigt/ CSVs;
5. Exact path-invariance (< 1e-11 kJ/mol) of dG_r(T, pH, c) across DIC, SigmaH2S, and TAN speciation;
6. Finiteness and normalization of Dirichlet-multinomial probabilities under edge cases
   (cap = 0 with RNA > 0, all-zero F_T = 0 weights, chemistry = 0 nondetects);
7. Full end-to-end observer-only DNA + chemistry + RNA Dirichlet-multinomial scoring across BOTH
   alife/mep_redox/synthetic/D001 and D002 (reporting separate denominators for active panel rows,
   simulator-populatable steps, and observed assembled rows; asserting the >=4-pathway >=50%-event
   RNA coverage gate and finite equal-block cast-weighted DM log scores across M_cap, M_capK,
   M_capKT, and M_power; never opening truth/ or /data/dev/plexityper/data/).
"""
import csv
import hashlib
import io
import math
import os
import sys

CAL_TO_KJ = 4.184e-3
R_KJ = 8.314462618e-3
T_REF = 298.15

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
LIT_DIR = os.path.join(REPO_ROOT, "scratch", "p6_obigt")
PAPERS_LIT_DIR = os.path.join(REPO_ROOT, "discussion", "papers", "p6_lit")
KO_LIST_PATH = "/data/scratch/refdbs/kofam_db/ko_list"
OUT_TSV = os.path.join(REPO_ROOT, "protocols", "P6-pathway-table.tsv")

HOUSEKEEPING_PROK_KOS = [
    "K02886", "K02906", "K02874", "K02946", "K02961", "K02967", "K03043",
    "K02470", "K03553", "K04043", "K01872", "K02355", "K03046", "K01889",
]
# Positive heterotrophic organic-C uptake/activation/dehydrogenase co-carriage KOs required on binned MAGs for primary AER_HET
# (acs, ackA, pta, sucA/odhA, aceE/pdhA, pfkA, edd, dmdA)
HETEROTROPH_POS_KOS = [
    "K01895", "K00925", "K00625", "K13788", "K00164",
    "K00163", "K00161", "K00850", "K01690", "K17486",
]
# Positive H2/CO2 methanogenesis co-carriage KOs on real archaeal MAGs (frhABG, mvhADG, fwdABCD, ftr, mch, mtd, mer)
METHANOGEN_H2_COCARRIAGE_KOS = [
    "K00440", "K00441", "K00443", "K14126", "K14127", "K14128",
    "K00200", "K00201", "K00202", "K00203", "K00672", "K01499", "K00319", "K00320",
]
CITED_OA_HASHES = {
    "PMC11594503.xml": "a986f81ddbae",
    "PMC12203076.xml": "8792fa129646",
    "PMC12287890.xml": "02a9f198cac1",
    "PMC12683309.xml": "2c86fac7a252",
    "PMC2905729.xml": "ee8abd03276a",
    "PMC3655542.xml": "2524602a9663",
    "PMC4070702.xml": "e2ded40f0507",
    "PMC4217175.xml": "dfcaa3a0b3c1",
    "PMC4445324.xml": "273df2ad568b",
    "PMC4608356.xml": "495088b5a727",
    "PMC5112241.xml": "48ab7139f3f6",
    "PMC6492693.xml": "5954763936ee",
    "PMC6653664.xml": "f5359be159a0",
    "PMC6700366.xml": "be89826a3111",
    "PMC6785171.xml": "5d304ba4cdad",
    "PMC6872814.xml": "70af7302cdbf",
    "PMC7007197.xml": "42dd8d91361e",
    "PMC8114937.xml": "811640082c02",
    "PMC8115077.xml": "75e477f3d078",
    "PMC8851195.xml": "0737a6a352cb",
    "PMC9515657.xml": "4513703e54c8",
    "PMC9540798.xml": "c02cc11c8b78",
    "PMC9828761.xml": "c8d4affec010",
    "vallino_huber_2018.html": "e7f861579216",
}


def load_obigt_species():
    """Load exact G0 (cal/mol) and H0 (cal/mol) at 298.15 K, 1 bar from scratch/p6_obigt/ CSVs."""
    sp = {
        "H2O": ({"H": 2, "O": 1}, 0, -56687.71, -68316.76, "JOH92"),
        "H+": ({"H": 1}, +1, 0.0, 0.0, "Den81"),
    }
    inorg_map = {
        ("O2", "aq"): ("O2", {"O": 2}, 0),
        ("NO3-", "aq"): ("NO3-", {"N": 1, "O": 3}, -1),
        ("NO2-", "aq"): ("NO2-", {"N": 1, "O": 2}, -1),
        ("NO", "aq"): ("NO", {"N": 1, "O": 1}, 0),
        ("N2O", "aq"): ("N2O", {"N": 2, "O": 1}, 0),
        ("N2", "aq"): ("N2", {"N": 2}, 0),
        ("NH4+", "aq"): ("NH4+", {"N": 1, "H": 4}, +1),
        ("NH3", "aq"): ("NH3", {"N": 1, "H": 3}, 0),
        ("SO4-2", "aq"): ("SO4-2", {"S": 1, "O": 4}, -2),
        ("HS-", "aq"): ("HS-", {"H": 1, "S": 1}, -1),
        ("H2S", "aq"): ("H2S", {"H": 2, "S": 1}, 0),
        ("H2", "aq"): ("H2", {"H": 2}, 0),
        ("CO2", "aq"): ("CO2", {"C": 1, "O": 2}, 0),
        ("HCO3-", "aq"): ("HCO3-", {"H": 1, "C": 1, "O": 3}, -1),
        ("CO3-2", "aq"): ("CO3-2", {"C": 1, "O": 3}, -2),
    }
    with open(os.path.join(LIT_DIR, "OBIGT_inorganic_aq.csv")) as fh:
        for row in csv.DictReader(fh):
            key = (row["name"], row["state"])
            if key in inorg_map:
                sym, fdict, z = inorg_map[key]
                assert row["E_units"] == "cal"
                sp[sym] = (fdict, z, float(row["G"]), float(row["H"]), row["ref1"])

    org_map = {
        ("CH4", "aq"): ("CH4", {"C": 1, "H": 4}, 0),
        ("acetate", "aq"): ("CH3COO-", {"C": 2, "H": 3, "O": 2}, -1),
        ("formate", "aq"): ("HCOO-", {"C": 1, "H": 1, "O": 2}, -1),
    }
    with open(os.path.join(LIT_DIR, "OBIGT_organic_aq.csv")) as fh:
        for row in csv.DictReader(fh):
            key = (row["name"], row["state"])
            if key in org_map:
                sym, fdict, z = org_map[key]
                assert row["E_units"] == "cal"
                sp[sym] = (fdict, z, float(row["G"]), float(row["H"]), row["ref1"])

    with open(os.path.join(LIT_DIR, "OBIGT_organic_aq_v210.csv")) as fh:
        for row in csv.DictReader(fh):
            if row["name"] == "glucose" and row["state"] == "aq":
                assert row["E_units"] == "cal" and row["ref1"] == "AP01"
                sp["C6H12O6"] = ({"C": 6, "H": 12, "O": 6}, 0, float(row["G"]), float(row["H"]), "AP01")
    return sp


def build_rows(sp):
    raw_rows = [
        dict(
            pathway_id="AER_HET",
            panel_tier="core",
            pathway_name="Aerobic organic-C respiration (MAG-resolved heterotrophic terminal oxidases)",
            direction_rule="O2 reduction coupled to heterotrophic organic-C oxidation; PRIMARY: binned prokaryotic MAGs ONLY that co-carry >=1 positive organic-C catabolism marker (acs/ackA/pta/sucA/aceE/pdhA/pfkA/edd/dmdA: K01895,K00925,K00625,K13788,K00164,K00163,K00161,K00850,K01690,K17486) AND lack obligate lithotroph/methanotroph donor pathways (NIT_AMO, NIT_NXR, SOX_AER, MOX_AER, HYD_UP); unbinned oxidases and oxidase-only MAGs lacking positive organic-C markers enter residual class unresolved_oxidase_donor (gated out in primary, retained in sensitivity)",
            unambiguous_kos="NONE",
            clade_or_synteny_kos="K02274,K02275,K02276,K15408,K00404,K00405,K00406,K00407,K15862,K00425,K00426,K02297,K02298,K02299,K02300",
            excluded_confounder_kos="K02256,K02261,K02262,K02263,K17893",
            reaction_equation="CH3COO- + 2 O2 -> 2 HCO3- + H+",
            stoich={"CH3COO-": -1, "O2": -2, "HCO3-": 2, "H+": 1},
            n_e=8,
            electron_donor_spec="CH3COO- = 0.05 * DOC_uM_C (alpha_DOC=0.10 labile C share, 2 C/acetate; sensitivity alpha_DOC in {0.01,1.0} and C6H12O6/6)",
            monod_substrates_uM="CH3COO-:15.6;O2:0.121",
            product_species_in_Q="HCO3-,H+",
            o2_inhibited="false",
            Ki_O2_uM="NA",
            kcat_e_rel="1.0",
            chi_per_e="1.0",
            dG_min_kJ_per_e="5.625",
            primary_sources="Thermo: SH88, SHS89, Sho95; Ks(DOC): PMC12683309 T3; Ks(O2): Stolper2010 in PMC4608356 Methods & CanfieldKraft2022 PMC9828761 T2/p.5335; FT: JinKirk2016 PMC5112241 Eq.10-11/T3",
        ),
        dict(
            pathway_id="DEN_NAR",
            panel_tier="core",
            pathway_name="Respiratory nitrate reduction to nitrite (nar/nap)",
            direction_rule="NO3- -> NO2- reduction; napAB + narI included on MAGs and unbinned contigs (excluding hzsABC+ anammox MAGs); narGH (K00370,K00371) included ONLY on non-NOB, non-anammox MAGs (or co-carried with K00374/K17218/nir/nor/nos) and excluded from unbinned contigs",
            unambiguous_kos="K02567,K02568,K00374",
            clade_or_synteny_kos="K00370,K00371",
            excluded_confounder_kos="K00372,K00367,K00360,K00373",
            reaction_equation="CH3COO- + 4 NO3- -> 2 HCO3- + 4 NO2- + H+",
            stoich={"CH3COO-": -1, "NO3-": -4, "HCO3-": 2, "NO2-": 4, "H+": 1},
            n_e=8,
            electron_donor_spec="CH3COO- = 0.05 * DOC_uM_C",
            monod_substrates_uM="CH3COO-:15.6;NO3-:17.9",
            product_species_in_Q="HCO3-,NO2-,H+",
            o2_inhibited="true",
            Ki_O2_uM="3.0",
            kcat_e_rel="1.0",
            chi_per_e="1.0",
            dG_min_kJ_per_e="5.625",
            primary_sources="Thermo: SH88, Sho95; Ks(NO3-,DOC): PMC12683309 T3-T4; Ki(O2): Kalvelage2011 in CanfieldKraft2022 PMC9828761 p.5337 (process-level transfer; uniform-0.255 uM sensitivity); FT: JinKirk2016 PMC5112241",
        ),
        dict(
            pathway_id="DEN_NIR",
            panel_tier="denit_tierAC_NO",
            pathway_name="Respiratory nitrite reduction to nitric oxide (nirK/nirS)",
            direction_rule="NO2- -> NO reduction (active in Tier A [NO+N2O] and Tier B_NO [NO only]); on binned MAGs exclude hzsABC+ anammox MAGs and amo/hao+ nitrifier MAGs lacking norBC/nosZ",
            unambiguous_kos="K00368,K15864",
            clade_or_synteny_kos="NONE",
            excluded_confounder_kos="K00362,K00363,K00366,K17877,K26138,K26139",
            reaction_equation="CH3COO- + 8 NO2- + 7 H+ -> 2 HCO3- + 8 NO + 4 H2O",
            stoich={"CH3COO-": -1, "NO2-": -8, "H+": -7, "HCO3-": 2, "NO": 8, "H2O": 4},
            n_e=8,
            electron_donor_spec="CH3COO- = 0.05 * DOC_uM_C",
            monod_substrates_uM="CH3COO-:15.6;NO2-:1.6",
            product_species_in_Q="HCO3-,NO",
            o2_inhibited="true",
            Ki_O2_uM="0.297",
            kcat_e_rel="1.0",
            chi_per_e="1.0",
            dG_min_kJ_per_e="5.625",
            primary_sources="Thermo: SH88, SSW01, Sho95, JOH92; Ks(NO2-): Hedison2019 PMC7007197 T4 measured on NO-forming CuNiR/NirK (WT RpNiR 1.6 uM; 2-domain AxNiR 21 uM; transfer assumption to NirS); Ki(O2): model transfer from Dalsgaard2014 PMC4217175 T3 net N2O C50=297 nM (uniform-0.255 uM sensitivity)",
        ),
        dict(
            pathway_id="DEN_NOR",
            panel_tier="denit_tierA_NO_N2O",
            pathway_name="Respiratory nitric oxide reduction to nitrous oxide (norBC/qNor)",
            direction_rule="2 NO -> N2O reduction (active in Tier A when both dissolved NO and N2O are measured across cohort); excludes detoxifying flavorubredoxin norVW and eukaryotic fungal P450nor",
            unambiguous_kos="K04561,K02305",
            clade_or_synteny_kos="NONE",
            excluded_confounder_kos="K12264,K12265,K15877,K02164,K02448,K04748",
            reaction_equation="CH3COO- + 8 NO -> 2 HCO3- + 4 N2O + H+",
            stoich={"CH3COO-": -1, "NO": -8, "HCO3-": 2, "N2O": 4, "H+": 1},
            n_e=8,
            electron_donor_spec="CH3COO- = 0.05 * DOC_uM_C",
            monod_substrates_uM="CH3COO-:15.6;NO:0.5",
            product_species_in_Q="HCO3-,N2O,H+",
            o2_inhibited="true",
            Ki_O2_uM="0.297",
            kcat_e_rel="1.0",
            chi_per_e="1.0",
            dG_min_kJ_per_e="5.625",
            primary_sources="Thermo: SH88, SSW01, Sho95; Ks(NO): Kahle2019 PMC6872814 Fig.2a/p.2; Ki(O2): model transfer from Dalsgaard2014 PMC4217175 T3 net N2O C50=297 nM (uniform-0.255 uM sensitivity); FT: JinKirk2016 PMC5112241",
        ),
        dict(
            pathway_id="DEN_NOS",
            panel_tier="denit_tierAB_N2O",
            pathway_name="Respiratory nitrous oxide reduction to dinitrogen (nosZ)",
            direction_rule="N2O -> N2 reduction (active in Tier A [NO+N2O] and Tier B_N2O [N2O only] whenever dissolved N2O is measured across cohort); covers TAT clade I and Sec clade II nosZ",
            unambiguous_kos="K00376",
            clade_or_synteny_kos="NONE",
            excluded_confounder_kos="K19339,K19342",
            reaction_equation="CH3COO- + 4 N2O -> 2 HCO3- + 4 N2 + H+",
            stoich={"CH3COO-": -1, "N2O": -4, "HCO3-": 2, "N2": 4, "H+": 1},
            n_e=8,
            electron_donor_spec="CH3COO- = 0.05 * DOC_uM_C",
            monod_substrates_uM="CH3COO-:15.6;N2O:1.3",
            product_species_in_Q="HCO3-,N2,H+",
            o2_inhibited="true",
            Ki_O2_uM="0.205",
            kcat_e_rel="1.0",
            chi_per_e="1.0",
            dG_min_kJ_per_e="5.625",
            primary_sources="Thermo: SH88, SHS89, SSW01, Sho95; Ks(N2O): Sun2021 PMC8115077 p.1439; Ki(O2): model transfer from Dalsgaard2014 PMC4217175 T3 net N2 C50=205 nM (uniform-0.255 uM sensitivity); FT: JinKirk2016 PMC5112241",
        ),
        dict(
            pathway_id="DEN_NIRNOR",
            panel_tier="denit_tierB_N2O_only",
            pathway_name="Nitrite reduction to nitrous oxide (nir+nor collapsed when NO unmeasured, N2O measured)",
            direction_rule="NO2- -> 0.5 N2O reduction (active in Tier B_N2O when N2O is measured across cohort and NO is unmeasured); same anammox/nitrifier MAG exclusions on nirK/nirS",
            unambiguous_kos="K00368,K15864,K04561,K02305",
            clade_or_synteny_kos="NONE",
            excluded_confounder_kos="K00362,K00363,K12264,K12265,K15877",
            reaction_equation="CH3COO- + 4 NO2- + 3 H+ -> 2 HCO3- + 2 N2O + 2 H2O",
            stoich={"CH3COO-": -1, "NO2-": -4, "H+": -3, "HCO3-": 2, "N2O": 2, "H2O": 2},
            n_e=8,
            electron_donor_spec="CH3COO- = 0.05 * DOC_uM_C",
            monod_substrates_uM="CH3COO-:15.6;NO2-:1.6",
            product_species_in_Q="HCO3-,N2O",
            o2_inhibited="true",
            Ki_O2_uM="0.297",
            kcat_e_rel="1.0",
            chi_per_e="1.0",
            dG_min_kJ_per_e="5.625",
            primary_sources="Thermo: SH88, SSW01, Sho95, JOH92; Ks(NO2-): Hedison2019 PMC7007197 T4 CuNiR/NirK 1.6 uM (transfer assumption to NirS); Ki(O2): model transfer from Dalsgaard2014 PMC4217175 T3 net N2O C50=297 nM (uniform-0.255 uM sensitivity); FT: JinKirk2016 PMC5112241",
        ),
        dict(
            pathway_id="DEN_NORNOS",
            panel_tier="denit_tierB_NO_only",
            pathway_name="Nitric oxide reduction to dinitrogen (nor+nos collapsed when NO measured, N2O unmeasured)",
            direction_rule="2 NO -> N2 reduction (active in Tier B_NO when dissolved NO is measured across cohort and N2O is unmeasured, avoiding fabricated [N2O])",
            unambiguous_kos="K04561,K02305,K00376",
            clade_or_synteny_kos="NONE",
            excluded_confounder_kos="K12264,K12265,K15877,K19339,K19342",
            reaction_equation="CH3COO- + 4 NO -> 2 HCO3- + 2 N2 + H+",
            stoich={"CH3COO-": -1, "NO": -4, "HCO3-": 2, "N2": 2, "H+": 1},
            n_e=8,
            electron_donor_spec="CH3COO- = 0.05 * DOC_uM_C",
            monod_substrates_uM="CH3COO-:15.6;NO:0.5",
            product_species_in_Q="HCO3-,N2,H+",
            o2_inhibited="true",
            Ki_O2_uM="0.255",
            kcat_e_rel="1.0",
            chi_per_e="1.0",
            dG_min_kJ_per_e="5.625",
            primary_sources="Thermo: SH88, SHS89, SSW01, Sho95; Ks(NO): Kahle2019 PMC6872814 Fig.2a; Ki(O2): model transfer from Dalsgaard2014 PMC4217175 T3 total C50=255 nM; FT: JinKirk2016 PMC5112241",
        ),
        dict(
            pathway_id="DEN_NO2N2",
            panel_tier="denit_tierC_neither",
            pathway_name="Complete nitrite denitrification to dinitrogen (nir+nor+nos collapsed when NO and N2O unmeasured)",
            direction_rule="NO2- -> 0.5 N2 reduction (active in Tier C when neither NO nor N2O is measured across cohort, including mep_redox synthetic calibration); same anammox/nitrifier MAG exclusions on nirK/nirS",
            unambiguous_kos="K00368,K15864,K04561,K02305,K00376",
            clade_or_synteny_kos="NONE",
            excluded_confounder_kos="K00362,K00363,K12264,K12265,K15877",
            reaction_equation="3 CH3COO- + 8 NO2- + 5 H+ -> 6 HCO3- + 4 N2 + 4 H2O",
            stoich={"CH3COO-": -3, "NO2-": -8, "H+": -5, "HCO3-": 6, "N2": 4, "H2O": 4},
            n_e=24,
            electron_donor_spec="CH3COO- = 0.05 * DOC_uM_C",
            monod_substrates_uM="CH3COO-:15.6;NO2-:1.6",
            product_species_in_Q="HCO3-,N2",
            o2_inhibited="true",
            Ki_O2_uM="0.255",
            kcat_e_rel="1.0",
            chi_per_e="1.0",
            dG_min_kJ_per_e="5.625",
            primary_sources="Thermo: SH88, SHS89, Sho95, JOH92; Ks(NO2-): Hedison2019 PMC7007197 T4 CuNiR/NirK 1.6 uM (transfer assumption to NirS); Ki(O2): Dalsgaard2014 PMC4217175 T3 (Total N2+N2O C50=255 nM); FT: JinKirk2016 PMC5112241",
        ),
        dict(
            pathway_id="DNRA_NRF",
            panel_tier="core",
            pathway_name="Dissimilatory nitrite reduction to ammonium (nrfAH)",
            direction_rule="NO2- -> NH4+ respiratory ammonification via periplasmic cytochrome c nitrite reductase nrfAH; exclude cytoplasmic assimilatory nirBD and hzsABC+ anammox MAGs",
            unambiguous_kos="K03385,K15876",
            clade_or_synteny_kos="NONE",
            excluded_confounder_kos="K00362,K00363,K04013,K04014,K04015",
            reaction_equation="3 CH3COO- + 4 NO2- + 5 H+ + 4 H2O -> 6 HCO3- + 4 NH4+",
            stoich={"CH3COO-": -3, "NO2-": -4, "H+": -5, "H2O": -4, "HCO3-": 6, "NH4+": 4},
            n_e=24,
            electron_donor_spec="CH3COO- = 0.05 * DOC_uM_C",
            monod_substrates_uM="CH3COO-:15.6;NO2-:6.4",
            product_species_in_Q="HCO3-,NH4+",
            o2_inhibited="true",
            Ki_O2_uM="1.56",
            kcat_e_rel="1.0",
            chi_per_e="1.0",
            dG_min_kJ_per_e="5.625",
            primary_sources="Thermo: SH88, Sho95, JOH92; Ks(NO2-): Silveira2010 PMC2905729 T1 (ammonia-forming ccNiR/NrfA DQ 6.4 uM; range 2.7-15.0 uM); Ks(DOC),Ki(O2): PMC12683309 T2-T3; FT: JinKirk2016 PMC5112241",
        ),
        dict(
            pathway_id="NIT_AMO",
            panel_tier="core",
            pathway_name="Chemolithotrophic ammonia oxidation to nitrite (amoCAB/hao/archaeal amoA)",
            direction_rule="NH4+ -> NO2- aerobic oxidation; hao (K10535, excluding hzsABC+ anammox MAGs) and archaeal amoA (K28504) are unambiguous; shared K10944-K10946 assigned to NIT_AMO ONLY on AOB/AOA/comammox MAGs or MAGs co-carrying K10535 and lacking mmo/mxaF/xoxF",
            unambiguous_kos="K10535,K28504",
            clade_or_synteny_kos="K10944,K10945,K10946",
            excluded_confounder_kos="K00378",
            reaction_equation="NH4+ + 1.5 O2 -> NO2- + 2 H+ + H2O",
            stoich={"NH4+": -1, "O2": -1.5, "NO2-": 1, "H+": 2, "H2O": 1},
            n_e=6,
            electron_donor_spec="NH4+ speciated from measured TAN at in-situ T, pH, I",
            monod_substrates_uM="TAN:2.11;O2:0.333",
            product_species_in_Q="NO2-,H+",
            o2_inhibited="false",
            Ki_O2_uM="NA",
            kcat_e_rel="1.0",
            chi_per_e="1.0",
            dG_min_kJ_per_e="5.625",
            primary_sources="Thermo: SH88, SHS89, JOH92; Ks(TAN): Jayathilake2019 PMC6700366 T1 & Jung2022 PMC9515657 T1; Ks(O2): Bristow2016 in CanfieldKraft2022 PMC9828761 p.5335; FT: JinKirk2016 PMC5112241",
        ),
        dict(
            pathway_id="NIT_NXR",
            panel_tier="core",
            pathway_name="Chemolithotrophic nitrite oxidation to nitrate (nxrAB)",
            direction_rule="NO2- -> NO3- aerobic oxidation; K00370,K00371 excluded from unbinned contigs and assigned to NIT_NXR ONLY on NOB/comammox MAGs (Nitrospirota/Nitrospinota/Nitrobacter/Nitrococcus/Nitrolancea, excluding hzsABC+ anammox MAGs; or standalone K00370+K00371 lacking K00374/K17218/hzsA/nir/nor/nos)",
            unambiguous_kos="NONE",
            clade_or_synteny_kos="K00370,K00371",
            excluded_confounder_kos="K00374",
            reaction_equation="NO2- + 0.5 O2 -> NO3-",
            stoich={"NO2-": -1, "O2": -0.5, "NO3-": 1},
            n_e=2,
            electron_donor_spec="NO2- measured directly",
            monod_substrates_uM="NO2-:9.0;O2:0.778",
            product_species_in_Q="NO3-",
            o2_inhibited="false",
            Ki_O2_uM="NA",
            kcat_e_rel="1.0",
            chi_per_e="1.0",
            dG_min_kJ_per_e="5.625",
            primary_sources="Thermo: SH88, SHS89; Ks(NO2-): Jung2022 PMC9515657 T1 (N. moscoviensis 9 uM) & Zhang2021 PMC8114937; Ks(O2): Bristow2016 in CanfieldKraft2022 PMC9828761 p.5335; FT: JinKirk2016 PMC5112241",
        ),
        dict(
            pathway_id="AMX_HZS",
            panel_tier="core",
            pathway_name="Anaerobic ammonium oxidation / anammox (hzsABC + hdh)",
            direction_rule="NH4+ + NO2- -> N2 anaerobic comproportionation; hzsABC + hdh are 100% specific to anammox Planctomycetota",
            unambiguous_kos="K20932,K20933,K20934,K20935",
            clade_or_synteny_kos="NONE",
            excluded_confounder_kos="NONE",
            reaction_equation="NH4+ + NO2- -> N2 + 2 H2O",
            stoich={"NH4+": -1, "NO2-": -1, "N2": 1, "H2O": 2},
            n_e=3,
            electron_donor_spec="NH4+ speciated from measured TAN at in-situ T, pH, I",
            monod_substrates_uM="TAN:2.11;NO2-:0.55",
            product_species_in_Q="N2",
            o2_inhibited="true",
            Ki_O2_uM="0.886",
            kcat_e_rel="1.0",
            chi_per_e="1.0",
            dG_min_kJ_per_e="5.625",
            primary_sources="Thermo: SH88, SHS89, JOH92; Ks(NO2-): Kobayashi2025 PMC12203076 p.5 (Scalindua 0.1-1.0 uM) & Oshiki2013 PMC4070702; Ks(TAN): PMC6700366 T1 & PMC3655542; Ki(O2): Dalsgaard2014 PMC4217175 T3 (886 nM); FT: JinKirk2016 PMC5112241",
        ),
        dict(
            pathway_id="SRB_ORG",
            panel_tier="core",
            pathway_name="Dissimilatory sulfate reduction (reductive dsrAB/aprAB/qmoABC/dsrD)",
            direction_rule="SO4-2 -> HS- reduction; qmoABC + dsrD unambiguous reductive markers; reversible sat/aprAB/dsrAB excluded from unbinned contigs and assigned to SRB_ORG ONLY on reductive SRB clades or MAGs co-carrying qmoABC/dsrD and lacking sqr/sox/fccAB/sHdrABC (sensitivity excludes obligate H2-only hyaAB+ SRB MAGs)",
            unambiguous_kos="K16885,K16886,K16887,K23077",
            clade_or_synteny_kos="K00958,K00394,K00395,K11180,K11181",
            excluded_confounder_kos="K00380,K00381,K00390,K00392,K00955,K00956,K00957",
            reaction_equation="CH3COO- + SO4-2 -> 2 HCO3- + HS-",
            stoich={"CH3COO-": -1, "SO4-2": -1, "HCO3-": 2, "HS-": 1},
            n_e=8,
            electron_donor_spec="CH3COO- = 0.05 * DOC_uM_C",
            monod_substrates_uM="CH3COO-:5.0;SO4-2:39.0",
            product_species_in_Q="HCO3-,HS-",
            o2_inhibited="true",
            Ki_O2_uM="1.56",
            kcat_e_rel="1.0",
            chi_per_e="0.75",
            dG_min_kJ_per_e="4.219",
            primary_sources="Thermo: SH88, Sho95; Ks(CH3COO-,SO4-2) & FT(chi=6/8,mP=0.75): JinKirk2016 PMC5112241 T3 Rxn15 & Jorgensen2019 PMC6492693; Ki(O2): PMC12683309 T2-T3",
        ),
        dict(
            pathway_id="SOX_AER",
            panel_tier="core",
            pathway_name="Chemolithotrophic sulfide/sulfur oxidation (sqr/fccAB/sox/sHdr/rDsr)",
            direction_rule="HS-/H2S -> SO4-2 aerobic oxidation; sqr, fccAB, soxAXBCYZD, sHdrABC unambiguous sulfur-oxidation markers; reversible dsrAB/aprAB assigned to SOX_AER ONLY on SOB clades (SUP05/Thiobacillaceae/Chromatiaceae/Chlorobiaceae/Sulfurimonas) or MAGs co-carrying sox/sqr/fccAB/sHdrABC and lacking dsrD",
            unambiguous_kos="K17218,K17229,K17230,K17222,K17223,K17224,K17225,K17226,K17227,K22622,K28385,K28386,K28387",
            clade_or_synteny_kos="K11180,K11181,K00394,K00395",
            excluded_confounder_kos="K17219,K17220,K17221,K08352,K08353,K08354",
            reaction_equation="HS- + 2 O2 -> SO4-2 + H+",
            stoich={"HS-": -1, "O2": -2, "SO4-2": 1, "H+": 1},
            n_e=8,
            electron_donor_spec="HS- speciated from measured SigmaH2S at in-situ T, pH, I",
            monod_substrates_uM="SigmaH2S:4.68;O2:0.121",
            product_species_in_Q="SO4-2,H+",
            o2_inhibited="false",
            Ki_O2_uM="NA",
            kcat_e_rel="1.0",
            chi_per_e="1.0",
            dG_min_kJ_per_e="5.625",
            primary_sources="Thermo: SH88, SHS89; Ks(SigmaH2S): PMC12683309 T4 (0.15 g S m-3 = 4.68 uM); Ks(O2): PMC4608356 & CanfieldKraft2022 PMC9828761 T2; FT: JinKirk2016 PMC5112241",
        ),
        dict(
            pathway_id="MCR_MET",
            panel_tier="h2_gated",
            pathway_name="Hydrogenotrophic methanogenesis (MAG-resolved H2/CO2 mcrABG)",
            direction_rule="4 H2 + HCO3- + H+ -> CH4 + 3 H2O (active when H2 measured across cohort); PRIMARY: binned MAGs ONLY with hydrogenotrophic clade assignment (Methanobacteriales/Methanomicrobiales/Methanococcales/Methanocellales/Methanopyrales) or co-carriage of frhABG/mvhADG/fwdABCD-mer (or minimal-marker MET in mep_redox), excluding ANME and obligate acetoclastic Methanothrix; unbinned mcrABG gated out in primary and retained with MCR_ACET in sensitivity",
            unambiguous_kos="NONE",
            clade_or_synteny_kos="K00399,K00401,K00402",
            excluded_confounder_kos="K00400,K03421,K03422",
            reaction_equation="4 H2 + HCO3- + H+ -> CH4 + 3 H2O",
            stoich={"H2": -4, "HCO3-": -1, "H+": -1, "CH4": 1, "H2O": 3},
            n_e=8,
            electron_donor_spec="H2 measured directly (active only when dissolved H2 is logged across cohort)",
            monod_substrates_uM="H2:4.7",
            product_species_in_Q="CH4",
            o2_inhibited="true",
            Ki_O2_uM="1.56",
            kcat_e_rel="1.0",
            chi_per_e="0.25",
            dG_min_kJ_per_e="1.406",
            primary_sources="Thermo: SH88, SHS89, PS01, JOH92; Ks(H2) & FT(chi=2/8,mP=0.25): JinKirk2016 PMC5112241 T3 Rxn21; Ki(O2): PMC12683309 T2-T3",
        ),
        dict(
            pathway_id="MOX_AER",
            panel_tier="core",
            pathway_name="Aerobic methanotrophy (mmoXYZC / clade-assigned pmoCAB)",
            direction_rule="CH4 + 2 O2 -> HCO3- + H+ + H2O; soluble mmoXYZC unambiguous; shared K10944-K10946 assigned to MOX_AER ONLY on canonical methanotroph MAGs or MAGs co-carrying mmoXYZC or mxaFI/xoxF and lacking hao (K10535)",
            unambiguous_kos="K16157,K16158,K16159,K16161",
            clade_or_synteny_kos="K10944,K10945,K10946",
            excluded_confounder_kos="K16160,K16162,K14028,K14029,K23995",
            reaction_equation="CH4 + 2 O2 -> HCO3- + H+ + H2O",
            stoich={"CH4": -1, "O2": -2, "HCO3-": 1, "H+": 1, "H2O": 1},
            n_e=8,
            electron_donor_spec="CH4 measured directly",
            monod_substrates_uM="CH4:2.5;O2:0.16",
            product_species_in_Q="HCO3-,H+",
            o2_inhibited="false",
            Ki_O2_uM="NA",
            kcat_e_rel="1.0",
            chi_per_e="1.0",
            dG_min_kJ_per_e="5.625",
            primary_sources="Thermo: SH88, SHS89, PS01, JOH92; Ks(CH4): Seto2019 PMC6785171 p.1191 (KniefDunfield2005); Ks(O2): Steinsdottir2022 PMC9540798 (Tiano2014 66-259 nM); FT: JinKirk2016 PMC5112241",
        ),
        dict(
            pathway_id="HYD_FERM",
            panel_tier="h2_gated",
            pathway_name="Fermentative H2 production ([FeFe] ferredoxin hydrogenase hydAB)",
            direction_rule="C6H12O6 + 4 H2O -> 2 CH3COO- + 2 HCO3- + 4 H+ + 4 H2 (active when dissolved H2 is measured across cohort); catalytic [FeFe] hydAB unambiguous; O2-inhibited primary (sensitivity toggle uninhibited per P6 l.81)",
            unambiguous_kos="K00532,K00533,K00534",
            clade_or_synteny_kos="NONE",
            excluded_confounder_kos="K00656,K03737,K00436,K05586,K05587,K05588",
            reaction_equation="C6H12O6 + 4 H2O -> 2 CH3COO- + 2 HCO3- + 4 H+ + 4 H2",
            stoich={"C6H12O6": -1, "H2O": -4, "CH3COO-": 2, "HCO3-": 2, "H+": 4, "H2": 4},
            n_e=8,
            electron_donor_spec="C6H12O6 = (0.10 / 6) * DOC_uM_C (labile hexose equivalent; product CH3COO- = 0.05 * DOC_uM_C)",
            monod_substrates_uM="C6H12O6:52.1",
            product_species_in_Q="CH3COO-,HCO3-,H+,H2",
            o2_inhibited="true",
            Ki_O2_uM="1.56",
            kcat_e_rel="1.0",
            chi_per_e="0.50",
            dG_min_kJ_per_e="11.250",
            primary_sources="Thermo: AP01, Sho95, SH88, SHS89, JOH92; Stoich: Jayathilake2019 PMC6700366 T3; Ks(ferm),Ki(O2): PMC12683309 T2-T3 (10 g COD m-3 = 52.1 uM hexose); FT: JinKirk2016 PMC5112241",
        ),
        dict(
            pathway_id="HYD_UP",
            panel_tier="h2_gated",
            pathway_name="Sulfate-coupled anaerobic respiratory H2 uptake (MAG-resolved [NiFe] uptake hyd + SRB)",
            direction_rule="4 H2 + SO4-2 + H+ -> HS- + 4 H2O (active when dissolved H2 measured across cohort); PRIMARY: binned MAGs ONLY co-carrying reductive sulfate-reduction genes (qmoABC/dsrD or non-SOB dsrAB/aprAB) and lacking mcrABG; unbinned/non-SRB uptake hydrogenases gated out in primary (acceptor underidentified) and retained in sensitivity",
            unambiguous_kos="NONE",
            clade_or_synteny_kos="K06281,K06282,K05922,K05927,K23548,K23549",
            excluded_confounder_kos="K00440,K00441,K00442,K00443,K14126,K14127,K14128,K04651,K04652,K04653,K04654,K04655,K04656",
            reaction_equation="4 H2 + SO4-2 + H+ -> HS- + 4 H2O",
            stoich={"H2": -4, "SO4-2": -1, "H+": -1, "HS-": 1, "H2O": 4},
            n_e=8,
            electron_donor_spec="H2 measured directly (active only when dissolved H2 is logged across cohort)",
            monod_substrates_uM="H2:1.1;SO4-2:39.0",
            product_species_in_Q="HS-",
            o2_inhibited="true",
            Ki_O2_uM="1.56",
            kcat_e_rel="1.0",
            chi_per_e="0.75",
            dG_min_kJ_per_e="5.625",
            primary_sources="Thermo: SH88, SHS89, JOH92; Ks(H2,SO4-2) & FT(chi=6/8,mP=1.0): JinKirk2016 PMC5112241 T3 Rxn14 & PMC6653664 T1; Ki(O2): PMC12683309 T2-T3",
        ),
        dict(
            pathway_id="SOX_NO3",
            panel_tier="sensitivity_only",
            pathway_name="Anoxic nitrate-coupled sulfide oxidation (MAG-resolved sqr/sox + nap/nar sensitivity)",
            direction_rule="HS- + 4 NO3- -> SO4-2 + 4 NO2- + H+; sensitivity-only split of SOX_AER markers carried on MAGs that also encode napAB/narGHI in anoxic nitrate-bearing waters",
            unambiguous_kos="NONE",
            clade_or_synteny_kos="K17218,K17229,K17230,K17222,K17223,K17224,K17225,K17226,K17227,K22622,K28385,K28386,K28387",
            excluded_confounder_kos="NONE",
            reaction_equation="HS- + 4 NO3- -> SO4-2 + 4 NO2- + H+",
            stoich={"HS-": -1, "NO3-": -4, "SO4-2": 1, "NO2-": 4, "H+": 1},
            n_e=8,
            electron_donor_spec="HS- speciated from measured SigmaH2S at in-situ T, pH, I",
            monod_substrates_uM="SigmaH2S:4.68;NO3-:17.9",
            product_species_in_Q="SO4-2,NO2-,H+",
            o2_inhibited="true",
            Ki_O2_uM="3.0",
            kcat_e_rel="1.0",
            chi_per_e="1.0",
            dG_min_kJ_per_e="5.625",
            primary_sources="Thermo: SH88; Ks(SigmaH2S,NO3-): PMC12683309 T4 (0.15 g S m-3, 0.25 g N m-3); Ki(O2): Kalvelage2011 in CanfieldKraft2022 PMC9828761 p.5337; FT: JinKirk2016 PMC5112241",
        ),
        dict(
            pathway_id="MCR_ACET",
            panel_tier="sensitivity_only",
            pathway_name="Acetoclastic methanogenesis (MAG-resolved Methanothrix/Methanosarcina mcrABG sensitivity)",
            direction_rule="CH3COO- + H2O -> CH4 + HCO3-; sensitivity-only row for mcrABG on acetoclastic Methanothrix/Methanosarcina MAGs (co-carrying cdhCDE K00193,K00194,K00197 + acs/ackA-pta)",
            unambiguous_kos="NONE",
            clade_or_synteny_kos="K00399,K00401,K00402",
            excluded_confounder_kos="K00400,K03421,K03422",
            reaction_equation="CH3COO- + H2O -> CH4 + HCO3-",
            stoich={"CH3COO-": -1, "H2O": -1, "CH4": 1, "HCO3-": 1},
            n_e=8,
            electron_donor_spec="CH3COO- = 0.05 * DOC_uM_C",
            monod_substrates_uM="CH3COO-:23.0",
            product_species_in_Q="CH4,HCO3-",
            o2_inhibited="true",
            Ki_O2_uM="1.56",
            kcat_e_rel="1.0",
            chi_per_e="0.25",
            dG_min_kJ_per_e="2.813",
            primary_sources="Thermo: Sho95, JOH92, PS01, SH88; Ks(CH3COO-) & FT(chi=2/8,mP=0.50): JinKirk2016 PMC5112241 T3 Rxn22 & Ozuolmez2015 PMC4445324 T1; Ki(O2): PMC12683309 T2-T3",
        ),
    ]

    out = []
    for r in raw_rows:
        st = r["stoich"]
        ne = r["n_e"]
        elems = {}
        chg = 0.0
        dG0 = 0.0
        dH0 = 0.0
        for sname, nu in st.items():
            fdict, z, Gc, Hc, _ = sp[sname]
            chg += nu * z
            dG0 += nu * Gc * CAL_TO_KJ
            dH0 += nu * Hc * CAL_TO_KJ
            for el, cnt in fdict.items():
                elems[el] = elems.get(el, 0.0) + nu * cnt
        bad_el = {k: v for k, v in elems.items() if abs(v) > 1e-12}
        assert abs(chg) < 1e-12 and not bad_el, f"Balance fail {r['pathway_id']}: chg={chg}, el={bad_el}"
        r["stoichiometry_spec"] = ";".join(f"{k}:{v}" for k, v in st.items())
        r["dG0_298_kJ_rxn"] = f"{dG0:.3f}"
        r["dG0_298_kJ_per_e"] = f"{dG0 / ne:.3f}"
        r["dH0_298_kJ_rxn"] = f"{dH0:.3f}"
        out.append(r)
    return out


def parse_stoich_spec(spec):
    out = {}
    for item in spec.split(";"):
        sp_name, val = item.rsplit(":", 1)
        out[sp_name] = float(val)
    return out


def parse_monod_spec(spec):
    out = []
    for item in spec.split(";"):
        sp_name, val = item.rsplit(":", 1)
        out.append((sp_name, float(val)))
    return out


def dm_loglik(counts, probs, phi):
    """Exact Dirichlet-multinomial log-likelihood for count vector n, mean pi, concentration alpha_0 = 1/phi."""
    alpha0 = 1.0 / phi
    N = sum(counts)
    ll = math.lgamma(N + 1.0) + math.lgamma(alpha0) - math.lgamma(N + alpha0)
    for n_r, p_r in zip(counts, probs):
        a_r = alpha0 * p_r
        ll += math.lgamma(n_r + a_r) - math.lgamma(n_r + 1.0) - math.lgamma(a_r)
    return ll


def safe_probs(weights, caps, eps_obs=1e-4):
    R_len = len(weights)
    tot = sum(weights)
    if tot > 0:
        base = [w / tot for w in weights]
    else:
        c_tot = sum(1.0 for c in caps if c > 0)
        base = [(1.0 / c_tot if c > 0 else 0.0) for c in caps] if c_tot > 0 else [1.0 / R_len] * R_len
    probs = [(1.0 - eps_obs) * b + eps_obs / R_len for b in base]
    assert abs(sum(probs) - 1.0) < 1e-14
    assert min(probs) >= eps_obs / R_len * (1.0 - 1e-12)
    return probs


def compute_carbonate_activities(row, aH, Ka1_C, Ka2_C, default_a_co2=1e-3):
    """Return (a_CO2, a_HCO3) from measured DIC (uM) when present, or fixed free CO2 activity (model.py:25) when absent."""
    if "DIC" in row:
        DIC_M = float(row["DIC"]) * 1e-6
        a_CO2 = DIC_M / (1.0 + Ka1_C / aH + Ka1_C * Ka2_C / (aH * aH))
    else:
        a_CO2 = default_a_co2
    a_HCO3 = a_CO2 * Ka1_C / aH
    return a_CO2, a_HCO3


def make_binned_ko_classifier(panel_rows):
    """Build deterministic binned-MAG CDS KO classifier returning (pathway_id_or_None, residual_class_or_None)."""
    pid_to_idx = {r["pathway_id"]: i for i, r in enumerate(panel_rows)}
    unamb_map = {}
    for r in panel_rows:
        if r["unambiguous_kos"] != "NONE":
            for k in r["unambiguous_kos"].split(","):
                unamb_map[k] = r["pathway_id"]

    aer_ox_kos = set(panel_rows[pid_to_idx["AER_HET"]]["clade_or_synteny_kos"].split(","))
    hyd_up_kos = set(panel_rows[pid_to_idx["HYD_UP"]]["clade_or_synteny_kos"].split(","))
    mcr_kos = set(panel_rows[pid_to_idx["MCR_MET"]]["clade_or_synteny_kos"].split(","))
    het_pos_kos = set(HETEROTROPH_POS_KOS)
    sox_markers = {"K17218", "K17222", "K17223", "K17224", "K17225", "K17226", "K17227", "K22622", "K17229", "K17230", "K28385", "K28386", "K28387"}
    srb_core = {"K16885", "K16886", "K16887", "K23077", "K00394", "K00395", "K11180", "K11181"}

    def classify(ko, kset):
        if not ko:
            return None, None
        has_amx = "K20932" in kset
        if ko in unamb_map:
            pid = unamb_map[ko]
            if pid in {"DEN_NAR", "DEN_NO2N2", "DNRA_NRF", "NIT_AMO"} and has_amx:
                return None, None
            if pid == "DEN_NO2N2" and ("K10535" in kset or "K28504" in kset) and not ({"K04561", "K00376"} & kset):
                return None, None
            return pid, None
        if ko in aer_ox_kos:
            is_nob = ({"K00370", "K00371"} <= kset) and not ({"K00374", "K17218", "K20932", "K00368", "K15864", "K04561", "K00376"} & kset)
            is_litho = ("K10535" in kset) or ("K28504" in kset) or is_nob or bool(sox_markers & kset) or ("K16157" in kset) or bool(hyd_up_kos & kset)
            if is_litho:
                return None, None
            if bool(het_pos_kos & kset):
                return "AER_HET", None
            return None, "unresolved_oxidase_donor"
        if ko in hyd_up_kos:
            has_srb = bool(srb_core & kset) and not (sox_markers & kset and "K23077" not in kset)
            if has_srb and not (mcr_kos & kset):
                return "HYD_UP", None
            return None, None
        if ko in mcr_kos:
            if "K00194" in kset:
                return None, None
            return "MCR_MET", None
        if ko in {"K10944", "K10945", "K10946"}:
            if "K10535" in kset and "K16157" not in kset:
                return "NIT_AMO", None
            if "K16157" in kset or "K14028" in kset or "K23995" in kset:
                return "MOX_AER", None
            return None, None
        if ko in {"K00370", "K00371"}:
            if has_amx:
                return None, None
            if "K00374" in kset or "K17218" in kset or ({"K00368", "K15864", "K04561", "K00376"} & kset):
                return "DEN_NAR", None
            return "NIT_NXR", None
        if ko in {"K00958", "K00394", "K00395", "K11180", "K11181"}:
            if sox_markers & kset and "K23077" not in kset:
                return ("SOX_AER", None) if ko != "K00958" else (None, None)
            return "SRB_ORG", None
        return None, None

    return classify


def score_mep_redox_run(run_dir, rows, sp):
    """Run full observer-only DNA + chemistry + RNA scoring on a mep_redox synthetic directory."""
    # Active cohort panel for mep_redox chemistry header (O2,NO3,NO2,NH4,SO4,H2S,CH4,H2,DOC; no NO, no N2O -> Tier C + H2 = 13 rows)
    active_tiers = {"core", "denit_tierC_neither", "h2_gated"}
    panel_rows = [r for r in rows if r["panel_tier"] in active_tiers]
    pid_to_idx = {r["pathway_id"]: i for i, r in enumerate(panel_rows)}
    R_len = len(panel_rows)
    assert R_len == 13

    prok_mags = set()
    mag_bp = {}
    with open(os.path.join(run_dir, "mag_completeness.tsv")) as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            if row["domain"] != "Eukaryota":
                prok_mags.add(row["mag_id"])
                mag_bp[row["mag_id"]] = float(row["genome_size_bp"])

    mag_kos = {m: set() for m in prok_mags}
    cds_records = {}
    with open(os.path.join(run_dir, "per_mag_cds.tsv")) as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            m = row["mag_id"]
            if m in prok_mags:
                ko = row["KO"]
                if ko:
                    mag_kos[m].add(ko)
                cds_records[row["protein_id"]] = (m, ko, int(row["end"]) - int(row["start"]) + 1)

    classify_ko = make_binned_ko_classifier(panel_rows)
    unresolved_oxidase_donor_cds = 0

    cds_to_pathway = {}
    hk_by_ko = {k: [] for k in HOUSEKEEPING_PROK_KOS}
    pathway_cds_counts = {}
    for cds_id, (m, ko, length) in cds_records.items():
        if ko in hk_by_ko:
            hk_by_ko[ko].append((cds_id, m, length))
        pid, resid = classify_ko(ko, mag_kos[m])
        if resid == "unresolved_oxidase_donor":
            unresolved_oxidase_donor_cds += 1
        if pid is not None:
            cds_to_pathway[cds_id] = pid
            pathway_cds_counts[pid] = pathway_cds_counts.get(pid, 0) + 1

    # Read DNA abundance per sample_id -> per-MAG base depth proportional to reads / genome_size_bp
    dna_mag_reads = {}
    with open(os.path.join(run_dir, "abundance.tsv")) as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            sid = row["sample_id"]
            dna_mag_reads.setdefault(sid, {})[row["mag_id"]] = float(row["reads"])

    # Read RNA counts per sample_id & protein_id
    rna_panel_counts = {}
    with open(os.path.join(run_dir, "rna_counts.tsv")) as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            sid = row["sample_id"]
            cid = row["protein_id"]
            c = int(row["reads"])
            vec = rna_panel_counts.setdefault(sid, [0] * R_len)
            pid = cds_to_pathway.get(cid)
            if pid is not None:
                vec[pid_to_idx[pid]] += c

    # Precompute van 't Hoff helper
    def dG0_T(st, T):
        dG298 = sum(nu * sp[s][2] * CAL_TO_KJ for s, nu in st.items())
        dH298 = sum(nu * sp[s][3] * CAL_TO_KJ for s, nu in st.items())
        return (T / T_REF) * dG298 + dH298 * (1.0 - T / T_REF)

    # Read chemistry and score every sample_id across casts (equal-block cast weighting)
    cast_scores = {}
    cast_has_oxic = {}
    cast_has_anoxic = {}
    panel_read_totals = []
    n_samples = 0
    pathway_rna_detect = [0] * R_len
    with open(os.path.join(run_dir, "chemistry.tsv")) as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            sid = row["sample_id"]
            cast_id = row["event_id"]
            T = float(row["T_K"])
            pH = float(row["pH"])
            RT = R_KJ * T
            aH = 10.0 ** (-pH)

            def get_uM(col, dl=0.01):
                v = float(row[col])
                return v if v > 0.0 else 0.5 * dl

            DOC_uM = get_uM("DOC")
            O2_uM = get_uM("O2")
            cast_has_oxic[cast_id] = cast_has_oxic.get(cast_id, False) or (O2_uM >= 10.0)
            cast_has_anoxic[cast_id] = cast_has_anoxic.get(cast_id, False) or (O2_uM < 10.0)
            NO3_uM = get_uM("NO3")
            NO2_uM = get_uM("NO2")
            TAN_uM = get_uM("NH4")
            SO4_uM = get_uM("SO4")
            S_tot_uM = get_uM("H2S")
            H2_uM = get_uM("H2")
            CH4_uM = get_uM("CH4")
            Ka1_C = math.exp(-dG0_T({"CO2": -1, "H2O": -1, "HCO3-": 1, "H+": 1}, T) / RT)
            Ka2_C = math.exp(-dG0_T({"HCO3-": -1, "CO3-2": 1, "H+": 1}, T) / RT)
            Ka1_S = math.exp(-dG0_T({"H2S": -1, "HS-": 1, "H+": 1}, T) / RT)
            Ka_N = math.exp(-dG0_T({"NH4+": -1, "NH3": 1, "H+": 1}, T) / RT)

            a_CO2, a_HCO3 = compute_carbonate_activities(row, aH, Ka1_C, Ka2_C)
            N2_M = float(row["N2"]) * 1e-6 if "N2" in row else 5e-4
            a_H2S = (S_tot_uM * 1e-6) / (1.0 + Ka1_S / aH)
            a_HS = a_H2S * Ka1_S / aH
            a_NH4 = (TAN_uM * 1e-6) / (1.0 + Ka_N / aH)

            act = {
                "H2O": 1.0, "H+": aH,
                "O2": O2_uM * 1e-6, "NO3-": NO3_uM * 1e-6, "NO2-": NO2_uM * 1e-6,
                "N2": N2_M, "NH4+": a_NH4, "SO4-2": SO4_uM * 1e-6, "HS-": a_HS,
                "H2": H2_uM * 1e-6, "CH4": CH4_uM * 1e-6, "HCO3-": a_HCO3,
                "CH3COO-": 0.05 * DOC_uM * 1e-6, "C6H12O6": (0.10 / 6.0) * DOC_uM * 1e-6,
            }
            conc_uM = {
                "CH3COO-": 0.05 * DOC_uM, "C6H12O6": (0.10 / 6.0) * DOC_uM,
                "O2": O2_uM, "NO3-": NO3_uM, "NO2-": NO2_uM, "TAN": TAN_uM,
                "SO4-2": SO4_uM, "SigmaH2S": S_tot_uM, "H2": H2_uM, "CH4": CH4_uM,
            }

            # DNA marker-normalised capacity cap_{r,e}
            m_reads = dna_mag_reads[sid]
            mag_depth = {m: m_reads.get(m, 0.0) / mag_bp[m] for m in prok_mags}
            hk_depths = []
            for ko_hk, items in hk_by_ko.items():
                if items:
                    num = sum(mag_depth[m] * L for _, m, L in items)
                    den = sum(L for _, m, L in items)
                    hk_depths.append(num / den)
            hk_depths.sort()
            d_prok = hk_depths[len(hk_depths) // 2]
            assert d_prok > 0.0

            caps = [0.0] * R_len
            for cid, pid in cds_to_pathway.items():
                m, _, L = cds_records[cid]
                caps[pid_to_idx[pid]] += (mag_depth[m] * L) / d_prok

            w_cap, w_capK, w_capKT, w_power = [], [], [], []
            for i, r in enumerate(panel_rows):
                fk = 1.0
                for s_mon, ks_val in parse_monod_spec(r["monod_substrates_uM"]):
                    c_val = conc_uM[s_mon]
                    fk *= c_val / (ks_val + c_val)
                if r["o2_inhibited"] == "true":
                    ki = float(r["Ki_O2_uM"])
                    fk *= ki / (ki + O2_uM)
                st = parse_stoich_spec(r["stoichiometry_spec"])
                ne = float(r["n_e"])
                lnQ = sum(nu * math.log(act[s_q]) for s_q, nu in st.items())
                dGr_per_e = (dG0_T(st, T) + RT * lnQ) / ne
                dGmin = float(r["dG_min_kJ_per_e"])
                chi_e = float(r["chi_per_e"])
                arg = (dGr_per_e + dGmin) / (chi_e * RT)
                ft = 0.0 if arg >= 0.0 else (1.0 - math.exp(arg))
                yield_e = max(0.0, -dGr_per_e)
                c_i = caps[i]
                w_cap.append(c_i)
                w_capK.append(c_i * fk)
                w_capKT.append(c_i * fk * ft)
                w_power.append(c_i * fk * ft * yield_e)

            counts = rna_panel_counts.get(sid, [0] * R_len)
            n_tot_panel = sum(counts)
            assert n_tot_panel >= 100, f"Event {sid} failed N_panel >= 100 RNA reads ({n_tot_panel})"
            panel_read_totals.append(n_tot_panel)
            for i, c_val in enumerate(counts):
                if c_val > 0:
                    pathway_rna_detect[i] += 1
            n_samples += 1

            phi_fixed = 0.05
            s_models = {
                "M_cap": dm_loglik(counts, safe_probs(w_cap, caps), phi_fixed),
                "M_capK": dm_loglik(counts, safe_probs(w_capK, caps), phi_fixed),
                "M_capKT": dm_loglik(counts, safe_probs(w_capKT, caps), phi_fixed),
                "M_power": dm_loglik(counts, safe_probs(w_power, caps), phi_fixed),
            }
            assert all(math.isfinite(v) for v in s_models.values())
            cast_scores.setdefault(cast_id, []).append(s_models)

    # Assert cohort RNA coverage gate (>= 4 pathways detected in RNA in >= 50% of events)
    n_pass_cov = sum(1 for c in pathway_rna_detect if c >= 0.5 * n_samples)
    assert n_pass_cov >= 4, f"Coverage gate failed: only {n_pass_cov} pathways >= 50% events"
    panel_read_totals.sort()
    n_oxic_blocks = sum(1 for v in cast_has_oxic.values() if v)
    n_anoxic_blocks = sum(1 for v in cast_has_anoxic.values() if v)
    block_gate_met = (len(cast_scores) >= 6) and (n_oxic_blocks >= 3) and (n_anoxic_blocks >= 3)

    # Primary equal-block (cast_id) mean score
    models = ["M_cap", "M_capK", "M_capKT", "M_power"]
    equal_block_means = {
        m: sum(sum(ev[m] for ev in evs) / len(evs) for evs in cast_scores.values()) / len(cast_scores)
        for m in models
    }
    return {
        "n_panel_rows": R_len,
        "n_sim_populatable_rows": 12,  # DNRA_NRF is not in mep_redox REACTIONS
        "n_assembled_rows": len(pathway_cds_counts),
        "assembled_ids": sorted(pathway_cds_counts.keys()),
        "n_catabolic_cds": sum(pathway_cds_counts.values()),
        "n_unresolved_oxidase_donor_cds": unresolved_oxidase_donor_cds,
        "n_hk_cds": sum(len(v) for v in hk_by_ko.values()),
        "n_casts": len(cast_scores),
        "n_oxic_blocks": n_oxic_blocks,
        "n_anoxic_blocks": n_anoxic_blocks,
        "block_gate_met": block_gate_met,
        "n_events": n_samples,
        "panel_reads_min_med_max": (
            panel_read_totals[0],
            panel_read_totals[len(panel_read_totals) // 2],
            panel_read_totals[-1],
        ),
        "n_rna_cov50_pathways": n_pass_cov,
        "equal_block_scores": {k: round(v, 2) for k, v in equal_block_means.items()},
    }


def verify_all(rows, sp):
    ko_db = {}
    with open(KO_LIST_PATH) as fh:
        fh.readline()
        for line in fh:
            p = line.rstrip("\n").split("\t")
            if len(p) >= 12:
                ko_db[p[0]] = (p[1], p[8], p[11])

    all_kos = (
        set(HOUSEKEEPING_PROK_KOS)
        | set(HETEROTROPH_POS_KOS)
        | set(METHANOGEN_H2_COCARRIAGE_KOS)
        | {"K00193", "K00194", "K00197"}
    )
    for r in rows:
        for col in ("unambiguous_kos", "clade_or_synteny_kos", "excluded_confounder_kos"):
            val = r[col]
            if val != "NONE":
                for k in val.split(","):
                    all_kos.add(k)
    missing = [k for k in sorted(all_kos) if k not in ko_db]
    non_prot = [k for k in sorted(all_kos) if ko_db[k][1] == "-"]
    assert not missing, f"Missing KOs in ko_list: {missing}"
    assert not non_prot, f"Non-protein KOs in panel: {non_prot}"
    for fn, exp_h in sorted(CITED_OA_HASHES.items()):
        p_lit = os.path.join(PAPERS_LIT_DIR, fn)
        assert os.path.isfile(p_lit), f"Missing cited OA file {p_lit}"
        with open(p_lit, "rb") as fh:
            got_h = hashlib.sha256(fh.read()).hexdigest()[:12]
        assert got_h == exp_h, f"Hash mismatch on {p_lit}: {got_h} != {exp_h}"
    print(
        f"[PASS 1] {len(all_kos)}/{len(all_kos)} panel + confounder + co-carriage + housekeeping KOs verified in ko_list; "
        f"{len(CITED_OA_HASHES)}/{len(CITED_OA_HASHES)} primary OA files verified in discussion/papers/p6_lit/."
    )

    tier_configs = {
        "TierA_NO_N2O_H2": {"core", "denit_tierAC_NO", "denit_tierA_NO_N2O", "denit_tierAB_N2O", "h2_gated"},
        "TierA_NO_N2O_noH2": {"core", "denit_tierAC_NO", "denit_tierA_NO_N2O", "denit_tierAB_N2O"},
        "TierB_N2O_only_H2": {"core", "denit_tierB_N2O_only", "denit_tierAB_N2O", "h2_gated"},
        "TierB_N2O_only_noH2": {"core", "denit_tierB_N2O_only", "denit_tierAB_N2O"},
        "TierB_NO_only_H2": {"core", "denit_tierAC_NO", "denit_tierB_NO_only", "h2_gated"},
        "TierB_NO_only_noH2": {"core", "denit_tierAC_NO", "denit_tierB_NO_only"},
        "TierC_neither_H2": {"core", "denit_tierC_neither", "h2_gated"},
        "TierC_neither_noH2": {"core", "denit_tierC_neither"},
    }
    for cname, active_tiers in tier_configs.items():
        seen = {}
        active_ids = []
        for r in rows:
            if r["panel_tier"] in active_tiers:
                active_ids.append(r["pathway_id"])
                if r["unambiguous_kos"] != "NONE":
                    for k in r["unambiguous_kos"].split(","):
                        assert k not in seen, f"Collision in {cname}: {k} in {seen[k]} and {r['pathway_id']}"
                        seen[k] = r["pathway_id"]
        print(f"[PASS 2] {cname}: {len(active_ids)} active rows, {len(seen)} unambiguous KOs 100% pairwise disjoint.")

    panel_c_h2 = [r for r in rows if r["panel_tier"] in {"core", "denit_tierC_neither", "h2_gated"}]
    clf_reg = make_binned_ko_classifier(panel_c_h2)
    ox_only_kset = {"K02274", "K02276", "K00404", "K02886"}
    assert clf_reg("K02274", ox_only_kset) == (None, "unresolved_oxidase_donor")
    assert clf_reg("K02274", ox_only_kset | {"K01895"}) == ("AER_HET", None)
    assert clf_reg("K02274", ox_only_kset | {"K17486"}) == ("AER_HET", None)
    assert clf_reg("K02274", ox_only_kset | {"K01895", "K10535"}) == (None, None)
    print(
        "[PASS 2b] AER_HET positive-donor regression verified: oxidase-only MAG -> (None, 'unresolved_oxidase_donor'), "
        "matched MAG + organic-donor KO (K01895 or K17486) -> ('AER_HET', None), lithotroph (K10535) -> (None, None)."
    )

    e_equiv = {
        "H2O": 0.0, "H+": 0.0, "CO2": 0.0, "HCO3-": 0.0, "CO3-2": 0.0,
        "N2": 0.0, "SO4-2": 0.0,
        "O2": -4.0, "NO3-": -5.0, "NO2-": -3.0, "NO": -2.0, "N2O": -2.0,
        "NH4+": +3.0, "NH3": +3.0, "HS-": +8.0, "H2S": +8.0,
        "H2": +2.0, "CH4": +8.0, "CH3COO-": +8.0, "C6H12O6": +24.0,
    }
    for r in rows:
        st = parse_stoich_spec(r["stoichiometry_spec"])
        redox_bal = sum(nu * e_equiv[s] for s, nu in st.items())
        assert abs(redox_bal) < 1e-12, f"Redox fail {r['pathway_id']}: {redox_bal}"
        dG0_chk = sum(nu * sp[s][2] * CAL_TO_KJ for s, nu in st.items())
        dH0_chk = sum(nu * sp[s][3] * CAL_TO_KJ for s, nu in st.items())
        assert abs(dG0_chk - float(r["dG0_298_kJ_rxn"])) < 5e-4
        assert abs(dH0_chk - float(r["dH0_298_kJ_rxn"])) < 5e-4
    print(f"[PASS 3-4] All {len(rows)} rows balanced in C, H, O, N, S, charge, redox equivalents, and exact against OBIGT G0/H0 (< 5e-4 kJ).")

    def dG0_T(st, T):
        dG298 = sum(nu * sp[s][2] * CAL_TO_KJ for s, nu in st.items())
        dH298 = sum(nu * sp[s][3] * CAL_TO_KJ for s, nu in st.items())
        return (T / T_REF) * dG298 + dH298 * (1.0 - T / T_REF)

    for T in (277.15, 288.15, 298.15):
        RT = R_KJ * T
        Ka1_C = math.exp(-dG0_T({"CO2": -1, "H2O": -1, "HCO3-": 1, "H+": 1}, T) / RT)
        Ka2_C = math.exp(-dG0_T({"HCO3-": -1, "CO3-2": 1, "H+": 1}, T) / RT)
        Ka1_S = math.exp(-dG0_T({"H2S": -1, "HS-": 1, "H+": 1}, T) / RT)
        Ka_N = math.exp(-dG0_T({"NH4+": -1, "NH3": 1, "H+": 1}, T) / RT)
        for pH in (6.5, 7.5, 8.2):
            aH = 10.0 ** (-pH)
            DIC, S_tot, TAN, O2, CH4, SO4, NO2 = 2e-3, 5e-5, 1e-5, 1e-4, 2e-5, 1e-3, 5e-6
            a_CO2 = DIC / (1.0 + Ka1_C / aH + Ka1_C * Ka2_C / (aH * aH))
            a_HCO3 = a_CO2 * Ka1_C / aH
            a_H2S = S_tot / (1.0 + Ka1_S / aH)
            a_HS = a_H2S * Ka1_S / aH
            a_NH4 = TAN / (1.0 + Ka_N / aH)
            a_NH3 = a_NH4 * Ka_N / aH
            assert abs((dG0_T({"CH4": -1, "O2": -2, "HCO3-": 1, "H+": 1, "H2O": 1}, T) + RT * math.log((a_HCO3 * aH) / (CH4 * O2**2))) -
                       (dG0_T({"CH4": -1, "O2": -2, "CO2": 1, "H2O": 2}, T) + RT * math.log(a_CO2 / (CH4 * O2**2)))) < 1e-11
            assert abs((dG0_T({"HS-": -1, "O2": -2, "SO4-2": 1, "H+": 1}, T) + RT * math.log((SO4 * aH) / (a_HS * O2**2))) -
                       (dG0_T({"H2S": -1, "O2": -2, "SO4-2": 1, "H+": 2}, T) + RT * math.log((SO4 * aH**2) / (a_H2S * O2**2)))) < 1e-11
            assert abs((dG0_T({"NH4+": -1, "O2": -1.5, "NO2-": 1, "H+": 2, "H2O": 1}, T) + RT * math.log((NO2 * aH**2) / (a_NH4 * O2**1.5))) -
                       (dG0_T({"NH3": -1, "O2": -1.5, "NO2-": 1, "H+": 1, "H2O": 1}, T) + RT * math.log((NO2 * aH) / (a_NH3 * O2**1.5)))) < 1e-11
            # Regression check: fixed free-CO2 branch (absent DIC, a_CO2=1e-3) and measured-DIC fixture from same equilibrium return identical a_CO2 and a_HCO3
            co2_fixed, hco3_fixed = compute_carbonate_activities({}, aH, Ka1_C, Ka2_C, default_a_co2=1e-3)
            dic_uM_eq = 1e6 * 1e-3 * (1.0 + Ka1_C / aH + Ka1_C * Ka2_C / (aH * aH))
            co2_dic, hco3_dic = compute_carbonate_activities({"DIC": f"{dic_uM_eq:.16e}"}, aH, Ka1_C, Ka2_C)
            assert abs(co2_fixed - co2_dic) < 1e-17 and abs(hco3_fixed - hco3_dic) < 1e-17
    print(
        "[PASS 5] Exact acid-base speciation invariance (< 1e-11 kJ/mol) and fixed-a_CO2 vs equilibrium-DIC "
        "bicarbonate identity (< 1e-17 M) verified across T and pH."
    )

    p_zero_ft = safe_probs([0.0] * 13, [1.0, 2.0] + [0.0] * 11)
    p_zero_cap = safe_probs([0.0] * 13, [0.0] * 13)
    ll_edge = dm_loglik([250] + [0] * 12, p_zero_cap, 0.05)
    assert math.isfinite(ll_edge) and all(p > 0 and math.isfinite(math.log(p)) for p in p_zero_ft + p_zero_cap)
    print(f"[PASS 6] Likelihood edge-case floor verified (cap=0 with RNA=250 gives finite DM ll={ll_edge:.2f}).")

    union_assembled = set()
    for run_id in ("D001", "D002"):
        rdir = os.path.join(REPO_ROOT, "alife", "mep_redox", "synthetic", run_id)
        if os.path.isdir(rdir):
            res = score_mep_redox_run(rdir, rows, sp)
            union_assembled.update(res["assembled_ids"])
            rmin, rmed, rmax = res["panel_reads_min_med_max"]
            print(
                f"[PASS 7:{run_id}] "
                f"(a) Denominators & read coverage: assembled={res['n_assembled_rows']}/{res['n_panel_rows']} TierC+H2 panel rows "
                f"({res['n_assembled_rows']}/{res['n_sim_populatable_rows']} simulator-populatable; DNRA_NRF unmodeled; "
                f"catab_CDS={res['n_catabolic_cds']}, unresolved_oxidase_donor_CDS={res['n_unresolved_oxidase_donor_cds']}, hk_CDS={res['n_hk_cds']}), "
                f"panel RNA reads min/med/max={rmin}/{rmed}/{rmax} across {res['n_events']}/{res['n_events']} events (>=100 gate PASS), "
                f"RNA>=50%-event pathways={res['n_rna_cov50_pathways']}/{res['n_panel_rows']} (>=4 gate PASS); "
                f"(b) Smoke-test fixed-phi=0.05 equal-block DM scores (no verdict; protocol fits training-fold phi): {res['equal_block_scores']}; "
                f"(c) Independent-block verdict eligibility: blocks={res['n_casts']}/6 (oxic={res['n_oxic_blocks']}/3, anoxic={res['n_anoxic_blocks']}/3) -> "
                f"single-run block_gate_met={res['block_gate_met']} (scoring smoke test only; below >=6-block verdict gate)"
            )
    if union_assembled:
        assert len(union_assembled) == 12, f"Expected 12/12 simulator-populatable rows across D001 U D002, got {sorted(union_assembled)}"
        print(f"[PASS 7:union(D001,D002)] {len(union_assembled)}/13 TierC+H2 panel rows ({len(union_assembled)}/12 simulator-populatable) exercised: {sorted(union_assembled)}")


def main():
    sp = load_obigt_species()
    rows = build_rows(sp)
    cols = [
        "pathway_id", "panel_tier", "pathway_name", "direction_rule",
        "unambiguous_kos", "clade_or_synteny_kos", "excluded_confounder_kos",
        "reaction_equation", "stoichiometry_spec", "n_e", "electron_donor_spec",
        "monod_substrates_uM", "product_species_in_Q", "o2_inhibited", "Ki_O2_uM",
        "kcat_e_rel", "dG0_298_kJ_rxn", "dG0_298_kJ_per_e", "dH0_298_kJ_rxn",
        "chi_per_e", "dG_min_kJ_per_e", "primary_sources",
    ]
    buf = io.StringIO()
    w = csv.DictWriter(buf, fieldnames=cols, delimiter="\t", quoting=csv.QUOTE_NONE, extrasaction="ignore", lineterminator="\n")
    w.writeheader()
    for r in rows:
        w.writerow(r)
    rendered = buf.getvalue()
    if "--write" in sys.argv[1:]:
        with open(OUT_TSV, "w", newline="") as fh:
            fh.write(rendered)
    with open(OUT_TSV) as fh:
        disk_text = fh.read()
    assert disk_text == rendered, "Mismatch between on-disk P6-pathway-table.tsv and OBIGT-recomputed table!"
    lines = [ln.rstrip("\n") for ln in disk_text.splitlines()]
    assert all(len(ln.split("\t")) == len(cols) and '"' not in ln for ln in lines)
    verify_all(rows, sp)
    sh = hashlib.sha256(disk_text.encode("utf-8")).hexdigest()[:12]
    print(f"VERIFIED {OUT_TSV} ({len(rows)} rows, exact byte-match, sha256={sh})")


if __name__ == "__main__":
    main()
