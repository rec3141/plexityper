"""P6 field-data extractor for the approved P6-A1 Saanich cohort (protocols/P6-A1-saanich-cohort.md).

Differences from the synthetic extractor (p6_scorer.extract_run), all required by P6 / P6-A1 and
the codex-plexi run-spec review (251098919257686016):
- Panel tier from the cohort's measured analytes: N2O measured, NO and H2 not -> core +
  Tier B_N2O (DEN_NIRNOR, DEN_NOS) = 11 rows; h2_gated rows inactive.
- Capacity per CDS: cap_r,e = sum over CDS in r of overlap_bp_cds,e / d_prok,e, with d_prok,e the
  median over the 14 prokaryotic single-copy marker families of (sum overlap_bp / sum CDS length)
  over prokaryotic marker CDSs (table notes §2.2). overlap_bp comes from p6_field_count.py, the
  same assignment used for RNA.
- Binned CDS: frozen make_binned_ko_classifier with the MAG's KO set. The Tier-C nirK/norB
  exclusions (anammox carriers; AOB hao/K28504 without norB/nosZ) apply to DEN_NIRNOR only;
  DEN_NOS (K00376) follows its frozen row without an anammox exclusion.
- Unbinned CDS: unambiguous_kos only (P6 'KO and direction collisions': ambiguous KOs are not
  assigned on unbinned contigs).
- Eukaryotic CDS (Tiara eukarya / organelle) are excluded from panel and marker depth.
- RNA n_r,e: assigned fragments per CDS summed over the event's libraries (A/B replicates).
"""
import csv
import math
import os
import statistics
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import build_and_verify_p6_table as bv  # noqa: E402

ACTIVE_TIERS = {"core", "denit_tierAB_N2O", "denit_tierB_N2O_only"}
TIERB_ROWS = {"DEN_NIRNOR", "DEN_NOS"}
MODELS = ["M_cap", "M_capK", "M_capKT", "M_power"]
EPS_OBS = 1e-4


def panel_rows(all_rows):
    rows = [r for r in all_rows if r["panel_tier"] in ACTIVE_TIERS]
    assert len(rows) == 11, len(rows)
    return rows


def make_classifier(all_rows, prows):
    """Frozen binned classifier, given the rows it needs for its KO-set lookups, restricted to
    the active panel. Returns classify(ko, mag_kset or None for unbinned) -> pathway_id|None."""
    need = {r["pathway_id"]: r for r in all_rows}
    crow = list(prows) + [need[p] for p in ("HYD_UP", "MCR_MET") if p not in {x["pathway_id"] for x in prows}]
    binned = bv.make_binned_ko_classifier(crow)
    active = {r["pathway_id"] for r in prows}
    unamb = {}
    for r in prows:
        if r["unambiguous_kos"] != "NONE":
            for k in r["unambiguous_kos"].split(","):
                unamb[k] = r["pathway_id"]

    def classify(ko, kset):
        if not ko:
            return None
        if kset is None:                         # unbinned: directional markers only
            return unamb.get(ko)
        pid, _ = binned(ko, kset)
        if pid == "DEN_NIRNOR":
            # the frozen DEN_NO2N2 exclusions apply to the nir/nor KOs it pooled; the frozen
            # DEN_NOS row (K00376) carries no anammox exclusion (codex-plexi 251099827689406464)
            if "K20932" in kset:
                return None
            if ("K10535" in kset or "K28504" in kset) and not ({"K04561", "K00376"} & kset):
                return None
        return pid if pid in active else None
    return classify


def marker_depth(overlap_bp, cds_len, cds_marker_family):
    fam_num, fam_den = {}, {}
    for cid, fam in cds_marker_family.items():
        fam_num[fam] = fam_num.get(fam, 0.0) + overlap_bp.get(cid, 0.0)
        fam_den[fam] = fam_den.get(fam, 0.0) + cds_len[cid]
    d = [fam_num[f] / fam_den[f] for f in fam_den if fam_den[f] > 0]
    assert d, "no marker CDS"
    return statistics.median(d)   # true median (codex-plexi 251099827689406464)


def n2_air_sat_uM(T_K, S):
    """Weiss (1970) N2 solubility, ml/L at 1 atm moist air -> uM (22.391 L/mol)."""
    A1, A2, A3, A4, B1, B2, B3 = -172.4965, 248.4262, 143.0738, -21.7120, -0.049781, 0.025018, -0.0034861
    T = T_K / 100.0
    ln_c = A1 + A2 / T + A3 * math.log(T) + A4 * T + S * (B1 + B2 * T + B3 * T * T)
    return math.exp(ln_c) / 22.391 * 1000.0


def event_weights(prows, sp, chem):
    """chem: dict with T_K, S, rho (kg/L), pH, DIC_umolkg, DOC_uM, O2, NO3, NO2, TAN, SO4_uM,
    H2S, CH4_uM, N2O_uM, N2_uM (all uM unless noted). Returns per-row (fk, ft, yield_e)."""
    T = chem["T_K"]; RT = bv.R_KJ * T; aH = 10.0 ** (-chem["pH"])

    def dG0_T(st):
        dG298 = sum(nu * sp[s][2] * bv.CAL_TO_KJ for s, nu in st.items())
        dH298 = sum(nu * sp[s][3] * bv.CAL_TO_KJ for s, nu in st.items())
        return (T / bv.T_REF) * dG298 + dH298 * (1.0 - T / bv.T_REF)

    Ka1_C = math.exp(-dG0_T({"CO2": -1, "H2O": -1, "HCO3-": 1, "H+": 1}) / RT)
    Ka2_C = math.exp(-dG0_T({"HCO3-": -1, "CO3-2": 1, "H+": 1}) / RT)
    Ka1_S = math.exp(-dG0_T({"H2S": -1, "HS-": 1, "H+": 1}) / RT)
    Ka_N = math.exp(-dG0_T({"NH4+": -1, "NH3": 1, "H+": 1}) / RT)
    DIC_M = chem["DIC_umolkg"] * chem["rho"] * 1e-6
    a_CO2 = DIC_M / (1.0 + Ka1_C / aH + Ka1_C * Ka2_C / (aH * aH))
    a_HCO3 = a_CO2 * Ka1_C / aH
    a_H2S = chem["H2S"] * 1e-6 / (1.0 + Ka1_S / aH)
    a_HS = a_H2S * Ka1_S / aH
    a_NH4 = chem["TAN"] * 1e-6 / (1.0 + Ka_N / aH)
    act = {"H2O": 1.0, "H+": aH, "O2": chem["O2"] * 1e-6, "NO3-": chem["NO3"] * 1e-6,
           "NO2-": chem["NO2"] * 1e-6, "N2": chem["N2_uM"] * 1e-6, "N2O": chem["N2O_uM"] * 1e-6,
           "NH4+": a_NH4, "SO4-2": chem["SO4_uM"] * 1e-6, "HS-": a_HS, "CH4": chem["CH4_uM"] * 1e-6,
           "HCO3-": a_HCO3, "CH3COO-": 0.05 * chem["DOC_uM"] * 1e-6}
    conc = {"CH3COO-": 0.05 * chem["DOC_uM"], "O2": chem["O2"], "NO3-": chem["NO3"], "NO2-": chem["NO2"],
            "TAN": chem["TAN"], "SO4-2": chem["SO4_uM"], "SigmaH2S": chem["H2S"], "CH4": chem["CH4_uM"],
            "N2O": chem["N2O_uM"]}
    out = []
    for r in prows:
        fk = 1.0
        for s_mon, ks in bv.parse_monod_spec(r["monod_substrates_uM"]):
            c = conc[s_mon]; fk *= c / (ks + c)
        if r["o2_inhibited"] == "true":
            ki = float(r["Ki_O2_uM"]); fk *= ki / (ki + chem["O2"])
        st = bv.parse_stoich_spec(r["stoichiometry_spec"])
        ne = float(r["n_e"])
        lnQ = sum(nu * math.log(act[s]) for s, nu in st.items())
        dGe = (dG0_T(st) + RT * lnQ) / ne
        arg = (dGe + float(r["dG_min_kJ_per_e"])) / (float(r["chi_per_e"]) * RT)
        ft = 0.0 if arg >= 0.0 else 1.0 - math.exp(arg)
        out.append((fk, ft, max(0.0, -dGe)))
    return out


def model_probs(caps, w):
    W = {"M_cap": caps,
         "M_capK": [c * fk for c, (fk, ft, y) in zip(caps, w)],
         "M_capKT": [c * fk * ft for c, (fk, ft, y) in zip(caps, w)],
         "M_power": [c * fk * ft * y for c, (fk, ft, y) in zip(caps, w)]}
    return {m: bv.safe_probs(W[m], caps, EPS_OBS) for m in MODELS}


# ---------------------------------------------------------------- P6-A1 chemistry supply
FALLBACK_DL_UM = 0.01   # table notes §3.3


def _half_dl(v):
    """nondetect ('0' or blank-as-zero) -> 1/2 DL with 0.01 uM fallback DL."""
    return v if v > 0 else 0.5 * FALLBACK_DL_UM


def scenario_grid():
    """P6-A1 §3: 4 pH x 3 DIC x 3 DOC = 36 scenarios. pH 'lin' = 8.0 at 10 m -> 7.4 at 200 m."""
    out = []
    for ph in (7.4, 7.7, 8.0, 'lin'):
        for dic in (1900.0, 2100.0, 2300.0):
            for doc in (50.0, 75.0, 100.0):
                out.append({'pH': ph, 'DIC': dic, 'DOC': doc})
    assert len(out) == 36
    return out


CENTRAL = {'pH': 'lin', 'DIC': 2100.0, 'DOC': 75.0}


def ph_at(ph, depth_m):
    if ph != 'lin':
        return float(ph)
    d = min(max(depth_m, 10.0), 200.0)
    return 8.0 + (7.4 - 8.0) * (d - 10.0) / 190.0


def build_chem(raw, ctd, scen, sulfide_so4=0.0):
    """raw: Dryad Chemical row values (uM; CH4/N2O in nM). ctd: dict T_C, S, sigma_theta.
    sulfide_so4: coefficient for the S-SO4 sensitivity (0 primary; 0.5 approved; 1.0 diagnostic)."""
    T_K = ctd['T_C'] + 273.15
    S = ctd['S']
    rho = (1000.0 + ctd['sigma_theta']) / 1000.0
    H2S = _half_dl(raw['H2S'])
    so4 = 28.24 * S / 35.0 * rho * 1000.0 - sulfide_so4 * H2S
    return dict(T_K=T_K, S=S, rho=rho, pH=ph_at(scen['pH'], raw['depth']), DIC_umolkg=scen['DIC'],
                DOC_uM=scen['DOC'], O2=_half_dl(raw['O2']), NO3=_half_dl(raw['NO3']), NO2=_half_dl(raw['NO2']),
                TAN=_half_dl(raw['NH4']), SO4_uM=so4, H2S=H2S,
                CH4_uM=_half_dl(raw['CH4_nM'] / 1000.0), N2O_uM=_half_dl(raw['N2O_nM'] / 1000.0),
                N2_uM=n2_air_sat_uM(T_K, S))


# ---------------------------------------------------------------- event assembly
def build_events(cds, run_counts, run_event, chem_by_event, scen, all_rows, sp, sulfide_so4=0.0,
                 expected_runs=None):
    """cds: {cds_id: dict(length, mag (None if unbinned), ko, euk(bool), marker_family or None)}
    mag_kos derived from cds. run_counts: {run: {'kind': 'DNA'|'RNA', 'counts': {cds: (frag, bp)}}}.
    run_event: {run: event_id}. chem_by_event: {event: (raw, ctd, block)}.
    Returns events for p6_scorer.Cohort: sample_id, cast_id (=block), O2_uM, N_panel, counts, probs."""
    if expected_runs is not None:
        # production: supplied runs must equal the 93-run manifest exactly, so a missing A/B
        # library cannot yield a partial event (codex-plexi 251100396642549760)
        got, exp = set(run_counts), set(expected_runs)
        if got != exp or set(run_event) != exp:
            raise ValueError(f'run set mismatch: missing {sorted(exp - got)[:5]}, extra {sorted(got - exp)[:5]}')
    prows = panel_rows(all_rows)
    pidx = {r['pathway_id']: i for i, r in enumerate(prows)}
    classify = make_classifier(all_rows, prows)
    mag_kos = {}
    for c in cds.values():
        if c['mag'] and not c['euk'] and c['ko']:
            mag_kos.setdefault(c['mag'], set()).add(c['ko'])
    cds_pid = {}
    for cid, c in cds.items():
        if c['euk']:
            continue
        pid = classify(c['ko'], mag_kos.get(c['mag'], set()) if c['mag'] else None)
        if pid is not None:
            cds_pid[cid] = pid
    markers = {cid: c['marker_family'] for cid, c in cds.items() if c.get('marker_family') and not c['euk']}
    lengths = {cid: c['length'] for cid, c in cds.items()}
    dna_bp, rna_n = {}, {}
    for run, rc in run_counts.items():
        ev = run_event[run]
        tgt = dna_bp if rc['kind'] == 'DNA' else rna_n
        acc = tgt.setdefault(ev, {})
        for cid, (frag, bp) in rc['counts'].items():
            acc[cid] = acc.get(cid, 0) + (bp if rc['kind'] == 'DNA' else frag)
    events = []
    for ev in sorted(chem_by_event):
        if ev not in dna_bp or ev not in rna_n:
            continue
        raw, ctd, block = chem_by_event[ev]
        d = marker_depth(dna_bp[ev], lengths, markers)
        caps = [0.0] * len(prows); counts = [0] * len(prows)
        for cid, pid in cds_pid.items():
            caps[pidx[pid]] += dna_bp[ev].get(cid, 0.0) / d
            counts[pidx[pid]] += rna_n[ev].get(cid, 0)
        chem = build_chem(raw, ctd, scen, sulfide_so4)
        w = event_weights(prows, sp, chem)
        events.append(dict(sample_id=ev, cast_id=block, O2_uM=chem['O2'], N_panel=sum(counts),
                           counts=counts, probs=model_probs(caps, w), caps=caps, d_prok=d))
    return events
