"""Fixtures for scratch/p6_field_extract.py (P6-A1 field path; codex-plexi 251098919257686016)."""
import csv, math, os, sys
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import p6_field_extract as fx
import build_and_verify_p6_table as bv
ROOT = os.path.dirname(HERE)
ROWS = list(csv.DictReader(open(os.path.join(ROOT, 'protocols', 'P6-pathway-table.tsv')), delimiter='\t'))
SP = bv.load_obigt_species()
PR = fx.panel_rows(ROWS)

def test_panel_tierB():
    ids = [r['pathway_id'] for r in PR]
    assert len(ids) == 11 and 'DEN_NIRNOR' in ids and 'DEN_NOS' in ids
    assert not {'MCR_MET', 'HYD_FERM', 'HYD_UP', 'DEN_NO2N2', 'DEN_NIR', 'DEN_NOR'} & set(ids)

def test_unbinned_only_unambiguous():
    cl = fx.make_classifier(ROWS, PR)
    assert cl('K00370', None) is None            # narG/nxrA collision: not assigned unbinned
    assert cl('K10944', None) is None            # amoA/pmoA collision
    nos = [r for r in PR if r['pathway_id'] == 'DEN_NOS'][0]['unambiguous_kos'].split(',')[0]
    assert cl(nos, None) == 'DEN_NOS'

def test_binned_tierB_anammox_excluded():
    cl = fx.make_classifier(ROWS, PR)
    nos = [r for r in PR if r['pathway_id'] == 'DEN_NOS'][0]['unambiguous_kos'].split(',')[0]
    assert cl(nos, {nos}) == 'DEN_NOS'
    assert cl(nos, {nos, 'K20932'}) is None

def test_marker_depth_and_capacity_ratio():
    # 14 families, depth 2.0 each -> d_prok 2.0; CDS a: 300 bp overlap, b: 900 bp -> cap 150 : 450
    fam = {f'm{i}': f'F{i}' for i in range(14)}; ln = {f'm{i}': 1000 for i in range(14)}
    ov = {f'm{i}': 2000.0 for i in range(14)}
    d = fx.marker_depth(ov, ln, fam); assert d == 2.0
    assert (300 / d, 900 / d) == (150.0, 450.0)

def test_n2_air_sat():
    v = fx.n2_air_sat_uM(283.15, 30.0)
    assert 450 < v < 650, v

def test_weights_finite_saanich_like():
    chem = dict(T_K=282.0, S=31.0, rho=1.0245, pH=7.6, DIC_umolkg=2100, DOC_uM=75, O2=2.0, NO3=5.0,
                NO2=0.1, TAN=3.0, SO4_uM=28.24 * 31 / 35 * 1.0245 * 1000, H2S=2.0, CH4_uM=0.3,
                N2O_uM=0.01, N2_uM=fx.n2_air_sat_uM(282.0, 31.0))
    w = fx.event_weights(PR, SP, chem)
    assert len(w) == 11 and all(all(math.isfinite(x) for x in t) for t in w)
    caps = [1.0] * 11
    P = fx.model_probs(caps, w)
    assert all(abs(sum(P[m]) - 1) < 1e-12 for m in fx.MODELS)

if __name__ == '__main__':
    n = 0
    for k, v in sorted(globals().items()):
        if k.startswith('test_'):
            v(); n += 1
    print(f'OK {n} fixtures')
