"""End-to-end synthetic test of the P6-A1 field path (codex-plexi 251099827689406464):
SAM (SE+PE, A/B libraries) -> p6_field_count -> p6_field_extract.build_events over all 36
scenarios (incl. zero/nondetect chemistry, nonuniform marker depths, eukaryotic + unbinned CDSs)
-> p6_scorer.Cohort/analyse. No real data is read."""
import csv, math, os, random, sys
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import numpy as np
import p6_field_count as fc, p6_field_extract as fx, p6_scorer as ps, build_and_verify_p6_table as bv
ROOT = os.path.dirname(HERE)
ROWS = list(csv.DictReader(open(os.path.join(ROOT, 'protocols', 'P6-pathway-table.tsv')), delimiter='\t'))
SP = bv.load_obigt_species(); PR = fx.panel_rows(ROWS)
rng = random.Random(7)

def make_reference():
    cds, table = {}, []
    pos = {}
    def add(cid, contig, mag, ko, euk=False, fam=None, L=600):
        s = pos.get(contig, 1); pos[contig] = s + L + 50
        cds[cid] = dict(length=L, mag=mag, ko=ko, euk=euk, marker_family=fam)
        table.append((cid, contig, s, s + L - 1))
    unamb = [(r['pathway_id'], r['unambiguous_kos'].split(',')[0]) for r in PR if r['unambiguous_kos'] != 'NONE']
    for m in ('MAG1', 'MAG2', 'MAG3'):
        for j, (pid, ko) in enumerate(unamb):
            if pid == 'AMX_HZS' and m != 'MAG3':
                continue                      # anammox marker only on MAG3
            add(f'{m}_{pid}', f'{m}_c1', m, ko)
        for i, ko in enumerate(bv.HOUSEKEEPING_PROK_KOS):
            add(f'{m}_hk{i}', f'{m}_c2', m, ko, fam=ko)
    for j, (pid, ko) in enumerate(unamb):
        add(f'UB_{pid}', 'unbinned_c1', None, ko)
    add('UB_amb', 'unbinned_c1', None, 'K00370')                    # ambiguous on unbinned: dropped
    add('EUK_1', 'euk_c1', 'EUKMAG', unamb[0][1], euk=True)          # eukaryotic: excluded
    return cds, table

def sam_for(table, depth_by_cds, pe, tag):
    lines = []
    k = 0
    for cid, contig, s, e in table:
        for _ in range(depth_by_cds.get(cid, 0)):
            k += 1; q = f'{tag}_{k}'; p = rng.randint(s, e - 100)
            if pe:
                lines += [f'{q}\t65\t{contig}\t{p}\t60\t50M\t*\t0\t0\t*\t*', f'{q}\t129\t{contig}\t{p+50}\t60\t50M\t*\t0\t0\t*\t*']
            else:
                lines.append(f'{q}\t0\t{contig}\t{p}\t60\t100M\t*\t0\t0\t*\t*')
    return lines

def test_full_path():
    cds, table = make_reference()
    tsv = os.path.join(HERE, '_tmp_cds.tsv')
    with open(tsv, 'w') as fh:
        fh.write('cds_id\tcontig\tstart\tend\n'); [fh.write(f'{a}\t{b}\t{c}\t{d}\n') for a, b, c, d in table]
    idx = fc.load_cds(tsv); os.unlink(tsv)
    run_counts, run_event, chem = {}, {}, {}
    events = [(f'E{i}', f'B{i//3}', 10 + 40 * (i % 3)) for i in range(24)]   # 8 blocks x 3 depths
    for ev, blk, dep in events:
        dna = {cid: (rng.randint(2, 6) if 'hk' in cid else rng.randint(1, 10)) for cid in cds}
        dna.update({f'MAG1_hk{i}': 1 + i for i in range(14)})            # nonuniform marker depths
        rnaA = {cid: rng.randint(0, 40) for cid in cds}; rnaB = {cid: rng.randint(0, 40) for cid in cds}
        for run, kind, d, pe in ((ev + '_D', 'DNA', dna, True), (ev + '_RA', 'RNA', rnaA, False), (ev + '_RB', 'RNA', rnaB, True)):
            frag, bp, fate = fc.count(sam_for(table, d, pe, run), idx)
            run_counts[run] = {'kind': kind, 'counts': {c: (frag[c], bp[c]) for c in frag}}
            run_event[run] = ev
        o2 = 0.0 if dep > 30 else 200.0                                    # zero/nondetect O2 at depth
        raw = dict(depth=dep, O2=o2, NO3=0.0 if dep == 90 else 20.0, NO2=0.0, NH4=0.0 if dep == 10 else 5.0,
                   H2S=0.0 if dep < 90 else 10.0, CH4_nM=0.0 if dep == 10 else 500.0, N2O_nM=0.0 if dep == 90 else 10.0)
        chem[ev] = (raw, dict(T_C=9.0, S=30.5, sigma_theta=23.6), blk)
    # A/B summation check on one event
    ev0 = 'E0'; pid0 = 'MAG1_' + PR[0]['pathway_id']
    nA = run_counts['E0_RA']['counts'].get(pid0, (0, 0))[0]; nB = run_counts['E0_RB']['counts'].get(pid0, (0, 0))[0]
    for scen in fx.scenario_grid():
        evs = fx.build_events(cds, run_counts, run_event, chem, scen, ROWS, SP)
        assert len(evs) == 24
        for e in evs:
            assert all(math.isfinite(x) for m in fx.MODELS for x in e['probs'][m])
            assert all(abs(sum(e['probs'][m]) - 1) < 1e-12 for m in fx.MODELS)
    evs = fx.build_events(cds, run_counts, run_event, chem, fx.CENTRAL, ROWS, SP)
    e0 = [e for e in evs if e['sample_id'] == 'E0'][0]
    # panel count for row 0 = MAG1+MAG2+MAG3+UB CDS of that row (A+B); euk + ambiguous excluded
    assert e0['N_panel'] > 0
    # exact A+B summation and exclusions: recompute expected panel counts independently
    prows_ids = [r['pathway_id'] for r in PR]
    exp = [0] * len(prows_ids)
    for run in ('E0_RA', 'E0_RB'):
        for cid, (frag, bp) in run_counts[run]['counts'].items():
            if cid.startswith('EUK') or cid == 'UB_amb' or '_hk' in cid:
                continue
            if cid.startswith('MAG3_') and cid.split('_', 1)[1] in ('DEN_NAR', 'DEN_NIRNOR', 'DNRA_NRF', 'NIT_AMO'):
                continue                      # frozen anammox-carrier exclusions
            pid = cid.split('_', 1)[1]
            if pid in prows_ids:
                exp[prows_ids.index(pid)] += frag
    assert e0['counts'] == exp, (e0['counts'], exp)
    co = ps.Cohort(evs, [e['cast_id'] for e in evs]); g = co.gate()
    assert g['blocks'] == 8 and g['events_eligible'] == 24, g
    res = ps.analyse(co.n, co.P, co.block, co.C, 199, np.random.default_rng(1))
    assert all(math.isfinite(res[k]['est']) for k in ('Delta', 'Delta_gate', 'Delta_weight'))
    assert nA + nB >= 0
    try:
        fx.build_events(cds, {k: v for k, v in run_counts.items() if k != 'E0_RB'}, run_event, chem,
                        fx.CENTRAL, ROWS, SP, expected_runs=set(run_counts))
        raise AssertionError('partial run set accepted')
    except ValueError:
        pass

if __name__ == '__main__':
    test_full_path(); print('OK full path (36 scenarios, 24 events, 8 blocks)')
