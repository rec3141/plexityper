"""Run the whole multiscale_learning sweep: simulate worlds, write ground truth and synthetic
observations, run the inference tests, write result tables.

Usage:  python3 -I run_all.py [--quick]
--quick shrinks lineage counts and durations for a smoke test (results not meaningful).
"""
import os, sys, time, json, hashlib
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from msl import model as Mdl, truth as Tr, observe as Ob, analysis as An  # noqa: E402
from msl.model import M, R, L, DT  # noqa: E402

QUICK = '--quick' in sys.argv
REVISION = '--revision' in sys.argv   # rerun only truth/KO/P3/bloom/mechinterp/P1
OUT = os.path.join(HERE, 'results_quick' if QUICK else 'results')
os.makedirs(OUT, exist_ok=True)
SEED = 101
N_LIN = 300 if QUICK else 1200
OBS_Y = 6
BURN_Y = 2 if QUICK else 4
CORE_N = 200
T0 = time.time()
LOG = open(os.path.join(OUT, 'run_log_revision.txt' if '--revision' in sys.argv else 'run_log.txt'), 'w')


def log(*a):
    s = '[%6.0fs] ' % (time.time() - T0) + ' '.join(str(x) for x in a)
    print(s, flush=True); LOG.write(s + '\n'); LOG.flush()


def write_tsv(name, rows, keys=None):
    path = os.path.join(OUT, name)
    os.makedirs(os.path.dirname(path), exist_ok=True)
    if not rows:
        open(path, 'w').close(); return
    if keys is None:
        keys = []
        for r in rows:
            keys += [k for k in r if k not in keys]
    with open(path, 'w') as fh:
        fh.write('\t'.join(keys) + '\n')
        for r in rows:
            fh.write('\t'.join(_fmt(r.get(k, '')) for k in keys) + '\n')


def _fmt(v):
    if isinstance(v, (float, np.floating)):
        return 'nan' if not np.isfinite(v) else ('%.4g' % v)
    if isinstance(v, (list, tuple)):
        return ','.join(_fmt(x) for x in v)
    return str(v)


DESIGNS = Ob.make_designs(OBS_Y)
OFFS = Ob.union_offsets(DESIGNS)
DESIGN_ROWS = {d: np.searchsorted(OFFS, [Ob.step_offset(a, b) for a, b in v]) for d, v in DESIGNS.items()}
CAMP_STARTS = [2 * 365 + 20, 3 * 365 + 20]
PERT_DAY = 14
SNAP_DAYS = [200, 500, 800, 1100, 1400, 1700]


def run_world(name, pool, over, seed=SEED, n_steps=None, snaps=True):
    cfg = dict(n_lin=N_LIN, burn_years=BURN_Y, obs_years=OBS_Y)
    cfg.update(over)
    t = time.time()
    rec = Mdl.simulate(pool, cfg, seed=seed, rec_offsets=OFFS,
                       snap_offsets=[d * 10 + 5 for d in SNAP_DAYS] if snaps else (),
                       n_steps=n_steps)
    log(f'world {name}: simulated in {time.time() - t:.0f}s, records {len(rec["t"])}')
    rec['name'] = name
    return rec


def events_for(rec, design, station=0, arm='ctl'):
    out = []
    for (d, h), row in zip(DESIGNS[design], DESIGN_ROWS[design]):
        if row >= len(rec['t']):
            continue
        out.append(dict(event_id=f"{rec['name']}-st{station:02d}-{arm}-d{d:04d}h{h}",
                        station=station, arm=arm, day=d, hour=round(float(Ob.HOURS[h]), 1),
                        year=d // 365, row=int(row)))
    return out


# ============================================================== ground truth outputs
def truth_outputs(rec, pool, tag):
    t = time.time()
    per, edges = Tr.truth_series(rec_daily(rec), pool, rec['cfg'], every_days=14, core_n=300)
    write_tsv(f'truth/{tag}/d_cycle_per_lineage.tsv', per,
              ['t', 'lineage', 'B', 'viable', 'd_h0', 'd_h1', 'd_hi_h1', 'd_unresolved_h1', 'layer',
               'kappa_any_support', 'kappa_all_supports', 'supports_complete', 'supports_unknown',
               'n_deficient', 'deficient', 'fan_class', 'plexitrophic_h0',
               'plexitrophic_h1', 'n_required', 'n_necessary', 'n_realised', 'n_realised_not_required',
               'euk_required', 'genes_lost'])
    write_tsv(f'truth/{tag}/edges.tsv',
              [dict(t=e[0], donor=('EUK' if e[1] == Tr.EUK else e[1]), recipient=e[2], kind=e[3], value=e[4])
               for e in edges])
    log(f'truth {tag}: {len(per)} lineage-times, {len(edges)} edges in {time.time() - t:.0f}s')
    v = [p for p in per if p['viable']]
    B = np.array([p['B'] for p in v])
    def frac(key, cond=lambda x: x >= 1):
        x = np.array([cond(p[key]) for p in v], float)
        return float(x.mean()), float((x * B).sum() / B.sum())
    summ = dict(world=tag, n_lineage_times=len(per), n_viable=len(v))
    for label, key, c in [('genome_dep_d_h1>=1', 'd_h1', lambda x: x >= 1),
                          ('plexi_regulated_d_h1>=2', 'd_h1', lambda x: x >= 2),
                          ('any_dep_frozen_d_h0>=1', 'd_h0', lambda x: x >= 1),
                          ('plexi_frozen_d_h0>=2', 'd_h0', lambda x: x >= 2),
                          ('kappa_all_supports', 'kappa_all_supports', lambda x: x == 1),
                          ('kappa_any_support', 'kappa_any_support', lambda x: x == 1),
                          ('d_hi_unresolved', 'd_unresolved_h1', lambda x: x == 1),
                          ('plexitrophic_h1', 'plexitrophic_h1', lambda x: x == 1),
                          ('plexitrophic_h0', 'plexitrophic_h0', lambda x: x == 1),
                          ('euk_required', 'euk_required', lambda x: x == 1),
                          ('realised_not_required>=1', 'n_realised_not_required', lambda x: x >= 1)]:
        f, fb = frac(key, c)
        summ[label + '_frac'] = round(f, 4); summ[label + '_Bweighted'] = round(fb, 4)
    fc = np.array([p['fan_class'] for p in v]); gl = np.array([p['genes_lost'] for p in v])
    d1 = np.array([p['d_h1'] for p in v])
    for c in ('functional_fan_in', 'dilution_witnessed', 'dilution_unresolved', 'other'):
        summ[f'd>=2_{c}_frac'] = round(float((fc == c).mean()), 4)
    summ['supports_unknown_rows'] = int(sum(p['supports_unknown'] for p in v))
    summ['rows_with_d>=1'] = int((d1 >= 1).sum())
    for lo, hi, nm in ((0, 0, '0'), (1, 1, '1'), (2, 4, '2-4'), (5, 5, '5')):
        msk = (gl >= lo) & (gl <= hi)
        summ[f'genes_lost_{nm}_frac'] = round(float(msk.mean()), 4)
        summ[f'd>=2_given_genes_lost_{nm}'] = round(float((d1[msk] >= 2).mean()) if msk.any() else float('nan'), 4)
    dm = d1 >= 2
    summ['d>=2_with_m0_deficient_frac'] = round(float(np.mean(['0' in p['deficient'].split(',') for p, x in zip(v, dm) if x])) if dm.any() else float('nan'), 4)
    lay = np.array([p['layer'] for p in v])
    for l_ in ('genome', 'expression', 'none'):
        summ[f'layer_{l_}_frac'] = round(float((lay == l_).mean()), 4)
    nreal = np.array([p['n_realised'] for p in v]); nreq = np.array([p['n_required'] for p in v])
    summ['mean_n_realised'] = round(float(nreal.mean()), 3)
    summ['mean_n_required'] = round(float(nreq.mean()), 3)
    return per, edges, summ


def rec_daily(rec):
    """View of rec restricted to the daily_ideal rows (noon each day)."""
    rows = DESIGN_ROWS['daily_ideal']
    rows = rows[rows < len(rec['t'])]
    out = {k: (v[rows] if isinstance(v, np.ndarray) and v.ndim >= 1 and len(v) == len(rec['t']) else v)
           for k, v in rec.items()}
    return out


def ko_validation(rec, pool, tag, n_donor=4, n_keep=8):
    """Validate the snapshot sensitivity against re-simulated removal trajectories (20 d, same
    forcing, from full-state snapshots).
    (a) single-donor removal: donors called necessary at h1 vs not required;
    (b) keep-only-one-supplier: for recipients with d = 1 (a single rescuing supplier exists) vs
        d >= 2 (no single supplier rescues, exhaustively tested), remove the leak of every
        lineage except the recipient and its best single supplier, and switch off eukaryote
        metabolite release unless EUK is that supplier. Prediction: d = 1 recipients persist,
        d >= 2 recipients decline relative to control.
    """
    rows, rows_b = [], []
    cfg = rec['cfg']; env = rec['env']; bs = rec['burn_steps']
    n_h = int(20 / DT)
    for k_step, snap in sorted(rec['snaps'].items()):
        T = env[1][k_step]
        r1 = {kk: np.array([snap[kk]]) for kk in ('B', 'S', 'P', 'Esyn', 'Eup', 'Ecat', 'F')}
        r1.update(T=np.array([T]), hc=np.array([snap['hc']]), A=np.array([snap['A']]),
                  Z=np.array([snap['Z']]), idx=rec['idx'], t=np.array([0.0]))
        recips = np.argsort(-snap['B'])[:150]
        res, g = Tr.necessity_at(r1, 0, pool, cfg, recips)
        nec = {}
        for r in res:
            if r['viable'] and r['minsets_h1'] and r['d_lo_h1'] >= 1:
                for dnr in set.intersection(*[set(ms) for ms in r['minsets_h1']]):
                    if dnr != Tr.EUK:
                        nec.setdefault(dnr, set()).add(r['i'])
        donors = sorted(nec, key=lambda d: -len(nec[d]))[:n_donor]
        sim = lambda **kw: Mdl.simulate(pool, kw.pop('cfg', cfg), seed=rec['seed'], idx=rec['idx'],
                                        start_state=snap, t_start_step=k_step + 1, n_steps=n_h,
                                        env=env, record=True, rec_offsets=[k_step + n_h - bs], **kw)
        ctrl = sim()
        viable = [r['i'] for r in res if r['viable']]
        for dnr in donors:
            ko = sim(knockout=dnr)
            ratio = ko['B'][-1] / np.maximum(ctrl['B'][-1], 1e-12)
            for i in viable:
                if i == dnr:
                    continue
                cls = 'necessary_h1' if i in nec.get(dnr, ()) else 'not_required'
                rows.append(dict(test='single_donor_removal', day=(k_step - bs) * DT,
                                 donor=int(rec['idx'][dnr]), recipient=int(rec['idx'][i]),
                                 static_class=cls, ko_ratio_20d=float(ratio[i])))
        # (b) keep-only-one-supplier
        d1 = [r for r in res if r['viable'] and r['d_lo_h1'] == 1]
        d2 = [r for r in res if r['viable'] and r['d_lo_h1'] >= 2]
        N = len(rec['idx'])
        for cls, group in (('d=1', d1[:n_keep]), ('d>=2', d2[:n_keep])):
            for r in group:
                i = r['i']
                if cls == 'd=1':
                    keep = r['minsets_h1'][0][0]
                else:   # best single supplier by share of what i cannot cover (all failed alone)
                    keep = None
                    st_, g_, inflow, cons, euk = Tr.fluxes_at(r1, 0, pool, cfg)
                    sh = inflow / np.maximum(inflow.sum(0) + euk, 1e-15)
                    sc = sh.sum(1); sc[i] = -1
                    keep = int(np.argmax(sc))
                    if (euk / np.maximum(inflow.sum(0) + euk, 1e-15)).sum() > sc[keep]:
                        keep = Tr.EUK
                kill = np.array([x for x in range(N) if x != i and x != keep])
                cfg_k = dict(cfg) if keep == Tr.EUK else dict(cfg, relA=(0.0,) * M)
                ko = sim(cfg=cfg_k, knockout=kill)
                rows_b.append(dict(test='keep_only_one_supplier', day=(k_step - bs) * DT,
                                   recipient=int(rec['idx'][i]), static_class=cls,
                                   kept=('EUK' if keep == Tr.EUK else int(rec['idx'][keep])),
                                   ko_ratio_20d=float(ko['B'][-1][i] / max(ctrl['B'][-1][i], 1e-12))))
    write_tsv(f'truth/{tag}/knockout_validation.tsv', rows + rows_b,
              ['test', 'day', 'donor', 'recipient', 'static_class', 'kept', 'ko_ratio_20d'])
    summ = []
    for test, R_ in (('single_donor_removal', rows), ('keep_only_one_supplier', rows_b)):
        for cls in sorted(set(r['static_class'] for r in R_)):
            x = np.array([r['ko_ratio_20d'] for r in R_ if r['static_class'] == cls])
            summ.append(dict(world=tag, test=test, static_class=cls, n_pairs=len(x),
                             median_ko_ratio=float(np.median(x)),
                             frac_ratio_below_0_5=float((x < 0.5).mean()),
                             frac_ratio_below_0_8=float((x < 0.8).mean())))
    return summ


# ============================================================== observations
def observe_design(rec, pool, mags, design, rng, station=0, arm='ctl'):
    ev = events_for(rec, design, station, arm)
    rows = [e['row'] for e in ev]
    obs = Ob.observe_events(rec, rows, pool, rec['cfg'], mags, rng)
    obs['t'] = np.array([e['day'] + e['hour'] / 24 for e in ev])
    obs['events'] = ev
    return obs


# ============================================================== (A) layer series
def layer_series(obs, mags, top):
    """Observed series per layer for the `top` MAG indices. Returns {layer: list of 1-d arrays}."""
    tot = obs['mag_reads'].sum(1) + obs['unbinned_reads'] + obs['euk_reads']
    S = {}
    S['population'] = [np.log((obs['mag_reads'][:, j] + 0.5) / tot) for j in top]
    rm = obs['rna_mag']; hk = rm[:, :, 0] + 0.5
    ex = []
    for j in top:
        ex.append(np.log((rm[:, j, 1 + 2 * M] + 0.5) / hk[:, j]))       # cat0 (diel resource)
        ex.append(np.log((rm[:, j, 1 + M] + 0.5) / hk[:, j]))           # upt0
        for m in range(M):
            if mags[j]['present'].get('bio%d' % m):
                ex.append(np.log((rm[:, j, 1 + m] + 0.5) / hk[:, j])); break
    S['expression'] = ex
    S['genome'] = [obs['gene_cov_ratio'][:, j, m] for j in top for m in range(M)
                   if mags[j]['present'].get('bio%d' % m)]
    S['nucleotide'] = [obs['snv_freq'][:, j, s] for j in top for s in range(L)]
    S['eukaryote'] = [np.log(obs['chl'] + 1e-3)]
    return S


def truth_layer_series(recd, locs):
    S = {}
    S['population'] = [np.log(recd['B'][:, i].astype(float) + 1e-12) for i in locs]
    S['expression'] = [recd['Ecat'][:, i, 0].astype(float) for i in locs] + \
                      [recd['Eup'][:, i, 0].astype(float) for i in locs]
    S['genome'] = [recd['F'][:, i, m].astype(float) for i in locs for m in range(M)]
    S['nucleotide'] = [recd['Psnv'][:, i, s].astype(float) for i in locs for s in range(L)]
    S['eukaryote'] = [np.log(recd['A'].astype(float) + 1e-6)]
    return S


def clean(t, y, min_n=12):
    ok = np.isfinite(y)
    if ok.sum() < min_n:
        return None
    t2, y2 = t[ok], y[ok]
    if np.std(y2) < 1e-9:
        return None
    return t2, y2


def timescale_table(world, design, t, S, rng, max_series=25, null=True):
    rows = []
    span = t[-1] - t[0]
    dmed = float(np.median(np.diff(t)))
    for layer, series in S.items():
        tau, pv, fr = [], [], {'fast': [], 'mid': [], 'slow': []}
        picked = [s for s in series if clean(t, s) is not None]
        if len(picked) > max_series:
            sel = rng.choice(len(picked), max_series, replace=False)
            picked = [picked[k] for k in sel]
        for y in picked:
            tt, yy = clean(t, y)
            yd = An.deseason(tt, yy)
            if len(tt) > 400:   # long regular series: noise-robust ratio estimator r(2)/r(1)
                yc = yd - yd.mean()
                r1 = np.mean(yc[1:] * yc[:-1]); r2 = np.mean(yc[2:] * yc[:-2])
                phi = r2 / r1 if r1 > 0 else np.nan
                tau.append(-dmed / np.log(phi) if 0 < phi < 1 else (np.inf if phi >= 1 else dmed * 0.3))
                p = np.nan
            elif null:
                fit, p = An.ar1_null_test(tt, yd, rng)
                tau.append(fit['tau'])
            else:
                tau.append(An.fit_ou_noise(tt, yd)['tau']); p = np.nan
            pv.append(p)
            bf = An.band_fractions(tt, yy[None, :])
            for b in fr:
                fr[b].append(float(bf[b][0]))
        if not picked:
            continue
        tau = np.array(tau, float); pv = np.array(pv, float)
        rows.append(dict(world=world, design=design, layer=layer, n_series=len(picked),
                         n_times=len(t), span_d=round(span, 1), median_dt_d=round(dmed, 2),
                         tau_hat_median=float(np.nanmedian(np.where(np.isinf(tau), 1e5, tau))),
                         tau_hat_q25=float(np.nanpercentile(np.where(np.isinf(tau), 1e5, tau), 25)),
                         tau_hat_q75=float(np.nanpercentile(np.where(np.isinf(tau), 1e5, tau), 75)),
                         frac_tau_at_bounds=float(np.mean((tau < 0.5 * dmed) | (tau > span))),
                         ar1_null_reject_frac=float(np.mean(pv < 0.05)) if np.isfinite(pv).any() else np.nan,
                         power_fast_lt30d=float(np.mean(fr['fast'])), power_mid_30_180d=float(np.mean(fr['mid'])),
                         power_slow_gt180d=float(np.mean(fr['slow'])), _tau=tau))
    # adjacent-layer ordering AUC: P(tau_slower > tau_faster)
    order = ['expression', 'population', 'genome', 'nucleotide']
    byl = {r['layer']: r['_tau'] for r in rows}
    for r in rows:
        r['auc_vs_next_faster'] = np.nan
        if r['layer'] in order and order.index(r['layer']) > 0:
            prev = order[order.index(r['layer']) - 1]
            if prev in byl:
                a, b = np.where(np.isinf(r['_tau']), 1e5, r['_tau']), np.where(np.isinf(byl[prev]), 1e5, byl[prev])
                r['auc_vs_next_faster'] = float(np.mean(a[:, None] > b[None, :]) + 0.5 * np.mean(a[:, None] == b[None, :]))
    for r in rows:
        r.pop('_tau')
    return rows


def euk_dominance(world, recd, locs):
    """Daily noise-free truth: extra R^2 of lagged algal biomass over season + temperature."""
    t = recd['t']; w = 2 * np.pi * t / 365
    base = np.column_stack([np.ones_like(t), np.sin(w), np.cos(w), np.sin(2 * w), np.cos(2 * w), recd['T']])
    lA = np.log(recd['A'].astype(float) + 1e-6)
    lags = [0, 3, 7, 14, 21]
    lagged = np.column_stack([np.roll(lA, k) for k in lags])
    ok = np.arange(len(t)) >= max(lags)
    S = truth_layer_series(recd, locs)
    rows = []
    for layer, series in S.items():
        if layer == 'eukaryote':
            continue
        gain, r2b = [], []
        for y in series:
            y = y[ok]
            if y.std() < 1e-9:
                continue
            X0 = base[ok]; X1 = np.column_stack([X0, lagged[ok]])
            def r2(X):
                b, *_ = np.linalg.lstsq(X, y, rcond=None); e = y - X @ b
                return 1 - e.var() / y.var()
            a0, a1 = r2(X0), r2(X1)
            gain.append(a1 - a0); r2b.append(a0)
        bf = An.band_fractions(t, np.array([s for s in series if s.std() > 1e-9])[:, :])
        rows.append(dict(world=world, layer=layer, n_series=len(gain),
                         r2_season_T=float(np.median(r2b)) if gain else np.nan,
                         r2_gain_lagged_algae=float(np.median(gain)) if gain else np.nan,
                         r2_gain_q75=float(np.percentile(gain, 75)) if gain else np.nan,
                         power_fast=float(bf['fast'].mean()), power_mid=float(bf['mid'].mean()),
                         power_slow=float(bf['slow'].mean()),
                         algal_share_of_organic_C_supply=float(recd['euk_doc'].mean() / (recd['euk_doc'].mean() + recd['cfg']['D'] * recd['Sin'].sum(1).mean())),
                         algal_biomass_cv=float(recd['A'].std() / (recd['A'].mean() + 1e-12))))
    return rows


# ============================================================== (C) aliasing and closure
def aliasing_rows(world, rec, pool, mags, rng, top):
    rows = []
    for design in ('weekly_2y_noon', 'weekly_2y', 'weekly_2y_drift'):
        obs = observe_design(rec, pool, mags, design, rng)
        t = obs['t']
        rm = obs['rna_mag']
        tot_rna = rm.sum((1, 2)) + obs['rna_unbinned'] + obs['rna_euk']
        Ys = {'per_housekeeping': np.array([np.log((rm[:, j, 1 + 2 * M] + 0.5) / (rm[:, j, 0] + 0.5)) for j in top]),
              'per_total_RNA': np.array([np.log((rm[:, j, 1 + 2 * M] + 0.5) / tot_rna) for j in top]),
              'prokaryote_share_of_RNA': np.log((rm.sum((1, 2)) + obs['rna_unbinned']) / tot_rna)[None, :]}
        for norm, Y in Ys.items():
            rows.append(_alias_row(world, design, norm, t, Y, obs))
    return rows


def _alias_row(world, design, norm, t, Y, obs):
    Y = np.array([An.deseason(t, y) for y in Y])
    freqs = An.freq_grid(t, 80)
    P = An.ls_power(t, Y, freqs).mean(0)
    per = 1 / freqs
    alias_band = (per > 28) & (per < 42)            # weekly_2y_drift alias of diel: ~35 d
    hours = np.array([e['hour'] for e in obs['events']])
    # fraction of variance explained by sampling hour (truth: diel forcing)
    X = np.column_stack([np.ones_like(hours), np.sin(2 * np.pi * hours / 24), np.cos(2 * np.pi * hours / 24)])
    r2h = []
    for y in Y:
        b, *_ = np.linalg.lstsq(X, y, rcond=None); r2h.append(1 - (y - X @ b).var() / y.var())
    return dict(world=world, design=design, normalisation=norm, n=len(t),
                peak_period_d=float(per[np.argmax(P)]),
                power_frac_28_42d=float(P[alias_band].sum()),
                r2_by_sampling_hour=float(np.median(r2h)),
                tau_hat_median=float(np.median([An.fit_ou_noise(t, y)['tau'] for y in Y[:15]])))


def closure_rows(world, recd, obs, mags, top):
    """Relative (closed) vs absolute abundance: spurious correlations and spectra."""
    loc = np.array([mags[j]['local'] for j in top])
    absB = np.log(recd['B'][:, loc].astype(float) + 1e-12)
    tot = obs['mag_reads'].sum(1) + obs['unbinned_reads'] + obs['euk_reads']
    rel = np.log((obs['mag_reads'][:, top] + 0.5) / tot[:, None])
    relp = np.log((obs['mag_reads'][:, top] + 0.5) / (obs['mag_reads'].sum(1) + obs['unbinned_reads'])[:, None])
    t = recd['t']
    out = []
    for nm, X in (('rel_incl_euk_reads', rel), ('rel_prokaryote_only', relp)):
        r_self = [np.corrcoef(absB[:, k], X[:, k])[0, 1] for k in range(len(top))]
        Ca = np.corrcoef(absB.T); Cr = np.corrcoef(X.T)
        iu = np.triu_indices(len(top), 1)
        dc = Cr[iu] - Ca[iu]
        euk_frac = obs['euk_reads'] / tot
        r_euk = [np.corrcoef(euk_frac, X[:, k])[0, 1] for k in range(len(top))]
        r_euk_abs = [np.corrcoef(euk_frac, absB[:, k])[0, 1] for k in range(len(top))]
        bfa = An.band_fractions(t, absB.T); bfr = An.band_fractions(t, X.T)
        out.append(dict(world=world, measure=nm, n_taxa=len(top),
                        median_r_rel_vs_abs=float(np.median(r_self)),
                        mean_pairwise_corr_abs=float(Ca[iu].mean()), mean_pairwise_corr_rel=float(Cr[iu].mean()),
                        frac_pairs_sign_flip=float(np.mean(np.sign(Cr[iu]) != np.sign(Ca[iu]))),
                        frac_pairs_spurious_neg=float(np.mean((Cr[iu] < -0.3) & (Ca[iu] > -0.1))),
                        median_r_with_euk_read_frac_rel=float(np.median(r_euk)),
                        median_r_with_euk_read_frac_abs=float(np.median(r_euk_abs)),
                        slow_power_abs=float(bfa['slow'].mean()), slow_power_rel=float(bfr['slow'].mean()),
                        fast_power_abs=float(bfa['fast'].mean()), fast_power_rel=float(bfr['fast'].mean())))
    return out


# ============================================================== (D) P3 + core
def required_pairs(edges, t_lo, t_hi):
    req = set(); donors_of = {}
    for (t, d, r, kind, v) in edges:
        if kind == 'required' and t_lo <= t < t_hi and d != Tr.EUK:
            req.add((d, r)); donors_of.setdefault(r, set()).add(d)
    return req, donors_of


def p3_block(world, rec, pool, mags, edges, rng):
    rows, corerows = [], []
    idx = rec['idx']
    for design, tr_years, te_years in (('weekly_2y', [0], [1]), ('monthly_4y_gaps', [1, 2, 3], [4]),
                                       ('daily_ideal', [0, 1, 2], [3, 4, 5])):
        obs = observe_design(rec, pool, mags, design, rng)
        ev = obs['events']; years = np.array([e['year'] for e in ev]); tday = np.array([e['day'] for e in ev])
        tr = np.nonzero(np.isin(years, tr_years))[0]; te = np.nonzero(np.isin(years, te_years))[0]
        counts = obs['amp']; ntax = counts.shape[1] - 1
        env = An.env_matrix(tday, obs['env_T'], obs['env_S'], obs['chl'])
        rk = An.core_ranking(counts[tr], ntax)
        core = np.sort(np.argsort(-rk['read_share'])[:CORE_N])
        occ = An.seasonal_occupancy(counts[tr], tday[tr], ntax)
        corerows.append(dict(world=world, design=design, n_train=len(tr), core_n=CORE_N,
                             core_read_share=float(rk['read_share'][core].sum()),
                             core_detect_freq_median=float(np.median(rk['detect_freq'][core])),
                             core_seasonal_occ_median=float(np.median(occ[core])),
                             core_frac_occ_lt_0_5=float(np.mean(occ[core] < 0.5)),
                             noncore_detected_taxa=int(((counts[tr][:, :ntax] > 0).any(0)).sum() - CORE_N)))
        t_lo, t_hi = min(te_years) * 365, (max(te_years) + 1) * 365
        req, donors_of = required_pairs(edges, t_lo, t_hi)
        for universe_name, uni in (('all', np.arange(ntax)), ('core200', core)):
            res = An.p3_test(counts, tday, env, tr, te, uni, rng)
            for r in res:
                tg = int(idx[r['target']]); pr = [int(idx[p]) for p in r['preds']]
                direct = any((p, tg) in req or (tg, p) in req for p in pr)
                shared = any(len(donors_of.get(tg, set()) & donors_of.get(p, set())) > 0 for p in pr)
                rows.append(dict(world=world, design=design, universe=universe_name,
                                 target=tg, n_preds=len(pr), test_prev=round(r['test_prev'], 3),
                                 gain_me=round(r['gain_me'][0], 4), gain_me_lo=round(r['gain_me'][1], 4),
                                 gain_me_hi=round(r['gain_me'][2], 4),
                                 gain_mn=round(r['gain_mn'][0], 4), gain_mn_lo=round(r['gain_mn'][1], 4),
                                 decision=r['decision'], truth_direct_required_edge=int(direct),
                                 truth_shared_required_donor=int(shared),
                                 n_train=r['n_train'], n_test=r['n_test']))
    return rows, corerows


def p3_block_v2(world, rec, recd, pool, mags, edges, rng):
    """P3 PILOT ALIGNED WITH protocols/P3-plexitype-recurrence.md where feasible: rotated masked
    targets, delta = 0.02, 20 matched control sets (10 for the 2190-event daily design), Mn as the
    mean of log scores, background gradient from V minus (C u N_k) with the target masked from
    the depth denominator, shrinkage (lambda 5) only on added taxon slopes, Bernoulli detection with
    log depth, cluster bootstrap over 28-d event blocks. NOT the frozen P3: no JSDM conditional,
    discovery by residual-correlation cliques.

    Ground truth for the estimand: the SAME pipeline run on noise-free daily truth (target = the
    expected detection probability 1 - exp(-expected reads), scored as expected log score; predictors = true relative biomass with the target masked; true
    environment) in the same train/test years ('oracle out-of-year transfer gain'). A required
    exchange edge is reported separately as mechanism, not as the estimand.
    """
    rows, corerows = [], []
    idx = rec['idx']
    Bd = recd['B'].astype(float); Brel = Bd / Bd.sum(1, keepdims=True)
    tdd = recd['t']; yrd = (tdd // 365).astype(int)
    Yor = 1.0 - np.exp(-Brel * 4e4)     # expected detection probability at the amplicon depth (soft label)
    w_ = 2 * np.pi * tdd / 365
    Eor = np.column_stack([np.sin(w_), np.cos(w_), recd['T'], np.log(recd['S'].astype(float) + 1e-6),
                           np.log(recd['A'].astype(float) + 1e-3)])
    def feat_or(cols, target):
        den = 1.0 - (Brel[:, target] if target is not None else 0.0)
        return np.log(Brel[:, cols] / np.asarray(den)[..., None] + 1e-7)
    for design, tr_years, te_years in (('weekly_2y', [0], [1]), ('monthly_4y_gaps', [1, 2, 3], [4]),
                                       ('daily_ideal', [0, 1, 2], [3, 4, 5])):
        obs = observe_design(rec, pool, mags, design, rng)
        ev = obs['events']; years = np.array([e['year'] for e in ev]); tday = np.array([e['day'] for e in ev])
        tr = np.nonzero(np.isin(years, tr_years))[0]; te = np.nonzero(np.isin(years, te_years))[0]
        counts = obs['amp'][:, :-1]; ntax = counts.shape[1]
        depth = obs['amp'].sum(1).astype(float)
        Eb = np.column_stack([An.env_matrix(tday, obs['env_T'], obs['env_S'], obs['chl']), np.log(depth)])
        Y = (counts > 0).astype(float)
        def feat_obs(cols, target):
            dm = depth - (counts[:, target] if target is not None else 0)
            return np.log((counts[:, cols] + 0.5) / dm[:, None])
        rk = An.core_ranking(obs['amp'][tr], ntax)
        core = np.sort(np.argsort(-rk['read_share'])[:CORE_N])
        occ = An.seasonal_occupancy(obs['amp'][tr], tday[tr], ntax)
        corerows.append(dict(world=world, design=design, n_train=len(tr), core_n=CORE_N,
                             core_read_share=float(rk['read_share'][core].sum()),
                             core_detect_freq_median=float(np.median(rk['detect_freq'][core])),
                             core_seasonal_occ_median=float(np.median(occ[core])),
                             core_frac_occ_lt_0_5=float(np.mean(occ[core] < 0.5)),
                             noncore_detected_taxa=int(((counts[tr] > 0).any(0)).sum() - CORE_N)))
        t_lo, t_hi = min(te_years) * 365, (max(te_years) + 1) * 365
        req, donors_of = required_pairs(edges, t_lo, t_hi)
        trd = np.nonzero(np.isin(yrd, tr_years))[0]; ted = np.nonzero(np.isin(yrd, te_years))[0]
        n_ctrl = 10 if design == 'daily_ideal' else 20
        for universe_name, uni in (('all', np.arange(ntax)), ('core200', core)):
            elig, groups = An.p3_discover(Y, lambda cols: np.log((counts[:, cols] + 0.5) / depth[:, None]),
                                          Eb, tr, uni)
            if not groups:
                rows.append(dict(world=world, design=design, universe=universe_name, verdict='N_no_clique'))
                continue
            prev_or = Yor[trd].mean(0)
            elig_or = np.nonzero((prev_or >= 0.2) & (prev_or <= 0.8))[0]
            for grp in groups:
                r = An.p3_score_group(Y, feat_obs, Eb, tday, tr, te, elig, grp, uni, rng, n_ctrl=n_ctrl)
                ver = An.p3_verdict(r)
                elig_o = np.union1d(elig_or, grp)
                ro = An.p3_score_group(Yor, feat_or, Eor, tdd, trd, ted, elig_o, grp, np.arange(ntax), rng,
                                       n_ctrl=5, nboot=300)
                ver_o = An.p3_verdict(ro)
                lin = [int(idx[x]) for x in grp]
                direct = any((a, b) in req for a in lin for b in lin if a != b)
                row = dict(world=world, design=design, universe=universe_name, members=lin, n_members=len(lin),
                           verdict=ver, n_train=len(tr), n_test=len(te),
                           oracle_verdict=ver_o,
                           truth_direct_required_edge_within_group=int(direct))
                for nm, rr in (('', r), ('oracle_', ro)):
                    if rr is not None:
                        row.update({f'{nm}gain_me': round(rr['gain_me'][0], 4), f'{nm}gain_me_lo': round(rr['gain_me'][1], 4),
                                    f'{nm}gain_me_hi': round(rr['gain_me'][2], 4), f'{nm}gain_mn': round(rr['gain_mn'][0], 4),
                                    f'{nm}gain_mn_lo': round(rr['gain_mn'][1], 4), f'{nm}gain_mn_hi': round(rr['gain_mn'][2], 4)})
                rows.append(row)
    return rows, corerows


def bloom_rows(world, recd, per, core_lineages):
    B = recd['B'].astype(float)
    rel = B / B.sum(1, keepdims=True)
    bloom = (rel.max(0) > 1e-3) & (rel.min(0) < 1e-5)
    idx = recd['idx']
    dep = {}
    for p in per:
        if p['viable']:
            dep.setdefault(p['lineage'], []).append((p['d_h1'] >= 1, p['d_h1'] >= 1 and p['genes_lost'] >= 1,
                                                     p['d_h1'] >= 2 and p['fan_class'] == 'functional_fan_in'))
    out = []
    in_core = np.isin(idx, list(core_lineages))
    for nm, msk in (('bloomers_outside_core', bloom & ~in_core), ('bloomers_in_core', bloom & in_core),
                    ('core_nonbloomers', ~bloom & in_core)):
        L_ = idx[msk]
        evald = [l for l in L_ if l in dep]
        f = lambda k: float(np.mean([any(x[k] for x in dep[l]) for l in evald])) if evald else float('nan')
        out.append(dict(world=world, group=nm, n=int(msk.sum()), n_with_truth=len(evald),
                        frac_with_genome_dependency=f(0), frac_dependency_and_genes_lost_ge1=f(1),
                        frac_functional_fan_in=f(2)))
    return out


# ============================================================== (F) mechinterp
def mechinterp(world, rec, recd, pool, mags, edges, obs_daily, rng, edge_source='own'):
    rows = []
    idx = rec['idx']
    t = obs_daily['t']
    yr = (t // 365).astype(int)
    counts = obs_daily['amp'][:, :-1]
    rk = An.core_ranking(counts[yr == 0], counts.shape[1])
    core = np.sort(np.argsort(-rk['read_share'])[:CORE_N])
    prop = counts[:, core] / counts[:, core].sum(1, keepdims=True)
    clr = np.log(counts[:, core] + 0.5); clr -= clr.mean(1, keepdims=True)
    # true modules over core taxa
    gl = idx[core]
    Fm = recd['F'][:, core, :].astype(float).mean(0)
    mods = {f'guild_resource{r}': pool['guild'][gl] == r for r in range(R)}
    q = np.percentile(pool['Topt'][gl], [33, 67])
    mods.update({'season_cold': pool['Topt'][gl] < q[0], 'season_warm': pool['Topt'][gl] > q[1]})
    for m in range(M):
        mods[f'auxotroph_m{m}'] = Fm[:, m] < 0.5
    lA = np.log(recd['A'].astype(float) + 1e-6)
    rA = np.array([np.corrcoef(lA, np.log(recd['B'][:, i].astype(float) + 1e-12))[0, 1] for i in core])
    mods['algal_bloom_followers'] = rA > 0.4
    k = 8
    U, s, Vt = np.linalg.svd(clr - clr.mean(0), full_matrices=False)
    W, H = An.nmf(prop, k, rng)
    for method, Ld in (('PCA_clr', Vt[:k]), ('NMF_prop', H)):
        al = An.module_alignment(Ld, mods, rng, nnull=300)
        for name, v in al.items():
            rows.append(dict(world=world, probe=method, target=name, kind='feature',
                             score=round(v['score'], 3), null_mean=round(v['null_mean'], 3),
                             percentile_vs_random=round(v['pct'], 3), size=v['size']))
    # circuits: per-capita growth probe vs required edges (donor -> recipient), core pairs only
    req = {}
    for (tt, d, r, kind, v) in edges:
        if kind == 'required' and d != Tr.EUK:
            req[(d, r)] = req.get((d, r), 0) + 1
    pos = {l: n for n, l in enumerate(gl)}
    truthM = np.zeros((len(core), len(core)), bool)
    for (d, r), c in req.items():
        if d in pos and r in pos and c >= 2:
            truthM[pos[d], pos[r]] = True
    off = ~np.eye(len(core), dtype=bool)
    for design, step in (('daily_ideal', 1), ('weekly_from_daily', 7)):
        sel = np.arange(0, len(t), step)
        X = clr[sel]; dX = np.diff(X, axis=0)
        env = np.column_stack([obs_daily['env_T'][sel], np.log(obs_daily['chl'][sel] + 1e-3),
                               np.sin(2 * np.pi * t[sel] / 365), np.cos(2 * np.pi * t[sel] / 365)])[:-1]
        Xp = np.column_stack([X[:-1], env])
        Wt = An.ridge(Xp, dX, lam=10.0 * Xp.shape[0] ** 0.5)[:len(core)]   # (donor j, recipient i)
        C = np.corrcoef(X.T)
        lab = truthM[off]
        if lab.sum() == 0:
            continue
        for pname, sc in (('ridge_growth', Wt), ('abs_ridge_growth', np.abs(Wt)), ('clr_correlation', C)):
            a_obs = An.auroc(sc[off], lab)
            # donor-preserving null: permute recipients within each donor row
            nulls = []
            for _ in range(30):
                Tp = np.array([rng.permutation(row) for row in truthM])
                nulls.append(An.auroc(sc[off], Tp[off]))
            rows.append(dict(world=world, probe=f'{pname}_{design}', target=f'required_edges[{edge_source}]',
                             kind='circuit', score=round(a_obs, 3), null_mean=round(float(np.nanmean(nulls)), 3),
                             percentile_vs_random=round(float(np.mean(np.array(nulls) < a_obs)), 3),
                             size=int(lab.sum())))
    # hidden-state probe: decode unobserved pools P_m from RNA vs DNA features
    top = np.argsort(-obs_daily['mag_cov'].mean(0))[:60]
    rm = obs_daily['rna_mag']
    Xr = np.column_stack([np.log((rm[:, j, 1 + M + m] + 0.5) / (rm[:, j, 0] + 0.5)) for j in top for m in range(M)])
    tot = obs_daily['mag_reads'].sum(1) + obs_daily['unbinned_reads'] + obs_daily['euk_reads']
    Xd = np.log((obs_daily['mag_reads'][:, top] + 0.5) / tot[:, None])
    envb = np.column_stack([obs_daily['env_T'], np.log(obs_daily['chl'] + 1e-3), np.log(obs_daily['env_S']),
                            np.sin(2 * np.pi * t / 365), np.cos(2 * np.pi * t / 365)])
    for m in range(M):
        y = np.log(recd['P'][:, m].astype(float) + 1e-6)
        for nm, X in (('env_only', envb), ('env+RNA_uptake_genes', np.column_stack([envb, Xr])),
                      ('env+DNA_abundance', np.column_stack([envb, Xd]))):
            sse = sst = 0.0
            for yy in np.unique(yr):                   # leave-one-year-out
                trn = yr != yy; tst = yr == yy
                mu, sd = X[trn].mean(0), X[trn].std(0) + 1e-9
                Z = (X - mu) / sd
                w = np.linalg.solve(Z[trn].T @ Z[trn] + 200 * np.eye(Z.shape[1]), Z[trn].T @ (y[trn] - y[trn].mean()))
                pred = Z[tst] @ w + y[trn].mean()
                sse += np.sum((y[tst] - pred) ** 2); sst += np.sum((y[tst] - y[trn].mean()) ** 2)
            rows.append(dict(world=world, probe=f'linear_probe_{nm}', target=f'hidden_pool_m{m}',
                             kind='state', score=round(float(1 - sse / sst), 3), null_mean=0.0,
                             percentile_vs_random=np.nan, size=len(y)))
    return rows


# ============================================================== (B) perturbation
def perturbation_rows(rec_c, rec_t, pool, mags, rng):
    rows = []
    loc = np.array([m['local'] for m in mags])
    dc, dt_ = rec_daily(rec_c), rec_daily(rec_t)
    n = min(len(dc['t']), len(dt_['t']))
    top = loc[np.argsort(-dc['B'][:n][:, loc].mean(0))[:60]]
    for ci, cs in enumerate(CAMP_STARTS):
        t0 = cs + PERT_DAY
        # truth, daily, noise-free
        win = (dc['t'][:n] >= cs - 30) & (dc['t'][:n] < cs + 120)
        tt = dc['t'][:n][win]; pre = tt < t0
        feats = {
            'expression': lambda d: np.concatenate([d['Ecat'][:n][win][:, top, :].reshape(win.sum(), -1),
                                                    d['Esyn'][:n][win][:, top, :].reshape(win.sum(), -1)], 1),
            'population': lambda d: np.log(d['B'][:n][win][:, top].astype(float) + 1e-12),
            'genome': lambda d: d['F'][:n][win][:, top, :].reshape(win.sum(), -1),
            'nucleotide': lambda d: d['Psnv'][:n][win][:, top, :].reshape(win.sum(), -1),
        }
        for layer, f in feats.items():
            Xc, Xt = f(dc).astype(float), f(dt_).astype(float)
            dv = An.divergence(Xt, Xc, pre)
            t10, t50, tpk = An.response_times(tt, dv, t0, pre, noise_floor=False)
            rows.append(dict(campaign=ci + 1, source='truth_daily', layer=layer, t10_d=t10, t50_d=t50,
                             tpeak_d=tpk, peak_div=float(dv[tt >= t0].max())))
    # observed, campaign design (two arms, same event schedule, independent sequencing noise)
    oc = observe_design(rec_c, pool, mags, 'two_campaigns', rng, arm='ctl')
    ot = observe_design(rec_t, pool, mags, 'two_campaigns', rng, arm='trt')
    jtop = np.argsort(-oc['mag_cov'].mean(0))[:60]
    for ci, cs in enumerate(CAMP_STARTS):
        t0 = cs + PERT_DAY
        msk = (oc['t'] >= cs) & (oc['t'] < cs + 80)
        tt = oc['t'][msk]; pre = tt < t0
        def lay(o):
            tot = o['mag_reads'].sum(1) + o['unbinned_reads'] + o['euk_reads']
            rm = o['rna_mag']
            ex = np.log((rm[:, jtop, :] + 0.5) / (rm[:, jtop, :1] + 0.5))[:, :, 1:].reshape(len(rm), -1)
            return dict(expression=ex, population=np.log((o['mag_reads'][:, jtop] + 0.5) / tot[:, None]),
                        genome=o['gene_cov_ratio'][:, jtop, :].reshape(len(rm), -1),
                        nucleotide=o['snv_freq'][:, jtop, :].reshape(len(rm), -1))
        Lc, Lt = lay(oc), lay(ot)
        for layer in Lc:
            Xc, Xt = Lc[layer][msk], Lt[layer][msk]
            dv = An.divergence(Xt, Xc, pre)
            t10, t50, tpk = An.response_times(tt, dv, t0, pre, noise_floor=True)
            rows.append(dict(campaign=ci + 1, source='observed_campaign', layer=layer, t10_d=t10, t50_d=t50,
                             tpeak_d=tpk, peak_div=float(dv[tt >= t0].max()) if (tt >= t0).any() else np.nan,
                             pre_div_mean=float(dv[pre].mean()), pre_div_sd=float(dv[pre].std())))
    return rows, oc, ot


# ============================================================== (E) P1 stations
def p1_block(rng, n_st=20):
    """P1 PROXY STRESS TEST (not the reviewed P1 protocol): LOSO over stations, margin 0.1 nat,
    protocol replication gate reported (>= 20 stations and >= 2 dominant genera at >= 5 stations),
    dominant-genus transfer folds, bootstrap slope at the same penalty as the predictive model."""
    rows, details = [], []
    pool = Mdl.make_pool(tail='powerlaw')
    days = [15 + 30 * k for k in range(12)]
    offs = np.array([d * 10 + 5 for d in days])
    for world, vit in (('P1_mixed_vitamer_off', False), ('P1_pseudo_vitamer_on', True)):
        data = []
        for st in range(n_st):
            sd = 5000 + st
            dT = float(np.random.default_rng(sd).uniform(-3, 3))
            cfg = dict(n_lin=200 if QUICK else 400, burn_years=2 if QUICK else 3, obs_years=1,
                       vitamer=vit, station_dT=dT)
            rec = Mdl.simulate(pool, cfg, seed=sd, rec_offsets=offs)
            rec['name'] = world
            mags, _ = Ob.build_catalogue(rec, pool, np.random.default_rng(sd + 1))
            obs = Ob.observe_events(rec, np.arange(len(offs)), pool, rec['cfg'], mags, np.random.default_rng(sd + 2))
            tot = obs['mag_reads'].sum(1) + obs['unbinned_reads'] + obs['euk_reads']
            rel = obs['mag_reads'] / tot[:, None]
            has0 = np.array([m['present'].get('bio0', False) for m in mags])
            comp = np.array([m['completeness_est'] for m in mags])
            lin = np.array([m['lineage'] for m in mags])
            # true usable-fraction of supply (for the truth label) : form-A share of producers
            Bm = rec['B'].astype(float)
            prod = rec['F'][:, :, 0] > 0.5
            fb = pool['form_b'][rec['idx']]
            usable = (Bm * prod * ~fb).sum(1) / np.maximum((Bm * prod).sum(1), 1e-12) if vit else np.ones(len(offs))
            rk = rel.mean(0)
            core = np.argsort(-rk)[:CORE_N]
            for j in np.nonzero(~has0 & (comp >= 0.7))[0]:
                hk = obs['rna_mag'][:, j, 0]
                y = np.log((obs['rna_mag'][:, j, 1 + M] + 0.5) / (hk + 0.5))
                for e in range(len(offs)):
                    if hk[e] < 20:          # gate per recipient-event (Carbon 249001105845305344)
                        continue
                    data.append(dict(station=st, event=e, recipient=int(lin[j]), y=y[e],
                                     T=obs['env_T'][e], chl=np.log(obs['chl'][e] + 1e-3), doy=days[e], rec_ab=np.log(rel[e, j] + 1e-7),
                                     guild_all=np.log(rel[e, has0].sum() + 1e-7),
                                     guild_core=np.log(rel[e, np.intersect1d(np.nonzero(has0)[0], core)].sum() + 1e-7),
                                     sup=dict(zip(lin[has0], rel[e, has0])), usable=usable[e],
                                     dominant=(int(lin[has0][np.argmax(rel[e, has0])])
                                               if has0.any() and rel[e, has0].max() > 0.5 * rel[e, has0].sum() else -1),
                                     oracle_true_pool=np.log(float(rec['P'][e, 0]) + 1e-6),
                                     true_eup=np.log(float(rec['Eup'][e, mags[j]['local'], 0]) + 1e-6)))
            log(f'P1 {world} station {st}: dT={dT:.1f}, recipients rows so far {len(data)}')
        details.append((world, data))
        stations = np.array([d['station'] for d in data])
        if len(np.unique(stations)) < 3:
            rows.append(dict(world=world, guild_set='n/a', n_stations=len(np.unique(stations)), n_obs=len(data),
                             role_recurrence_decision='not_run_too_few_stations',
                             identity_decision='not_run_too_few_stations', truth_identity_matters=int(vit)))
            continue
        y = np.array([d['y'] for d in data])
        recs_ = sorted(set(d['recipient'] for d in data))
        FE = np.array([[float(d['recipient'] == r_) for r_ in recs_] for d in data])   # recipient intercepts
        envX = np.column_stack([[d['T'] for d in data], np.sin(2 * np.pi * np.array([d['doy'] for d in data]) / 365),
                                np.cos(2 * np.pi * np.array([d['doy'] for d in data]) / 365),
                                [d['rec_ab'] for d in data], [d['chl'] for d in data], FE])
        for guild_kind in ('guild_all', 'guild_core', 'oracle_true_pool'):
            G = np.array([d[guild_kind] for d in data])[:, None]
            lpd = {k: np.zeros(len(y)) for k in ('Me', 'Mg', 'Mgi')}
            slopes = []
            for s in np.unique(stations):
                trn = stations != s; tst = ~trn
                # identity features: top-6 suppliers by training mean share among producers
                tot_sup = {}
                for d, tr_ in zip(data, trn):
                    if tr_:
                        for l_, v in d['sup'].items():
                            tot_sup[l_] = tot_sup.get(l_, 0) + v
                topsup = sorted(tot_sup, key=lambda l_: -tot_sup[l_])[:3]
                I = np.array([[np.log(d['sup'].get(l_, 0) + 1e-7) for l_ in topsup] for d in data])
                for k, X in (('Me', envX), ('Mg', np.column_stack([envX, G])), ('Mgi', np.column_stack([envX, G, I]))):
                    lp, w = An.gauss_fit_lpd(X[trn], y[trn], X[tst], y[tst], lam=5.0)
                    lpd[k][tst] = lp
                    if k == 'Mg':
                        slopes.append(w[-1])
            g_guild = An.block_boot_ci(lpd['Mg'] - lpd['Me'], stations, rng)
            g_id = An.block_boot_ci(lpd['Mgi'] - lpd['Mg'], stations, rng)
            sl = np.array(slopes)
            # station-bootstrap of the guild slope (full-data refits on resampled stations)
            bs = []
            us = np.unique(stations)
            XG = np.column_stack([envX, G])
            for _ in range(300):
                pick = rng.choice(us, len(us))
                ii = np.concatenate([np.nonzero(stations == p)[0] for p in pick])
                _, w = An.gauss_fit_lpd(XG[ii], y[ii], XG[:1], y[:1], lam=5.0)
                bs.append(w[-1])
            sl_hi = float(np.percentile(bs, 97.5))
            delta = 0.1
            # replication gate (protocol): >= 20 stations; >= 2 dominant genera at >= 5 stations each
            dom_st = {}
            for d in data:
                if d['dominant'] >= 0:
                    dom_st.setdefault(d['dominant'], set()).add(d['station'])
            n_dom5 = sum(1 for v in dom_st.values() if len(v) >= 5)
            gate = (len(us) >= 20) and (n_dom5 >= 2)
            # dominant-genus transfer: hold out stations where genus dominates any sample
            n_tr_tested = n_tr_pass = 0
            for gen, sts in dom_st.items():
                if len(sts) < 2:
                    continue
                hold = np.isin(stations, list(sts)); trn = ~hold
                if len(np.unique(stations[trn])) < 3:
                    continue
                lpE, _ = An.gauss_fit_lpd(envX[trn], y[trn], envX[hold], y[hold], lam=5.0)
                lpG, _ = An.gauss_fit_lpd(XG[trn], y[trn], XG[hold], y[hold], lam=5.0)
                gtr = An.block_boot_ci(lpG - lpE, stations[hold], rng)
                n_tr_tested += 1; n_tr_pass += int(gtr[1] > delta)
            role = 'support' if (g_guild[1] > delta and -delta < g_id[1] and g_id[2] < delta and sl_hi < 0) else (
                'against' if g_guild[2] < delta else 'inconclusive')
            ident = 'support' if g_id[1] > delta else ('against' if g_id[2] < delta else 'inconclusive')
            usable_var = float(np.std([d['usable'] for d in data]))
            tp = np.array([d['oracle_true_pool'] for d in data]); te_ = np.array([d['true_eup'] for d in data])
            r_truth = float(np.corrcoef(tp, te_)[0, 1])                    # truth: supply -> uptake RNA
            r_guild_pool = float(np.corrcoef(G[:, 0], tp)[0, 1])            # observed guild vs true pool
            r_obs_true = float(np.corrcoef(y, te_)[0, 1]); sd_ratio = float(np.std(y) / (np.std(te_) + 1e-12))
            if not gate:
                role = role + '|gate_failed->inconclusive_under_protocol'
            rows.append(dict(study='P1_proxy_stress_test', world=world, guild_set=guild_kind, n_stations=len(us), n_obs=len(y),
                             gate_pass=int(gate), n_dominant_genera_ge5_stations=n_dom5,
                             transfer_genera_tested=n_tr_tested, transfer_genera_gain_lo_gt_margin=n_tr_pass,
                             guild_gain=round(g_guild[0], 4), guild_gain_lo=round(g_guild[1], 4), guild_gain_hi=round(g_guild[2], 4),
                             identity_gain=round(g_id[0], 4), identity_gain_lo=round(g_id[1], 4), identity_gain_hi=round(g_id[2], 4),
                             guild_slope_bootstrap_hi=round(sl_hi, 4), role_recurrence_decision=role,
                             identity_decision=ident,
                             truth_identity_matters=int(vit), truth_sd_usable_fraction=round(usable_var, 3),
                             truth_r_pool_vs_uptake_expr=round(r_truth, 3), r_guild_vs_true_pool=round(r_guild_pool, 3),
                             r_observed_y_vs_true_uptake=round(r_obs_true, 3), sd_ratio_observed_to_true=round(sd_ratio, 2)))
    return rows


# ============================================================== main
def main():
    rng = np.random.default_rng(SEED)
    if '--p1-only' in sys.argv:
        write_tsv(f"{NAMES['p1']}.tsv", p1_block(rng, n_st=4 if QUICK else 20))
        log('P1 only: done'); return
    summaries = dict(truth=[], ko=[], ts=[], euk=[], alias=[], closure=[], p3=[], core=[], bloom=[],
                     mech=[], pert=[], p1=[])
    pools = dict(powerlaw=Mdl.make_pool(tail='powerlaw'), lognormal=Mdl.make_pool(tail='lognormal'))
    worlds = [('W_powerlaw', 'powerlaw', {}), ('W_lognormal', 'lognormal', {}),
              ('W_flat_euk', 'powerlaw', dict(euk='flat')), ('W_no_exchange', 'powerlaw', dict(leak=False))]
    A_mean = None
    edges_pl = None
    for wname, tail, over in worlds:
        pool = pools[tail]
        if over.get('euk') == 'flat':
            over = dict(over, A_flat=A_mean[0], Z_flat=A_mean[1])
        rec = run_world(wname, pool, over, snaps=(wname == 'W_powerlaw'))
        recd = rec_daily(rec)
        if wname == 'W_powerlaw':
            A_mean = (float(recd['A'].mean()), float(recd['Z'].mean()))
            log('mean algal biomass', A_mean)
        per, edges, ts = truth_outputs(rec, pool, wname)
        summaries['truth'].append(ts)
        if wname == 'W_powerlaw':
            summaries['ko'] += ko_validation(rec, pool, wname)
            log('knockout validation done')
            # layer time series truth (daily, top 300 lineages by mean biomass)
            top300 = np.argsort(-recd['B'].mean(0))[:300]
            np.savez_compressed(os.path.join(OUT, 'truth', wname, 'layer_series_daily_top300.npz'),
                                t=recd['t'], T=recd['T'], S=recd['S'], P=recd['P'], A=recd['A'], Z=recd['Z'],
                                lineage=rec['idx'][top300], B=recd['B'][:, top300], F=recd['F'][:, top300],
                                Esyn=recd['Esyn'][:, top300], Eup=recd['Eup'][:, top300],
                                Ecat=recd['Ecat'][:, top300], Psnv=recd['Psnv'][:, top300],
                                Sin=recd['Sin'], storm=recd['storm'])
        rng_w = np.random.default_rng([SEED, len(wname)])
        mags, cds = Ob.build_catalogue(recd, pool, rng_w)
        log(f'{wname}: catalogue {len(mags)} MAGs, {len(cds)} CDS')
        # observations + timescale tables per design
        obs_daily = observe_design(rec, pool, mags, 'daily_ideal', rng_w)
        top = np.argsort(-obs_daily['mag_cov'].mean(0))[:40]
        if REVISION and wname == 'W_powerlaw':     # re-emit observer tables in the aligned format
            for design in ('weekly_2y', 'monthly_4y_gaps', 'two_campaigns', 'subdaily_burst'):
                o_ = obs_daily if design == 'daily_ideal' else observe_design(rec, pool, mags, design, rng_w)
                Ob.write_tables(os.path.join(OUT, 'observations', wname, design), design,
                                [{k: v for k, v in e.items() if k != 'row'} for e in o_['events']], o_, mags, cds)
        for design in ([] if REVISION else DESIGNS):
            if design in ('weekly_2y_noon', 'weekly_2y_drift') and wname != 'W_powerlaw':
                continue
            if wname in ('W_flat_euk', 'W_no_exchange') and design not in ('daily_ideal', 'weekly_2y', 'monthly_4y_gaps'):
                continue
            obs = obs_daily if design == 'daily_ideal' else observe_design(rec, pool, mags, design, rng_w)
            if wname == 'W_powerlaw' and design in ('weekly_2y', 'monthly_4y_gaps', 'two_campaigns', 'subdaily_burst'):
                Ob.write_tables(os.path.join(OUT, 'observations', wname, design), design,
                                [{k: v for k, v in e.items() if k != 'row'} for e in obs['events']],
                                obs, mags, cds)
            S = layer_series(obs, mags, top)
            summaries['ts'] += timescale_table(wname, design, obs['t'], S, rng_w,
                                               null=(wname in ('W_powerlaw', 'W_lognormal')))
            log(f'{wname} {design}: timescales done (n={len(obs["t"])})')
        locs = np.array([mags[j]['local'] for j in top])
        if not REVISION:
          summaries['ts'] += timescale_table(wname, 'truth_daily_noisefree', recd['t'],
                                           truth_layer_series(recd, locs), rng_w, null=False)
        if wname in ('W_powerlaw', 'W_flat_euk', 'W_lognormal') and not REVISION:
            summaries['euk'] += euk_dominance(wname, recd, locs)
        if wname == 'W_powerlaw' and not REVISION:
            summaries['alias'] += aliasing_rows(wname, rec, pool, mags, rng_w, top[:20])
            summaries['closure'] += closure_rows(wname, recd, obs_daily, mags, top)
            log('aliasing/closure done')
        if wname != 'W_flat_euk':
            p3rows, corerows = p3_block_v2(wname, rec, recd, pool, mags, edges, rng_w)
            summaries['p3'] += p3rows; summaries['core'] += corerows
            obs_w = observe_design(rec, pool, mags, 'weekly_2y', rng_w)
            tr0 = np.array([e['year'] == 0 for e in obs_w['events']])
            rk = An.core_ranking(obs_w['amp'][tr0], obs_w['amp'].shape[1] - 1)
            core_l = set(rec['idx'][np.argsort(-rk['read_share'])[:CORE_N]].tolist())
            summaries['bloom'] += bloom_rows(wname, recd, per, core_l)
            log(f'{wname}: P3 done ({len(p3rows)} groups)')
            if wname == 'W_powerlaw':
                edges_pl = edges
            if wname == 'W_no_exchange':   # negative control: no exchange, but score against W_powerlaw's edges
                summaries['mech'] += mechinterp(wname, rec, recd, pool, mags, edges_pl, obs_daily, rng_w,
                                                edge_source='W_powerlaw_edges_same_lineages')
            else:
                summaries['mech'] += mechinterp(wname, rec, recd, pool, mags, edges, obs_daily, rng_w)
            log(f'{wname}: mechinterp done')
        if wname == 'W_powerlaw' and not REVISION:
            pert = [dict(t0=BURN_Y * 365 + cs + PERT_DAY, dur=10.0, tox=0.6, hc_in=3.0) for cs in CAMP_STARTS]
            n_steps = int((BURN_Y * 365 + CAMP_STARTS[-1] + 130) / DT)
            rec_t = run_world('W_powerlaw_treated', pool, dict(perturb=pert), n_steps=n_steps, snaps=False)
            rec_t['name'] = 'W_powerlaw'
            prow, oc, ot = perturbation_rows(rec, rec_t, pool, mags, rng_w)
            summaries['pert'] += prow
            ev_t = [{k: v for k, v in e.items() if k != 'row'} for e in ot['events']]
            Ob.write_tables(os.path.join(OUT, 'observations', wname, 'two_campaigns_treated'), 'two_campaigns',
                            ev_t, ot, mags, cds)
            del rec_t
            log('perturbation done')
        del rec, recd, obs_daily
        write_all(summaries)
    summaries['p1'] = p1_block(rng, n_st=4 if QUICK else 20)
    log('P1 done')
    write_all(summaries)
    log('ALL DONE')


NAMES = dict(truth='truth_summary', ko='knockout_validation_summary', ts='timescales_by_design_layer',
             euk='eukaryote_bloom_dominance', alias='aliasing', closure='closure', p3='p3_groups',
             core='core_ranking', bloom='bloomers_vs_core', mech='mechinterp_probes',
             pert='perturbation_propagation', p1='p1_song_vs_singer')


def write_all(summaries):
    names = NAMES
    for k, v in summaries.items():
        if REVISION and not v:
            continue              # keep tables from the original sweep that were not rerun
        write_tsv(f'{names[k]}.tsv', v)
    # provenance
    prov = []
    for f in ['run_all.py', 'msl/model.py', 'msl/truth.py', 'msl/observe.py', 'msl/analysis.py']:
        prov.append(dict(file=f, sha256_12=hashlib.sha256(open(os.path.join(HERE, f), 'rb').read()).hexdigest()[:12]))
    write_tsv('provenance_revision.tsv' if REVISION else 'provenance.tsv', prov)
    json.dump(dict(revision=REVISION, seed=SEED, n_lin=N_LIN, burn_years=BURN_Y, obs_years=OBS_Y, quick=QUICK,
                   runtime_s=round(time.time() - T0), numpy=np.__version__,
                   designs={k: len(v) for k, v in DESIGNS.items()}),
              open(os.path.join(OUT, 'run_meta_revision.json' if REVISION else 'run_meta.json'), 'w'), indent=1)


if __name__ == '__main__':
    main()
