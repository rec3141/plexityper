"""Ground truth from recorded states: realised transfers, necessity by counterfactual removal,
d (minimal partner number), cycles, and per-lineage dependence layer.

Definitions (stated, graded where possible):
- Realised transfer j->i on metabolite m at time t: flux = share_j,m(t) * consumption_i,m(t), with
  share_j,m = j's leak inflow into the usable pool of m / total usable inflow (well-mixed pool).
  'EUK' (phototrophic eukaryotes, grazing/lysis release) is a donor node like any lineage.
- Necessity (counterfactual removal, same medium/endpoint/horizon): recipient i is evaluated
  only if viable now (endpoint: net growth mu_i - D - mort >= 0 under the current medium: the
  recorded resources, temperature and pools). For a donor subset S, every usable pool is rescaled
  to the inflow coming from S plus i itself (quasi-steady pool linear in inflow). Two horizons:
    h0 'frozen'  : expression held at its recorded value (minutes; no regulatory response),
    h1 'regulated': expression re-set to its regulatory targets under the counterfactual pools
                   (hours); genome (F) fixed in both.
  d_h(i) = min |S| such that i stays viable; 0 = needs no partner (singletons tested
  exhaustively; see necessity_at). Plexitrophic iff d >= 2 or kappa = 1 on every support.
  Dependence layer: 'genome' if d_h1 >= 1 (needed even after regulation), 'expression' if
  d_h0 >= 1 but d_h1 == 0 (facultative: removable by derepressing own synthesis).
"""
import itertools
import numpy as np
from .model import growth, pool_fluxes, M, R

EUK = -1   # donor id for the eukaryote compartment


def state_at(rec, j, idx):
    return dict(idx=idx, B=rec['B'][j].astype(float), S=rec['S'][j].astype(float),
                P=rec['P'][j].astype(float), Esyn=rec['Esyn'][j].astype(float),
                Eup=rec['Eup'][j].astype(float), Ecat=rec['Ecat'][j].astype(float),
                F=rec['F'][j].astype(float), hc=float(rec['hc'][j]))


def euk_release(rec, j, cfg):
    """Metabolite release by eukaryotes (grazing/lysis) at record j, per pool m (usable forms)."""
    if cfg['euk'] == 'none':
        return np.zeros(M)
    A = rec['A'][j] if cfg['euk'] == 'dynamic' else cfg['A_flat']
    Z = rec['Z'][j] if cfg['euk'] == 'dynamic' else cfg['Z_flat']
    graz = cfg['gmax'] * Z * A / (A + cfg['Kg'])
    return np.array(cfg['relA']) * graz


def fluxes_at(rec, j, pool, cfg):
    idx = rec['idx']
    st = state_at(rec, j, idx)
    g = growth(st, pool, cfg, float(rec['T'][j]))
    inflow, cons = pool_fluxes(st, g, pool, cfg)
    return st, g, inflow[:, :M], cons, euk_release(rec, j, cfg)


def _mu_one(i, st, g, pool, cfg, Pnew, T, horizon):
    """Growth of lineage i under counterfactual usable pools Pnew (K x M). Vectorised over K."""
    idx = st['idx']; gi = idx[i]
    lam = g['lam'][i]; Tfac = g['Tfac'][i]; F = st['F'][i]
    if horizon == 'h0':
        Esyn = np.broadcast_to(st['Esyn'][i], Pnew.shape); Eup = np.broadcast_to(st['Eup'][i], Pnew.shape)
    else:
        if cfg['regulation'] and pool['regulated'][gi]:
            Krep = pool['Krep'][gi] * cfg['KP']
            Esyn = cfg['basal'] + (1 - cfg['basal']) / (1 + (Pnew / Krep)**2)
        else:
            Esyn = np.ones_like(Pnew)
        Eup = 0.05 + 0.95 * (1 - 0.8 * F * Esyn) / (1 + Pnew / (cfg['Krs_mult'] * cfg['KP']))
    y = F * Esyn * cfg['sig'] * Tfac
    upt = Eup * cfg['vup'] * Pnew / (cfg['KP'] + Pnew)
    supply = (1 - lam) * y + upt
    cm = np.asarray(cfg['cost_mult'])
    muE = cfg['Y'] * g['J'][i] - cfg['csyn'] * (y * cm).sum(1) - cfg['cgene'] * (F * cm).sum() - cfg['maint']
    return np.minimum(muE, supply.min(1))


_SUBSETS = {}


def _subsets(C, kmax=4):
    key = (C, kmax)
    if key not in _SUBSETS:
        rows = [()]
        for k in range(1, min(kmax, C) + 1):
            rows += list(itertools.combinations(range(C), k))
        mask = np.zeros((len(rows), C), bool)
        for r, comb in enumerate(rows):
            mask[r, list(comb)] = True
        _SUBSETS[key] = (rows, mask)
    return _SUBSETS[key]


def necessity_at(rec, j, pool, cfg, recipients, n_pair=30, n_cand=8, kmax=4, share_arc=0.1):
    """Snapshot-sensitivity search for each recipient (local indices). Returns list of dicts.

    Revised after Astra's review (relay 248998526021189632):
    - the empty set and EVERY possible single supplier (every lineage with any release, plus EUK)
      are tested exhaustively, so d = 0, d = 1 and the claim d >= 2 are exact under the
      snapshot model;
    - if no single supplier rescues, pairs are searched among the 30 suppliers with the largest
      share of the metabolites the recipient cannot cover itself, then triples/quadruples among
      the top 8. A pair found gives d = 2 exactly. A larger set found gives only an upper bound
      (d_hi), and no rescue up to size 4 leaves d_hi censored (reported as -1, 'unresolved').
    - 'Snapshot': pools are the recorded stock rescaled by instantaneous release shares, with
      competitors' uptake held fixed. It is a sensitivity, validated separately against
      re-simulated removal trajectories (ko_validation in run_all.py).
    - kappa per support: for each rescuing set S found, X's community S u {X} gets an arc a -> b
      when a supplies >= share_arc of a metabolite b cannot cover alone; kappa_S = 1 if X lies on
      a simple directed cycle of length >= 3 within S u {X} (plexitrophy-computation.md, kappa).
    """
    st, g, inflow, cons, euk = fluxes_at(rec, j, pool, cfg)
    T = float(rec['T'][j])
    P = st['P'][:M]
    tot_in = inflow.sum(0) + euk                                          # (M,)
    thr = cfg['D'] + cfg['mort']
    share = inflow / np.maximum(tot_in, 1e-15)                            # (N, M)
    share_e = euk / np.maximum(tot_in, 1e-15)
    own = (1 - g['lam'])[:, None] * g['y']                                # own usable supply
    needs = own < thr                                                     # (N, M) cannot self-cover m
    out = []
    for i in recipients:
        if g['mu'][i] < thr:
            out.append(dict(i=i, viable=False)); continue
        self_share = share[i]
        donors = np.nonzero(share.sum(1) > 0)[0]
        donors = donors[donors != i]
        ids_all = list(donors) + [EUK]
        sh_all = np.vstack([share[donors], share_e[None, :]])               # (K, M)
        res = dict(i=i, viable=True)
        need_w = needs[i].astype(float) + 1e-3
        rank = np.argsort(-(sh_all * need_w).sum(1))
        for h in ('h0', 'h1'):
            mu0 = _mu_one(i, st, g, pool, cfg, (P * np.clip(self_share, 0, 1))[None, :], T, h)[0]
            if mu0 >= thr:
                res.update({f'd_lo_{h}': 0, f'd_hi_{h}': 0, f'minsets_{h}': [()],
                            f'minsets_complete_{h}': True}); continue
            Pn = P[None, :] * np.clip(sh_all + self_share[None, :], 0, 1)
            ok1 = _mu_one(i, st, g, pool, cfg, Pn, T, h) >= thr
            if ok1.any():
                res.update({f'd_lo_{h}': 1, f'd_hi_{h}': 1,
                            f'minsets_{h}': [(ids_all[k],) for k in np.nonzero(ok1)[0]],
                            f'minsets_complete_{h}': True}); continue
            # exhaustive singletons failed: d >= 2 exactly. Search larger sets among top suppliers.
            found = None
            top = rank[:n_pair]
            pairs = list(itertools.combinations(range(len(top)), 2))
            if pairs:
                A = np.array([sh_all[top[a]] + sh_all[top[b]] for a, b in pairs])
                okp = _mu_one(i, st, g, pool, cfg, P[None, :] * np.clip(A + self_share, 0, 1), T, h) >= thr
                if okp.any():
                    found = (2, [(ids_all[top[pairs[k][0]]], ids_all[top[pairs[k][1]]]) for k in np.nonzero(okp)[0]])
            if found is None:
                c8 = rank[:n_cand]
                rows_, mask = _subsets(len(c8), kmax)
                sizes = mask.sum(1)
                sel = sizes >= 3
                if sel.any():
                    A = mask[sel].astype(float) @ sh_all[c8]
                    oks = _mu_one(i, st, g, pool, cfg, P[None, :] * np.clip(A + self_share, 0, 1), T, h) >= thr
                    if oks.any():
                        ss = sizes[sel]; dmin = int(ss[oks].min())
                        rr = np.nonzero(sel)[0][oks & (ss == dmin)]
                        found = (dmin, [tuple(ids_all[c8[c]] for c in rows_[r]) for r in rr])
            if found is None:
                res.update({f'd_lo_{h}': 2, f'd_hi_{h}': -1, f'minsets_{h}': [],
                            f'minsets_complete_{h}': False})
            else:
                res.update({f'd_lo_{h}': 2, f'd_hi_{h}': found[0], f'minsets_{h}': found[1],
                            f'minsets_complete_{h}': False})
        # functional fan-in vs donor-count dilution (h1), for d >= 2 recipients
        res['n_deficient'] = 0; res['deficient'] = ''; res['fan_class'] = 'n/a'
        if res.get('d_lo_h1', 0) >= 2:
            deficient = []
            for m in range(M):
                Pm = np.tile(P, (1, 1)); Pm[0, m] = P[m] * min(self_share[m], 1.0)
                if _mu_one(i, st, g, pool, cfg, Pm, T, 'h1')[0] < thr:
                    deficient.append(m)
            single_ok, witnessed = {}, {}
            for m in deficient:
                Pm = np.tile(P, (len(ids_all), 1))
                Pm[:, m] = P[m] * np.clip(sh_all[:, m] + self_share[m], 0, 1)
                single_ok[m] = bool((_mu_one(i, st, g, pool, cfg, Pm, T, 'h1') >= thr).any())
                if not single_ok[m]:      # witnessed multi-producer rescue of this one metabolite?
                    c8 = np.argsort(-sh_all[:, m])[:n_cand]
                    rows_, mask = _subsets(len(c8), kmax)
                    sel = mask.sum(1) >= 2
                    A = mask[sel].astype(float) @ sh_all[c8, m]
                    Pm = np.tile(P, (sel.sum(), 1))
                    Pm[:, m] = P[m] * np.clip(A + self_share[m], 0, 1)
                    witnessed[m] = bool((_mu_one(i, st, g, pool, cfg, Pm, T, 'h1') >= thr).any())
            res['n_deficient'] = len(deficient)
            res['deficient'] = ','.join(str(m) for m in deficient)
            if any(not v for v in single_ok.values()):
                res['fan_class'] = ('dilution_witnessed' if all(witnessed.get(m, True) for m in single_ok if not single_ok[m])
                                    else 'dilution_unresolved')
            elif len(deficient) >= 2:
                res['fan_class'] = 'functional_fan_in'
            else:
                res['fan_class'] = 'other'
        # kappa per rescuing support (h1)
        kap = []
        for S in res.get('minsets_h1', []):
            if len(S) < 2:
                kap.append(0); continue
            nodes = [i] + [x for x in S]
            def sup(a, b):      # does a supply a metabolite b cannot cover alone?
                sa = share_e if a == EUK else share[a]
                nb = needs[i] if b == i else (np.zeros(M, bool) if b == EUK else needs[b])
                return bool(np.any((sa >= share_arc) & nb))
            arcs = [(a, b) for a in nodes for b in nodes if a != b and sup(a, b)]
            kap.append(int(i in cycles_through(arcs, nodes, min_len=3, max_len=6)))
        res['kappa_supports'] = kap
        # realised donors (flux-weighted, over ALL suppliers), contrasted with required partners
        real = cons[i][None, :] * sh_all
        tot_need = max(g['mu'][i], 1e-9) * st['B'][i]
        res['realised'] = {ids_all[c]: float(real[c].sum()) for c in range(len(ids_all))
                           if real[c].sum() > 0.02 * tot_need}
        out.append(res)
    return out, g


def cycles_through(edges, nodes, min_len=3, max_len=6):
    """Set of nodes lying on a simple directed cycle with length in [min_len, max_len]."""
    adj = {}
    for a, b in edges:
        adj.setdefault(a, set()).add(b)
    on = set()
    allnodes = sorted(set(adj) | set(x for v in adj.values() for x in v))
    for s in allnodes:             # each cycle is enumerated from its smallest node
        stack = [(s, [s])]
        while stack:
            v, path = stack.pop()
            for w in adj.get(v, ()):
                if w == s:
                    if len(path) >= min_len:
                        on.update(path)
                elif w > s and w not in path and len(path) < max_len:
                    stack.append((w, path + [w]))
    return on


def truth_series(rec, pool, cfg, every_days=14, core_n=200):
    """Biweekly ground truth for the currently most abundant `core_n` lineages.

    plexitrophic_h1 is evaluated per support (P-all over the rescuing sets found): X is
    plexitrophic iff no single supplier rescues it (d_lo >= 2, exact) or every support has
    kappa = 1. A singleton support {Y} can never carry a cycle of length >= 3, so whenever d = 1
    the cycle criterion cannot apply; cycles in a union of alternative supports are NOT counted.
    """
    t = rec['t']
    rows = np.nonzero(np.diff(np.floor(t / every_days), prepend=-1) > 0)[0]
    per_lineage = []; edges = []
    for j in rows:
        B = rec['B'][j]
        recips = np.argsort(-B)[:core_n]
        res, g = necessity_at(rec, j, pool, cfg, recips)
        for r in res:
            i = r['i']
            row = dict(t=float(t[j]), lineage=int(rec['idx'][i]), local=int(i), B=float(B[i]),
                       viable=r['viable'])
            if r['viable']:
                d0, d1 = r['d_lo_h0'], r['d_lo_h1']
                layer = 'genome' if d1 >= 1 else ('expression' if d0 >= 1 else 'none')
                kap = r['kappa_supports']
                kappa_all = int(len(kap) > 0 and all(k == 1 for k in kap) and r['minsets_complete_h1'])
                req = set(x for ms in r['minsets_h1'] for x in ms)
                nec = set.intersection(*[set(ms) for ms in r['minsets_h1']]) if r['minsets_h1'] else set()
                real = set(r['realised'])
                row.update(d_h0=d0, d_h1=d1, d_hi_h1=r['d_hi_h1'],
                           d_unresolved_h1=int(r['d_hi_h1'] == -1),
                           layer=layer, kappa_any_support=int(any(k == 1 for k in kap)),
                           kappa_all_supports=kappa_all,
                           plexitrophic_h1=int(d1 >= 2 or kappa_all),
                           plexitrophic_h0=int(d0 >= 2),
                           n_required=len(req), n_necessary=len(nec), n_realised=len(real),
                           n_realised_not_required=len(real - req),
                           supports_complete=int(r['minsets_complete_h1']),
                           supports_unknown=int(d1 >= 1 and len(r['minsets_h1']) == 0),
                           n_deficient=r['n_deficient'], deficient=r['deficient'], fan_class=r['fan_class'],
                           euk_required=int(EUK in req),
                           genes_lost=int((rec['F'][j][i] < 0.5).sum()))
                for dnr in real:
                    edges.append((float(t[j]), int(dnr if dnr == EUK else rec['idx'][dnr]),
                                  int(rec['idx'][i]), 'realised', r['realised'][dnr]))
                for dnr in req:
                    edges.append((float(t[j]), int(dnr if dnr == EUK else rec['idx'][dnr]),
                                  int(rec['idx'][i]), 'required', float(dnr in nec)))
            per_lineage.append(row)
    return per_lineage, edges
