"""Inference on synthetic observations, scored against ground truth.

Sections: (A) timescales / spectra / AR(1) null; (B) perturbation propagation; (C) aliasing and
closure artefacts; (D) P3-style masked-member recurrence; (E) P1-style song-vs-singer;
(F) mechinterp probes (PCA/NMF features vs true modules; linear probe circuits vs true edges).
"""
import numpy as np
from scipy.optimize import minimize
from scipy.signal import lombscargle

# ============================================================== (A) timescales and spectra
def ou_noise_nll(params, t, y):
    """Exact negative log-likelihood of OU (tau, q) + white noise (r) at irregular times (Kalman)."""
    tau, q, r = np.exp(params)
    m, P = 0.0, q
    nll = 0.0
    prev = t[0]
    for k in range(len(t)):
        if k:
            phi = np.exp(-(t[k] - prev) / tau)
            m = phi * m; P = phi * phi * P + q * (1 - phi * phi)
        S = P + r
        e = y[k] - m
        nll += 0.5 * (np.log(2 * np.pi * S) + e * e / S)
        K = P / S
        m = m + K * e; P = (1 - K) * P
        prev = t[k]
    return nll


def fit_ou_noise(t, y):
    y = y - y.mean()
    v = max(y.var(), 1e-12)
    dt = np.median(np.diff(t))
    best = None
    for tau0 in (dt, 5 * dt, 30 * dt):
        x0 = np.log([tau0, 0.7 * v, 0.3 * v])
        res = minimize(ou_noise_nll, x0, args=(t, y), method='Nelder-Mead',
                       options=dict(maxiter=400, xatol=1e-3, fatol=1e-4))
        if best is None or res.fun < best.fun:
            best = res
    tau, q, r = np.exp(best.x)
    return dict(tau=tau, q=q, r=r, nll=best.fun)


def simulate_ou_noise(t, tau, q, r, rng, nrep):
    n = len(t)
    x = np.empty((nrep, n))
    x[:, 0] = rng.normal(0, np.sqrt(q), nrep)
    for k in range(1, n):
        phi = np.exp(-(t[k] - t[k - 1]) / tau)
        x[:, k] = phi * x[:, k - 1] + rng.normal(0, np.sqrt(q * (1 - phi * phi)), nrep)
    return x + rng.normal(0, np.sqrt(r), (nrep, n))


def deseason(t, y):
    """Remove annual + semiannual harmonics when the series spans >= 300 d (OLS)."""
    if t[-1] - t[0] < 300:
        return y - y.mean()
    w = 2 * np.pi * t / 365.0
    X = np.column_stack([np.ones_like(t), np.sin(w), np.cos(w), np.sin(2 * w), np.cos(2 * w)])
    beta, *_ = np.linalg.lstsq(X, y, rcond=None)
    return y - X @ beta


def freq_grid(t, n=60):
    span = t[-1] - t[0]
    dtm = np.median(np.diff(t))
    fmin, fmax = 1.0 / span, 0.5 / dtm
    return np.exp(np.linspace(np.log(fmin), np.log(fmax), n))


def ls_power(t, Y, freqs):
    """Lomb-Scargle power for rows of Y (normalised to sum 1 per row)."""
    Y = np.atleast_2d(Y)
    out = np.empty((Y.shape[0], len(freqs)))
    for i, y in enumerate(Y):
        yy = y - y.mean()
        out[i] = lombscargle(t, yy, 2 * np.pi * freqs, precenter=False) if yy.std() > 0 else 0
    s = out.sum(1, keepdims=True)
    return out / np.where(s > 0, s, 1)


def band_fractions(t, Y, freqs=None):
    """Fractions of LS power in period bands <30 d, 30-180 d, >180 d (resolvable range only)."""
    if freqs is None:
        freqs = freq_grid(t)
    P = ls_power(t, Y, freqs)
    per = 1 / freqs
    bands = dict(fast=per < 30, mid=(per >= 30) & (per <= 180), slow=per > 180)
    return {b: P[:, msk].sum(1) for b, msk in bands.items()}


def lowhigh_stat(t, Y, freqs):
    P = ls_power(t, Y, freqs)
    k = len(freqs) // 3
    return np.log(P[:, :k].mean(1) + 1e-12) - np.log(P[:, -k:].mean(1) + 1e-12)


def ar1_null_test(t, y, rng, nsur=49):
    """Fit OU+noise to the series; p-value of its low/high spectral ratio vs surrogates."""
    y = (y - y.mean()) / (y.std() + 1e-12)
    fit = fit_ou_noise(t, y)
    freqs = freq_grid(t, 40)
    obs = lowhigh_stat(t, y[None, :], freqs)[0]
    sur = simulate_ou_noise(t, fit['tau'], fit['q'], fit['r'], rng, nsur)
    st = lowhigh_stat(t, sur, freqs)
    p = (1 + np.sum(np.abs(st - np.median(st)) >= abs(obs - np.median(st)))) / (nsur + 1)
    return fit, p


def phi_resolvable(tau, delta, span):
    phi = np.exp(-delta / tau)
    return phi, bool((0.43 * delta <= tau) and (tau <= min(9.5 * delta, span / 5)))


# ============================================================== (B) perturbation propagation
def divergence(Xt, Xc, pre_mask):
    """Per-time divergence between arms, standardised per feature by the control's sd over the
    window (floored at the median feature sd, so near-constant features cannot dominate).
    Features with >20% missing values are dropped; remaining NaNs are ignored."""
    ok = (np.isfinite(Xt).mean(0) > 0.8) & (np.isfinite(Xc).mean(0) > 0.8)
    Xt, Xc = Xt[:, ok], Xc[:, ok]
    sd = np.nanstd(Xc, 0)
    sd = np.maximum(sd, np.nanmedian(sd) + 1e-12)
    return np.nanmean(np.abs(Xt - Xc) / sd, 1)


def response_times(t, div, t0, pre_mask, noise_floor=True):
    """t10 / t50 / tpeak after t0: first times the divergence exceeds 10% and 50% of its
    post-t0 peak above the pre-perturbation baseline. With noise_floor, the 10% threshold is
    raised to baseline mean + 3 sd (observed data)."""
    base = div[pre_mask].mean() if pre_mask.any() else 0.0
    post = t >= t0
    if not post.any():
        return np.nan, np.nan, np.nan
    tp, dp = t[post], div[post]
    peak = dp.max() - base
    if peak <= 0:
        return np.nan, np.nan, np.nan
    thr10 = base + 0.1 * peak
    if noise_floor and pre_mask.sum() > 2:
        thr10 = max(thr10, base + 3 * div[pre_mask].std())
    def first(th):
        j = np.nonzero(dp >= th)[0]
        return tp[j[0]] - t0 if len(j) else np.nan
    return first(thr10), first(base + 0.5 * peak), tp[np.argmax(dp)] - t0


# ============================================================== shared regression helpers
def logit_fit(X, y, lam=1.0, iters=50):
    Xb = np.column_stack([np.ones(len(X)), X])
    w = np.zeros(Xb.shape[1])
    pen = np.full(Xb.shape[1], lam); pen[0] = 1e-6
    for _ in range(iters):
        p = 1 / (1 + np.exp(-np.clip(Xb @ w, -30, 30)))
        g = Xb.T @ (p - y) + pen * w
        H = (Xb * (p * (1 - p))[:, None]).T @ Xb + np.diag(pen)
        step = np.linalg.solve(H, g)
        w -= step
        if np.abs(step).max() < 1e-7:
            break
    return w


def logit_lpd(w, X, y):
    p = 1 / (1 + np.exp(-np.clip(np.column_stack([np.ones(len(X)), X]) @ w, -30, 30)))
    p = np.clip(p, 1e-4, 1 - 1e-4)
    return y * np.log(p) + (1 - y) * np.log(1 - p)


def standardize(Xtr, Xte):
    mu = Xtr.mean(0); sd = Xtr.std(0) + 1e-9
    return (Xtr - mu) / sd, (Xte - mu) / sd


def block_boot_ci(diff, blocks, rng, nboot=2000):
    ub = np.unique(blocks)
    means = np.array([diff[blocks == b].mean() for b in ub])
    bs = np.array([rng.choice(means, len(means)).mean() for _ in range(nboot)])
    return means.mean(), np.percentile(bs, 2.5), np.percentile(bs, 97.5)


def decide(lo, hi, delta):
    if lo > delta:
        return 'support'
    if hi < delta:
        return 'against'
    return 'inconclusive'


# ============================================================== (D) P3 masked-member recurrence
def env_matrix(tday, envT, envS, chl):
    w = 2 * np.pi * tday / 365.0
    return np.column_stack([np.sin(w), np.cos(w), envT, np.log(envS + 1e-6), np.log(chl + 1e-3)])


def core_ranking(counts_tr, ntaxa):
    """Freeze the ranking statistic on training events, weighting events equally."""
    prop = counts_tr[:, :ntaxa] / counts_tr.sum(1, keepdims=True)
    return dict(read_share=prop.mean(0), detect_freq=(counts_tr[:, :ntaxa] > 0).mean(0))


def seasonal_occupancy(counts, tday, ntaxa):
    mon = (np.asarray(tday) % 365 // 30.4).astype(int)
    ms = np.unique(mon)
    occ = np.zeros(ntaxa)
    for m in ms:
        occ += (counts[mon == m, :ntaxa] > 0).any(0)
    return occ / len(ms)


def p3_test(counts, tday, env, train, test, universe, rng, delta=0.05, max_groups=5,
            min_group=3, max_group=6, rcut=0.3):
    """P3-lite. counts: events x taxa (amplicon, last column chloroplast). Returns list of dicts.

    Discovery uses training events only: eligibility (20-80% detection, >=5 detections and
    non-detections), residual log-proportion correlations after env + background gradient,
    greedy cliques (all pairwise r > rcut), outcome = member with prevalence nearest 50%.
    Models: Me (env + background PC1 + log depth), Mm (+ members), Mn (+ matched nonmembers).
    The target is masked before predictor-side normalisation (its reads are removed from the
    depth used to normalise predictors).
    """
    det = counts > 0
    tr = np.asarray(train); te = np.asarray(test)
    prev = det[tr][:, universe].mean(0)
    nd = det[tr][:, universe].sum(0)
    elig = universe[(prev >= 0.2) & (prev <= 0.8) & (nd >= 5) & (len(tr) - nd >= 5)]
    if len(elig) < min_group + 2:
        return []
    depth = counts.sum(1).astype(float)
    logp = np.log((counts + 0.5) / depth[:, None])
    # background gradient: PC1 of universe taxa not eligible (fit on train only)
    bg = np.setdiff1d(universe, elig)
    Xbg = logp[:, bg]
    mu = Xbg[tr].mean(0)
    _, _, Vt = np.linalg.svd(Xbg[tr] - mu, full_matrices=False)
    pc1 = (Xbg - mu) @ Vt[0]
    E = np.column_stack([env, pc1, np.log(depth)])
    # residual correlations (train)
    Xe = np.column_stack([np.ones(len(tr)), E[tr]])
    beta, *_ = np.linalg.lstsq(Xe, logp[tr][:, elig], rcond=None)
    res = logp[tr][:, elig] - Xe @ beta
    C = np.corrcoef(res.T)
    np.fill_diagonal(C, 0)
    used = set(); groups = []
    order = np.dstack(np.unravel_index(np.argsort(-C, axis=None), C.shape))[0]
    for a, b in order:
        if C[a, b] < rcut or len(groups) >= max_groups:
            break
        if a in used or b in used or a >= b:
            continue
        grp = [a, b]
        cand = [c for c in np.argsort(-(C[a] + C[b])) if c not in used and c not in grp]
        for c in cand:
            if len(grp) >= max_group:
                break
            if all(C[c, g] > rcut for g in grp):
                grp.append(c)
        if len(grp) >= min_group:
            groups.append(grp); used.update(grp)
    peak_month = np.array([np.argmax(np.bincount((np.asarray(tday)[tr] % 365 // 30.4).astype(int),
                                                 weights=det[tr][:, x], minlength=13)) for x in elig])
    meanab = logp[tr][:, elig].mean(0)
    out = []
    for grp in groups:
        taxa = elig[grp]
        pv = det[tr][:, taxa].mean(0)
        o = int(np.argmin(np.abs(pv - 0.5)))
        target = taxa[o]; preds = np.delete(taxa, o)
        # matched nonmembers: nearest on (logit prev, circular peak month, mean log abundance)
        nonm = [k for k in range(len(elig)) if k not in grp]
        feats = np.column_stack([np.log(np.clip(det[tr][:, elig].mean(0), .01, .99) /
                                        (1 - np.clip(det[tr][:, elig].mean(0), .01, .99))),
                                 np.cos(2 * np.pi * peak_month / 12), np.sin(2 * np.pi * peak_month / 12),
                                 meanab])
        sd = feats.std(0) + 1e-9
        matched = []
        for pidx in [grp[k] for k in range(len(grp)) if k != o]:
            dd = (((feats[nonm] - feats[pidx]) / sd) ** 2).sum(1)
            for k in np.argsort(dd):
                if nonm[k] not in matched:
                    matched.append(nonm[k]); break
        matched = elig[matched]
        y = det[:, target].astype(float)
        depth_m = depth - counts[:, target]                   # mask target before normalising
        def lp(cols):
            return np.log((counts[:, cols] + 0.5) / depth_m[:, None])
        Xe_ = E.copy(); Xe_[:, -1] = np.log(depth)            # target detection depends on depth
        models = dict(Me=Xe_, Mm=np.column_stack([Xe_, lp(preds)]), Mn=np.column_stack([Xe_, lp(matched)]))
        lpd = {}
        for k, X in models.items():
            Xtr, Xte = standardize(X[tr], X[te])
            w = logit_fit(Xtr, y[tr], lam=2.0)
            lpd[k] = logit_lpd(w, Xte, y[te])
        blocks = (np.asarray(tday)[te] // 28).astype(int)
        g_me = block_boot_ci(lpd['Mm'] - lpd['Me'], blocks, rng)
        g_mn = block_boot_ci(lpd['Mm'] - lpd['Mn'], blocks, rng)
        if g_me[1] > delta and g_mn[1] > 0:
            dec = 'support'
        elif g_me[2] < delta:
            dec = 'against'
        else:
            dec = 'inconclusive'
        out.append(dict(target=int(target), preds=[int(x) for x in preds],
                        matched=[int(x) for x in matched], gain_me=g_me, gain_mn=g_mn,
                        decision=dec, n_train=len(tr), n_test=len(te),
                        test_prev=float(y[te].mean())))
    return out


# ============================================================== (E) P1 song vs singer
def gauss_fit_lpd(Xtr, ytr, Xte, yte, lam=1.0):
    Xtr, Xte = standardize(Xtr, Xte)
    Xb = np.column_stack([np.ones(len(Xtr)), Xtr]); Xt = np.column_stack([np.ones(len(Xte)), Xte])
    pen = np.eye(Xb.shape[1]) * lam; pen[0, 0] = 0
    w = np.linalg.solve(Xb.T @ Xb + pen, Xb.T @ ytr)
    s2 = np.mean((ytr - Xb @ w) ** 2) * len(ytr) / max(len(ytr) - Xb.shape[1], 1)
    e = yte - Xt @ w
    return -0.5 * (np.log(2 * np.pi * s2) + e * e / s2), w


# ============================================================== (F) mechinterp probes
def nmf(X, k, rng, iters=400):
    W = rng.random((X.shape[0], k)) + 0.1; H = rng.random((k, X.shape[1])) + 0.1
    for _ in range(iters):
        H *= (W.T @ X) / (W.T @ W @ H + 1e-12)
        W *= (X @ H.T) / (W @ H @ H.T + 1e-12)
    return W, H


def module_alignment(loadings, modules, rng, nnull=500):
    """loadings: k x taxa. modules: {name: bool mask over taxa}. Best |cosine| vs random sets."""
    Lz = loadings / (np.linalg.norm(loadings, axis=1, keepdims=True) + 1e-12)
    out = {}
    for name, msk in modules.items():
        if msk.sum() < 2 or msk.sum() > len(msk) - 2:
            continue
        def score(m):
            v = m.astype(float) - m.mean(); v /= np.linalg.norm(v) + 1e-12
            return np.abs(Lz @ v).max()
        s = score(msk)
        null = np.array([score(rng.permutation(msk)) for _ in range(nnull)])
        out[name] = dict(score=float(s), null_mean=float(null.mean()),
                         pct=float((null < s).mean()), size=int(msk.sum()))
    return out


def auroc(scores, labels):
    labels = np.asarray(labels, bool); scores = np.asarray(scores, float)
    if labels.sum() == 0 or (~labels).sum() == 0:
        return np.nan
    order = np.argsort(scores)
    ranks = np.empty(len(scores)); ranks[order] = np.arange(1, len(scores) + 1)
    # average ties
    _, inv, cnt = np.unique(scores, return_inverse=True, return_counts=True)
    sums = np.bincount(inv, ranks); ranks = (sums / cnt)[inv]
    n1 = labels.sum(); n0 = (~labels).sum()
    return (ranks[labels].sum() - n1 * (n1 + 1) / 2) / (n1 * n0)


def ridge(X, Y, lam):
    X = (X - X.mean(0)) / (X.std(0) + 1e-9)
    Yc = Y - Y.mean(0)
    return np.linalg.solve(X.T @ X + lam * np.eye(X.shape[1]), X.T @ Yc)


# ============================================================== (D') P3 aligned with the protocol
def logit_fit_pen(X, y, pen, iters=60):
    """Ridge-logistic with per-coefficient penalties (pen excludes the intercept)."""
    Xb = np.column_stack([np.ones(len(X)), X])
    P_ = np.concatenate([[1e-6], pen])
    w = np.zeros(Xb.shape[1])
    for _ in range(iters):
        p = 1 / (1 + np.exp(-np.clip(Xb @ w, -30, 30)))
        g = Xb.T @ (p - y) + P_ * w
        H = (Xb * (p * (1 - p))[:, None]).T @ Xb + np.diag(P_)
        step = np.linalg.solve(H, g)
        w -= step
        if np.abs(step).max() < 1e-7:
            break
    return w


def pc1_power(Xtr, Xall, iters=25, seed=0):
    mu = Xtr.mean(0)
    A = Xtr - mu
    v = np.random.default_rng(seed).normal(size=A.shape[1])
    for _ in range(iters):
        v = A.T @ (A @ v); v /= np.linalg.norm(v) + 1e-12
    return (Xall - mu) @ v


def p3_discover(Y, featfn_disc, E, tr, universe, min_group=3, max_group=6, rcut=0.3, max_groups=5):
    det = Y
    prev = det[tr][:, universe].mean(0)
    nd = det[tr][:, universe].sum(0)
    elig = universe[(prev >= 0.2) & (prev <= 0.8) & (nd >= 5) & (len(tr) - nd >= 5)]
    if len(elig) < min_group + 2:
        return elig, []
    L = featfn_disc(elig)
    Xe = np.column_stack([np.ones(len(tr)), E[tr]])
    beta, *_ = np.linalg.lstsq(Xe, L[tr], rcond=None)
    res = L[tr] - Xe @ beta
    C = np.corrcoef(res.T); np.fill_diagonal(C, 0); C = np.nan_to_num(C)
    used = set(); groups = []
    order = np.dstack(np.unravel_index(np.argsort(-C, axis=None), C.shape))[0]
    for a, b in order:
        if C[a, b] < rcut or len(groups) >= max_groups:
            break
        if a in used or b in used or a >= b:
            continue
        grp = [a, b]
        for c in np.argsort(-(C[a] + C[b])):
            if len(grp) >= max_group:
                break
            if c not in used and c not in grp and all(C[c, g_] > rcut for g_ in grp):
                grp.append(c)
        if len(grp) >= min_group:
            groups.append([int(elig[g_]) for g_ in grp]); used.update(grp)
    return elig, groups


def p3_score_group(Y, featfn, Ebase, tday, tr, te, elig, group, universe, rng, n_ctrl=20,
                   lam_taxa=5.0, nboot=1000):
    """Rotated targets; 20 matched control sets; background V minus (C u N_k), target-masked;
    shrinkage only on added taxon slopes. Returns per-event gains averaged over targets."""
    det = Y
    peak = np.array([np.argmax(np.bincount((np.asarray(tday)[tr] % 365 // 30.4).astype(int),
                                           weights=det[tr][:, x], minlength=13)) for x in elig])
    pv = np.clip(det[tr][:, elig].mean(0), .01, .99)
    feats = np.column_stack([np.log(pv / (1 - pv)), np.cos(2 * np.pi * peak / 12),
                             np.sin(2 * np.pi * peak / 12), featfn(elig, None)[tr].mean(0)])
    sd = feats.std(0) + 1e-9
    pos = {x: k for k, x in enumerate(elig)}
    nonm = [x for x in elig if x not in group]
    g_me = np.zeros(len(te)); g_mn = np.zeros(len(te))
    n_t = 0
    for t_ in group:
        preds = [x for x in group if x != t_]
        if len(nonm) < 2 * len(preds):
            return None          # specificity unresolved
        y = det[:, t_].astype(float)
        gm_k, gn_k = [], []
        for k in range(n_ctrl):
            order = rng.permutation(len(preds)); chosen = []
            for o in order:
                pf = feats[pos[preds[o]]]
                cand = [x for x in nonm if x not in chosen]
                dd = [(((feats[pos[x]] - pf) / sd) ** 2).sum() for x in cand]
                chosen.append(cand[int(np.argmin(dd))])
            bg = np.setdiff1d(universe, np.array(group + chosen))
            Fbg = featfn(bg, t_)
            pc = pc1_power(Fbg[tr], Fbg, seed=k)
            E = np.column_stack([Ebase, pc])
            Fm = featfn(np.array(preds), t_); Fn = featfn(np.array(chosen), t_)
            lp = {}
            for nm, X, ntax in (('Me', E, 0), ('Mm', np.column_stack([E, Fm]), Fm.shape[1]),
                                ('Mn', np.column_stack([E, Fn]), Fn.shape[1])):
                Xtr, Xte = standardize(X[tr], X[te])
                pen = np.concatenate([np.full(E.shape[1], 1e-3), np.full(ntax, lam_taxa)])
                w = logit_fit_pen(Xtr, y[tr], pen)
                lp[nm] = logit_lpd(w, Xte, y[te])
            gm_k.append(lp['Mm'] - lp['Me']); gn_k.append(lp['Mm'] - lp['Mn'])
        g_me += np.mean(gm_k, 0); g_mn += np.mean(gn_k, 0); n_t += 1
    g_me /= n_t; g_mn /= n_t
    blocks = (np.asarray(tday)[te] // 28).astype(int)
    return dict(gain_me=block_boot_ci(g_me, blocks, rng, nboot), gain_mn=block_boot_ci(g_mn, blocks, rng, nboot),
                n_targets=n_t, n_test=len(te))


def p3_verdict(r, delta=0.02, min_events=10):
    if r is None:
        return 'unresolved_specificity'
    if r['n_test'] < min_events:
        return 'inconclusive'
    if r['gain_me'][1] > delta and r['gain_mn'][1] > delta:
        return 'supported'
    if r['gain_me'][2] < delta or r['gain_mn'][2] < delta:
        return 'falsified'
    return 'inconclusive'
