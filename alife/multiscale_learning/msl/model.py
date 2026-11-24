"""Multiscale ecosystem model: environment -> expression -> population -> genome -> nucleotide.

Time unit: days. Explicit integration with dt = 0.1 d (2.4 h).
Layers and their characteristic response times (set by parameters, not hard-coded labels):
  expression  : relaxation tau_e = 0.2 d towards regulatory targets set by local pools
  population  : growth/dilution, ~1/(D) = 10 d
  genome      : within-lineage gene-carrier frequency F under selection/loss/HGT, months-years
  nucleotide  : SNV allele frequencies under drift (Ne scaled by abundance) + weak fluctuating
                selection + hitchhiking with gene sweeps, years
  environment : seasonal + AR(1) temperature, seasonal resource inflow, storm pulses, diel light

Nothing about who depends on whom is specified. Lineages start prototrophic (F ~ 1 for all
biosynthesis genes); dependencies arise only through (a) gene loss when leaked metabolites make a
gene's cost exceed its benefit (genome layer), and (b) regulated repression of own synthesis when
the external pool is high (expression layer).
"""
import numpy as np

DT = 0.1
SLOW_EVERY = 5           # slow layers integrated with step 0.5 d
R = 3                    # resources (0 = photosynthate with diel inflow)
M = 5                    # essential leaky metabolites (biosynthesis genes)
L = 8                    # SNV loci per lineage (0: linked to gene 0; 1-2 selected; 3-7 neutral)

# ---------------------------------------------------------------- species pool (shared traits)
def make_pool(pool_seed=7, n_pool=1500, tail='powerlaw', zipf=1.6, ln_sd=2.3):
    """Regional pool of lineage types with fixed traits, shared by all stations of one study."""
    rng = np.random.default_rng(pool_seed)
    guild = np.arange(n_pool) % R
    a = np.full((n_pool, R), 0.12)
    a[np.arange(n_pool), guild] = 1.0
    # regional (metacommunity) abundance: power law over a random rank order. This sets the
    # immigration supply, so locally unfit lineages persist as a rare tail (mass effect) and can
    # bloom from rarity when local conditions favour them.
    # tail='powerlaw': Zipf weights rank^-zipf; tail='lognormal': exp(N(0, ln_sd)). Both use the
    # same trait draws (traits are drawn before the tail), so only the abundance tail differs.
    rng_tail = np.random.default_rng([pool_seed, 99])
    if tail == 'powerlaw':
        rank = rng_tail.permutation(n_pool) + 1
        w_reg = rank.astype(float) ** (-zipf)
    else:
        w_reg = np.exp(rng_tail.normal(0, ln_sd, n_pool))
    w_reg /= w_reg.sum()
    return dict(
        n=n_pool, guild=guild, a=a, w_reg=w_reg, tail=tail,
        Topt=rng.uniform(3, 23, n_pool),          # temperature optimum (seasonal niche)
        lam=rng.uniform(0.05, 0.45, n_pool),      # leakiness of own synthesis
        regulated=rng.random(n_pool) < 0.5,       # represses own synthesis when pool is high
        Krep=rng.uniform(0.15, 0.6, n_pool),
        form_b=rng.random(n_pool) < 0.4,          # makes the 'pseudo' vitamer of metabolite 0
        tox_sens=rng.random(n_pool) < 0.5,        # sensitive to the perturbation toxicant
        hc_deg=rng.random(n_pool) < 0.25,         # can use the perturbation substrate
        vfast=np.exp(rng.normal(0, 0.35, n_pool)), # gleaner-opportunist axis: Vmax*v, KS*v^2.5
        gsize=rng.uniform(2.0, 6.0, n_pool),      # Mb, for read sampling
        snv_sign=rng.choice([-1.0, 1.0], (n_pool, L)),
        snv_Tref=rng.uniform(5, 20, (n_pool, L)),
    )


DEFAULTS = dict(
    n_lin=1200, burn_years=4, obs_years=6, D=0.1, mort=0.01, Vmax=1.6, KS=0.4, Y=0.5, maint=0.03,
    Tw=7.0, sig=1.6, csyn=0.008, cgene=0.007, vup=2.0, KP=0.3, Krs_mult=0.15, basal=0.05,
    tau_e=0.2, uloss=2e-5, uhgt=2e-5, Ne_g=1e5, Fmin=1e-8,
    Ne_snv=4000.0, kappa_snv=1/1500., s_snv=0.012, link=0.9,
    leak=True, regulation=True, vitamer=False, imm=0.02,
    cost_mult=(3.0, 1.0, 1.0, 1.0, 1.0),   # metabolite 0 (B12-like) is costly to make
    Sin=(3.0, 2.0, 2.0), Samp=(0.6, 0.5, 0.3), Sphase=(0.0, 1.2, 2.6), diel=0.9,
    T_mean=12.0, T_amp=10.0, T_ar_sd=1.5, T_ar_tau=5.0, storm_rate=1/30., storm_mult=3.0,
    station_dT=0.0, perturb=None,
    # phototrophic eukaryotes (algae A, grazers Z): 'dynamic' bloom-crash cycle, 'flat' (A, Z held
    # constant at A_flat/Z_flat with mean light and temperature: same mean fluxes, no bloom cycle),
    # or 'none' (exogenous diel photosynthate only).
    euk='dynamic', A_flat=1.0, Z_flat=0.5, muA_max=1.1, Acap=8.0, TA_opt=9.0, TA_w=9.0,
    gmax=1.0, Kg=2.0, eZ=0.25, mZ=0.12, immA=1e-3, immZ=2e-4, exA=0.35, sloppy=0.3,
    relA=(0.0, 0.10, 0.10, 0.08, 0.0), qB12=0.04, KB12=0.03, exo0_frac=0.25,
)


def make_env(cfg, rng, n_steps):
    t = np.arange(n_steps) * DT
    # AR(1) temperature anomaly with tau T_ar_tau
    phi = np.exp(-DT / cfg['T_ar_tau'])
    eps = rng.normal(0, cfg['T_ar_sd'] * np.sqrt(1 - phi**2), n_steps)
    anom = np.empty(n_steps); a = 0.0
    for k in range(n_steps):   # small loop, fine
        a = phi * a + eps[k]; anom[k] = a
    T = cfg['T_mean'] + cfg['station_dT'] + cfg['T_amp'] * np.sin(2*np.pi*(t - 110)/365.) + anom
    Sin = np.empty((n_steps, R))
    for r in range(R):
        Sin[:, r] = cfg['Sin'][r] * (1 + cfg['Samp'][r] * np.sin(2*np.pi*t/365. + cfg['Sphase'][r]))
    # diel light on photosynthate inflow, gated by season (more light in summer)
    light = np.clip(np.sin(2*np.pi*(t % 1.0 - 0.25)), 0, None) * np.pi  # mean 1 over a day
    diel = (1 - cfg['diel']) + cfg['diel'] * light
    season_light = 0.55 + 0.45 * np.sin(2*np.pi*(t - 80)/365.)
    if cfg['euk'] == 'none':
        Sin[:, 0] *= diel            # all photosynthate exogenous, diel
    else:
        Sin[:, 0] *= cfg['exo0_frac']  # allochthonous part; algal exudate added in the step
    # storm pulses on resource 1 (unmeasured stochastic forcing), exponential decay tau 3 d
    storms = np.zeros(n_steps)
    ev = np.nonzero(rng.random(n_steps) < cfg['storm_rate'] * DT)[0]
    kern = np.exp(-np.arange(int(15/DT)) * DT / 3.0)
    for e in ev:
        seg = storms[e:e+len(kern)]; seg += (cfg['storm_mult'] - 1) * kern[:len(seg)]
    Sin[:, 1] *= (1 + storms)
    return t, T, Sin, storms, season_light, diel


def _rowmin(X):
    m = X[:, 0].copy()
    for c in range(1, X.shape[1]):
        np.minimum(m, X[:, c], out=m)
    return m


def growth(st, pool, cfg, T, sel=False):
    """Instantaneous rates for all lineages. st: dict of state arrays. Returns dict."""
    B, S, P, Esyn, Eup, Ecat, F = (st[k] for k in ('B', 'S', 'P', 'Esyn', 'Eup', 'Ecat', 'F'))
    idx = st['idx']
    a = pool['a'][idx]; lam = pool['lam'][idx] if cfg['leak'] else np.zeros(len(idx))
    Tfac = np.exp(-((T - pool['Topt'][idx]) / cfg['Tw'])**2)
    v = pool['vfast'][idx]
    KSi = cfg['KS'] * v**2.5
    sat = S[None, :] / (KSi[:, None] + S[None, :])
    Jr = (cfg['Vmax'] * v * Tfac)[:, None] * Ecat * a * sat          # per-biomass uptake
    if st.get('hc') is not None:                                          # perturbation substrate
        Jhc = cfg['Vmax'] * v * Tfac * pool['hc_deg'][idx] * st['hc'] / (KSi + st['hc'])
    else:
        Jhc = 0.0
    J = Jr.sum(1) + Jhc
    y = F * Esyn * cfg['sig'] * Tfac[:, None]                              # synthesis capacity
    Puse = P[:M]
    upt = Eup * cfg['vup'] * Puse / (cfg['KP'] + Puse)
    supply = (1 - lam)[:, None] * y + upt
    cm = np.asarray(cfg['cost_mult'])
    muE = cfg['Y'] * J - cfg['csyn'] * (y * cm).sum(1) - cfg['cgene'] * (F * cm).sum(1) - cfg['maint']
    mu = np.minimum(muE, _rowmin(supply))
    out = dict(mu=mu, muE=muE, y=y, upt=upt, supply=supply, Jr=Jr, J=J, lam=lam, Tfac=Tfac,
               Jhc=Jhc)
    if sel:   # carrier vs non-carrier growth for each gene (selection coefficient)
        s = np.empty_like(F)
        for m in range(M):
            ycar = Esyn[:, m] * cfg['sig'] * Tfac
            sup_c = supply.copy(); sup_c[:, m] = (1 - lam) * ycar + upt[:, m]
            sup_n = supply.copy(); sup_n[:, m] = upt[:, m]
            base = muE + cm[m] * (cfg['csyn'] * y[:, m] + cfg['cgene'] * F[:, m])
            mu_c = np.minimum(base - cm[m] * (cfg['csyn'] * ycar + cfg['cgene']), _rowmin(sup_c))
            mu_n = np.minimum(base, _rowmin(sup_n))
            s[:, m] = mu_c - mu_n
        out['s'] = s
    return out


def pool_fluxes(st, g, pool, cfg):
    """Leak inflow per donor into each pool and per-recipient consumption from pools."""
    idx = st['idx']; B = st['B']
    leak = (g['lam'][:, None] * g['y']) * B[:, None]                     # (N, M)
    inflow = np.zeros((len(B), M + 1))
    inflow[:, :M] = leak
    if cfg['vitamer']:
        fb = pool['form_b'][idx]
        inflow[fb, M] = leak[fb, 0]; inflow[fb, 0] = 0.0                  # pseudo form unusable
    mup = np.clip(g['mu'], 0, None)
    use_frac = np.clip(mup[:, None] / np.maximum(g['supply'], 1e-12), 0, 1)
    cons = g['upt'] * use_frac * B[:, None]                               # (N, M) from pools
    return inflow, cons


def init_state(pool, cfg, rng, idx):
    n = len(idx)
    return dict(
        idx=idx, B=2.0 * pool['w_reg'][idx] / pool['w_reg'][idx].sum() * rng.uniform(0.5, 1.5, n), S=np.array(cfg['Sin'], float),
        P=np.full(M + 1, 0.05), Esyn=np.ones((n, M)), Eup=np.full((n, M), 0.5),
        Ecat=np.full((n, R), 1.0 / R), F=np.full((n, M), 1 - 1e-3),
        Psnv=rng.beta(2, 2, (n, L)), hc=0.0, A=0.5, Z=0.2,
    )


def simulate(pool, cfg_over=None, seed=1, idx=None, record=True, start_state=None,
             t_start_step=0, n_steps=None, env=None, knockout=None, rec_offsets=None,
             snap_offsets=()):
    """Run the model. Returns dict with recorded ground truth (observation window only).

    knockout: donor local index or array of indices -> their leak no longer reaches the pools
    (counterfactual removal; the leaked fraction is lost, the donors' own growth is unchanged).
    """
    cfg = dict(DEFAULTS); cfg.update(cfg_over or {})
    rng_env = np.random.default_rng([seed, 1])
    rng_dyn = np.random.default_rng([seed, 2])
    if idx is None:
        idx = np.sort(np.random.default_rng([seed, 3]).choice(pool['n'], cfg['n_lin'], replace=False))
    total_steps = int(round((cfg['burn_years'] + cfg['obs_years']) * 365 / DT))
    if env is None:
        env = make_env(cfg, rng_env, total_steps)
    t_all, T_all, Sin_all, storms, slight_all, diel_all = env
    burn_steps = int(round(cfg['burn_years'] * 365 / DT))
    st = init_state(pool, cfg, np.random.default_rng([seed, 4]), idx) if start_state is None \
        else {k: (v.copy() if isinstance(v, np.ndarray) else v) for k, v in start_state.items()}
    if n_steps is None:
        n_steps = total_steps - t_start_step
    n = len(idx)
    pert = cfg['perturb']
    Krep = pool['Krep'][idx] * cfg['KP']
    Krs = cfg['Krs_mult'] * cfg['KP']
    imm = cfg['imm'] * pool['w_reg'][idx] / pool['w_reg'][idx].sum()   # total immigration flux
    recs = {k: [] for k in ('t', 'T', 'S', 'P', 'B', 'Esyn', 'Eup', 'Ecat', 'F', 'Psnv', 'mu', 'A', 'Z', 'euk_doc',
                            'Sin', 'storm', 'hc')}
    snaps = {}
    snap_steps = set()
    # recording: by default daily at noon (offset 5 steps = 0.5 d) through the obs window; designs
    # add their own (day, hour) sample times so sub-daily states exist exactly where sampled.
    if rec_offsets is None:
        rec_offsets = np.arange(5, total_steps - burn_steps, 10)
    rec_mask = np.zeros(total_steps + 1, bool)
    rec_mask[burn_steps + np.asarray(rec_offsets, int)] = True
    snap_steps = set(int(burn_steps + o) for o in snap_offsets)
    for k in range(t_start_step, t_start_step + n_steps):
        T = T_all[k]; Sin = Sin_all[k]; tday = t_all[k]
        tox = np.zeros(n)
        for pe in (pert or []):        # list of perturbation events (CMO-like dosing)
            if pe['t0'] <= tday < pe['t0'] + pe['dur']:
                tox = pe['tox'] * pool['tox_sens'][idx]
                st['hc'] += pe['hc_in'] * DT / pe['dur']
        slow = (k % SLOW_EVERY == 0)          # genome/SNV updated every 0.5 d (alternating hours)
        g = growth(st, pool, cfg, T, sel=slow)
        inflow, cons = pool_fluxes(st, g, pool, cfg)
        if knockout is not None:       # removed donors' leak is lost (not returned to the donor)
            inflow = inflow.copy(); inflow[knockout] = 0.0
        # resources
        dS = cfg['D'] * (Sin - st['S']) - (g['Jr'] * st['B'][:, None]).sum(0)
        st['S'] = np.maximum(st['S'] + dS * DT, 1e-9)
        if st['hc'] > 0:
            st['hc'] = max(st['hc'] + (-cfg['D'] * st['hc'] - np.sum(g['Jhc'] * st['B'])) * DT, 0.0)
        # phototrophic eukaryotes: bloom-crash with grazers; DOC exudate (diel), metabolite release
        # on grazing/lysis, and B12 (metabolite 0, form A only) demand
        dP = inflow.sum(0) - cfg['D'] * st['P']
        dP[:M] -= cons.sum(0)
        if cfg['euk'] != 'none':
            dyn = cfg['euk'] == 'dynamic'
            A = st['A'] if dyn else cfg['A_flat']; Z = st['Z'] if dyn else cfg['Z_flat']
            sl = slight_all[k] if dyn else 0.55
            TA = T if dyn else cfg['T_mean'] + cfg['station_dT']
            limB = st['P'][0] / (cfg['KB12'] + st['P'][0])
            muA = cfg['muA_max'] * sl * np.exp(-((TA - cfg['TA_opt'])/cfg['TA_w'])**2) * limB \
                * max(1 - A / cfg['Acap'], 0.0)
            graz = cfg['gmax'] * Z * A / (A + cfg['Kg'])
            prod = muA * A
            st['euk_doc'] = cfg['exA'] * prod * diel_all[k] + cfg['sloppy'] * graz
            st['S'][0] += st['euk_doc'] * DT
            dP[:M] += np.array(cfg['relA']) * graz
            dP[0] -= cfg['qB12'] * prod
            if dyn:
                st['A'] = max(A + (prod - graz - cfg['D'] * A + cfg['immA']) * DT, 1e-6)
                st['Z'] = max(Z + (cfg['eZ'] * graz - cfg['mZ'] * Z + cfg['immZ']) * DT, 1e-6)
        st['P'] = np.maximum(st['P'] + dP * DT, 1e-9)
        # population
        r = g['mu'] * (1 - tox) - cfg['D'] - cfg['mort'] - tox * 0.5
        st['B'] = st['B'] * np.exp(r * DT) + imm * DT
        # expression (fast) : relax to regulatory targets set by local pools
        Puse = st['P'][:M]
        if cfg['regulation']:
            reg = pool['regulated'][idx]
            rep = 1.0 / (1.0 + (Puse[None, :] / Krep[:, None])**2)
            tgt_syn = np.where(reg[:, None], cfg['basal'] + (1 - cfg['basal']) * rep, 1.0)
        else:
            tgt_syn = np.ones((n, M))
        own = st['F'] * st['Esyn']
        tgt_up = 0.05 + 0.95 * (1 - 0.8 * own) / (1 + Puse[None, :] / Krs)
        KSi = cfg['KS'] * pool['vfast'][idx]**2.5
        w = pool['a'][idx] * st['S'][None, :] / (KSi[:, None] + st['S'][None, :]) + 0.02
        tgt_cat = w / w.sum(1, keepdims=True)
        fe = DT / cfg['tau_e']
        st['Esyn'] += (tgt_syn - st['Esyn']) * fe
        st['Eup'] += (tgt_up - st['Eup']) * fe
        st['Ecat'] += (tgt_cat - st['Ecat']) * fe
        if slow:
            DTs = DT * SLOW_EVERY
            # genome (slow): carrier frequency under selection, loss, HGT, drift
            F = st['F']
            Bsum = st['B'].sum()
            carriers = (st['B'][:, None] * F).sum(0) / max(Bsum, 1e-12)
            dF = (g['s'] * F * (1 - F) - cfg['uloss'] * F + cfg['uhgt'] * (1 - F) * carriers[None, :]) * DTs
            dF += np.sqrt(F * (1 - F) * DTs / cfg['Ne_g']) * rng_dyn.standard_normal(F.shape)
            Fn = np.clip(F + dF, cfg['Fmin'], 1 - cfg['Fmin'])
            # nucleotide (slowest): drift scaled by abundance, fluctuating selection, hitchhiking
            p = st['Psnv']
            Ne = cfg['Ne_snv'] * np.clip(st['B'] / 0.01, 0.01, 50.0)
            drift = np.sqrt(p * (1 - p) * DTs / Ne[:, None]) * rng_dyn.standard_normal(p.shape)
            sel = np.zeros_like(p)
            sel[:, 1:3] = cfg['s_snv'] * pool['snv_sign'][idx][:, 1:3] * \
                (T - pool['snv_Tref'][idx][:, 1:3]) / 10.0
            dp = cfg['kappa_snv'] * (0.5 - p) * DTs + sel * p * (1 - p) * DTs + drift
            dp[:, 0] += cfg['link'] * (Fn[:, 0] - F[:, 0])
            st['Psnv'] = np.clip(p + dp, 1e-4, 1 - 1e-4)
            st['F'] = Fn
        if record and rec_mask[k]:
            recs['t'].append(tday - burn_steps * DT); recs['T'].append(T)
            recs['S'].append(st['S'].copy()); recs['P'].append(st['P'].copy())
            recs['B'].append(st['B'].copy()); recs['Esyn'].append(st['Esyn'].copy())
            recs['Eup'].append(st['Eup'].copy()); recs['Ecat'].append(st['Ecat'].copy())
            recs['F'].append(st['F'].copy()); recs['Psnv'].append(st['Psnv'].copy())
            recs['mu'].append(g['mu'].copy()); recs['Sin'].append(Sin.copy())
            recs['storm'].append(storms[k]); recs['hc'].append(st['hc'])
            recs['A'].append(st['A']); recs['Z'].append(st['Z']); recs['euk_doc'].append(st.get('euk_doc', 0.0))
        if k in snap_steps:
            snaps[k] = {kk: (v.copy() if isinstance(v, np.ndarray) else v) for kk, v in st.items()}
    out = {k: (np.array(v, dtype=np.float64) if k == 't' else np.array(v, dtype=np.float32)) for k, v in recs.items()}
    out.update(cfg=cfg, idx=idx, seed=seed, final=st, snaps=snaps, env=env,
               burn_steps=burn_steps)
    return out
