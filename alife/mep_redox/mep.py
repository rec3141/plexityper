"""Dissipation yardsticks (revised after Astra 248992826624753664: 'best-found reference', not
'MEP optimum') on the same system: enzyme fields (the 'song') optimised for
entropy production at a fixed catalyst budget, over different horizons, plus a greedy
instantaneous (short-horizon) policy and random allocations."""
import numpy as np
from scipy.optimize import minimize
from model import Model, FixedEnzyme, enzyme_field, NR, NS, S, T

BUDGET_FLOOR = 1e-3


def fixed_inputs(m, y):
    C, B = m.unpack(y)
    P = m.tendencies(0.0, y, parts=True)
    E = enzyme_field(m, C, B)
    necro = m.p["f_necro"] * P["mortB"][m.ialg]      # (v3) algal necromass only
    return C.copy(), E, B[m.ialg].copy(), necro


def family(E0, theta, z):
    """E_r(z) = E0_r(z - shift_r) * exp(a_r) + floor; renormalised to E0's total."""
    a, sh = theta[:NR], theta[NR:]
    tot = E0.sum()
    floor = BUDGET_FLOOR * tot / E0.size
    E = np.empty_like(E0)
    for r in range(NR):
        E[r] = np.interp(z - sh[r], z, E0[r]) * np.exp(a[r]) + floor
    return E * tot / E.sum()


def sigma_fixed(p, E, Balg, necro, C0, t_end, n_eval=60, greedy=None):
    """Mean chemotrophic sigma over [0, t_end] and at t_end (J K^-1 m^-2 d^-1)."""
    f = FixedEnzyme(p, E, Balg, necro, greedy_beta=greedy)
    if t_end >= 50:
        te = np.unique(np.concatenate([[0], np.geomspace(0.05, t_end, n_eval)]))
    else:
        te = np.linspace(0, t_end, n_eval)
    try:
        sol = f.run(np.asarray(C0, float).ravel(), t_end, t_eval=te, rtol=1e-4, atol=1e-6)
    except RuntimeError:
        return -1.0, -1.0, None
    sig = np.array([f.sigma_parts(sol.y[:, k], sol.t[k])[0] for k in range(len(sol.t))])
    mean = np.trapezoid(sig, sol.t) / t_end
    return float(mean), float(sig[-1]), (sol.t, sig, sol.y[:, -1])


def optimise(p, E0, Balg, necro, C0, horizon, maxfev=250, z=None, steady=False,
             recheck_fev=60):
    """Bounded Powell search over the 28-parameter family around E0.  This is a best-found
    reference under stated constraints (fixed catalyst total, fixed algae, necromass and B12,
    no channelling), not an MEP optimum.  Returns (E_best, sigma_best, sigma_start, info) with
    convergence diagnostics: Powell success/message/nfev, a second Powell pass restarted at
    the best point (gain reported), and the budget check."""
    z = np.arange(E0.shape[1]) * p["dz"] + 0.5 * p["dz"] if z is None else z
    cache = {}

    def obj(th):
        key = tuple(np.round(th, 4))
        if key in cache:
            return cache[key]
        E = family(E0, th, z)
        mean, last, _ = sigma_fixed(p, E, Balg, necro, C0, horizon)
        v = -(last if steady else mean)
        cache[key] = v
        return v
    th0 = np.zeros(2 * NR)
    bounds = [(-4, 4)] * NR + [(-10, 10)] * NR
    s_start = -obj(th0)
    res = minimize(obj, th0, method="Powell", bounds=bounds,
                   options=dict(maxfev=maxfev, xtol=1e-2, ftol=1e-4))
    best1 = -res.fun
    res2 = minimize(obj, res.x, method="Powell", bounds=bounds,
                    options=dict(maxfev=recheck_fev, xtol=1e-2, ftol=1e-4))
    xb = res2.x if -res2.fun >= best1 else res.x
    best = max(best1, -res2.fun)
    Eb = family(E0, xb, z)
    n_at_bound = int(np.sum(np.isclose(xb, [b[0] for b in bounds], atol=1e-3)
                            | np.isclose(xb, [b[1] for b in bounds], atol=1e-3)))
    info = dict(powell_success=bool(res.success), powell_message=str(res.message).replace("\t", " "),
                nfev=int(res.nfev + res2.nfev), second_pass_gain=float(best - best1),
                second_pass_rel_gain=float((best - best1) / max(abs(best1), 1e-12)),
                params_at_bound=n_at_bound,
                budget_rel_err=float(abs(Eb.sum() - E0.sum()) / E0.sum()))
    return Eb, best, s_start, info


def random_fields(E0, rng, n=20):
    Ez = E0.sum(0)
    out = []
    for _ in range(n):
        w = rng.dirichlet(np.ones(NR), size=E0.shape[1]).T
        out.append(w * Ez[None])
    return out
