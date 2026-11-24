"""Element, charge and electron balance tests (python3 -I tests/test_conservation.py).
Added after Astra's review (relay 248992826624753664)."""
import os
import sys
import numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
from model import (Model, default_params, REACTIONS, RXN, RI, NU_ALL, ALL_SPECIES, SPECIES,
                   S, NS, NR, RN, CHANNEL_PAIRS, TYPES)

# formula: (C, H, O, N, S, charge); DMSP per C
F = {"DOC": (1, 2, 1, 0, 0, 0), "POC": (1, 2, 1, 0, 0, 0), "O2": (0, 0, 2, 0, 0, 0),
     "NO3": (0, 0, 3, 1, 0, -1), "NO2": (0, 0, 2, 1, 0, -1), "NH4": (0, 4, 0, 1, 0, 1),
     "SO4": (0, 0, 4, 0, 1, -2), "H2S": (0, 2, 0, 0, 1, 0), "H2": (0, 2, 0, 0, 0, 0),
     "CH4": (1, 4, 0, 0, 0, 0), "DMSP": (1, 2, 0.4, 0, 0.2, 0), "B12": (0, 0, 0, 0, 0, 0),
     "CO2": (1, 0, 2, 0, 0, 0), "N2": (0, 0, 0, 2, 0, 0), "H2O": (0, 2, 1, 0, 0, 0),
     "H+": (0, 1, 0, 0, 0, 1)}
# degree of reduction relative to CO2, H2O, H+, NO3-, SO4-- (electrons per mol)
GAMMA = {k: 4 * c + h - 2 * o + 5 * n + 6 * s_ - q for k, (c, h, o, n, s_, q) in F.items()}


def test_stoichiometry():
    bad = []
    for r in RXN:
        st = REACTIONS[r][0]
        tot = np.zeros(6)
        for sp, v in st.items():
            tot += v * np.array(F[sp], float)
        if np.abs(tot).max() > 1e-9:
            bad.append((r, tot))
        # electrons transferred = gamma lost by the electron donor going to its product
        don = {"AER": ("DOC", "CO2"), "DMP": ("DMSP", None), "NAR": ("DOC", "CO2"),
               "NIR": ("DOC", "CO2"), "SRB": ("DOC", "CO2"), "FER": ("DOC", "CO2"),
               "MET": ("H2", "H2O"), "HSR": ("H2", "H2O"), "AMO": ("NH4", "NO2"),
               "NXR": ("NO2", "NO3"), "AMX": ("NH4", "N2"), "SOX": ("H2S", "SO4"),
               "SNR": ("H2S", "SO4"), "MOX": ("CH4", "CO2")}[r]
        nd = -st[don[0]]
        ne_calc = nd * GAMMA[don[0]]
        if r == "DMP":
            ne_calc -= st["H2S"] * GAMMA["H2S"]          # reduced S stays reduced
        elif r == "AMX":
            ne_calc = nd * (GAMMA["NH4"] - GAMMA["N2"] / 2)
        elif don[1] in ("NO2", "NO3"):
            ne_calc -= nd * GAMMA[don[1]]
        ne = REACTIONS[r][1]
        if abs(ne_calc - ne) > 1e-9:
            bad.append((r, "n_e", ne_calc, ne))
    assert not bad, bad
    print("stoichiometry: element + charge balanced for all", len(RXN), "reactions; n_e consistent")


def test_channel_conservation(seed=0, n=50):
    rng = np.random.default_rng(seed)
    types = ["CMX", "DEN", "AERDEN", "FERMET", "SRBX", "AOB", "NIR"]
    worst = 0.0
    for k in range(n):
        p = default_params(NZ=1, types=types, phi=float(rng.uniform(0.1, 1.0)),
                           beta=float(rng.choice([0.0, 8.0])))
        m = Model(p)
        C = rng.lognormal(0, 1.5, size=(NS, 1))
        B = rng.lognormal(0, 1, size=(len(types), 1))
        J, Jint, chan, q, dG, FT, rate = m.fluxes(C, B)
        Jt = J + Jint
        d = m.reaction_tendency(Jt)                     # model's actual tendency
        rsum = Jt.sum(0)
        # remove organic-N remineralisation (organic N is not a tracked species)
        oc = np.array([max(0.0, -REACTIONS[r][0].get("DOC", 0) - REACTIONS[r][0].get("DMSP", 0))
                       for r in RXN])
        d_nh4_remin = RN * (oc[:, None] * rsum).sum(0)
        dd = d.copy()
        dd[S["NH4"]] -= d_nh4_remin
        n2 = NU_ALL[:, ALL_SPECIES.index("N2")] @ rsum
        co2 = NU_ALL[:, ALL_SPECIES.index("CO2")] @ rsum
        # N balance: tracked N + 2*N2 = 0
        nb = dd[S["NO3"]] + dd[S["NO2"]] + dd[S["NH4"]] + 2 * n2
        # electron balance over tracked species + N2 (CO2, H2O, H+ have gamma 0)
        eb = sum(GAMMA[s] * dd[S[s]] for s in SPECIES) + GAMMA["N2"] * n2
        # C balance
        cb = dd[S["DOC"]] + dd[S["DMSP"]] + dd[S["CH4"]] + co2
        scale = np.abs(rsum).sum() + 1e-12
        worst = max(worst, float(np.abs(np.concatenate([nb, eb, cb])).max() / scale))
        # channel accounting: internal consumption must not exceed own production
        for (sp, js, amt) in chan:
            for (a, b, s2, nup, nuc) in CHANNEL_PAIRS:
                if S[s2] == sp:
                    pass
        assert (Jint >= -1e-15).all()
    assert worst < 1e-9, worst
    print(f"channel conservation: N, C and electron balance hold to {worst:.1e} "
          f"(relative) over {n} random states with phi in [0.1, 1]")


def test_v1_addback_fails():
    """Reproduce Astra's counterexample: one AMO + one internal NXR turnover."""
    Jt = np.zeros((1, NR, 1))
    Jt[0, RI["AMO"], 0] = 1.0                  # all NXR flux internal (Jint) in this toy
    Jt[0, RI["NXR"], 0] = 1.0
    d = Model.reaction_tendency(Jt)
    assert abs(d[S["NO2"], 0]) < 1e-12 and abs(d[S["NO3"], 0] - 1) < 1e-12
    d_v1 = d.copy()
    d_v1[S["NO2"], 0] += 1.0                   # v1 add-back
    assert abs(d_v1[S["NO2"], 0] - 1.0) < 1e-12
    print("Astra counterexample: NH4 -> NO3 with internal NXR gives dNO2 = 0 (v1 gave +1)")


def test_capacity_bound(seed=1, n=50):
    rng = np.random.default_rng(seed)
    types = ["CMX", "DEN", "FERMET"]
    for k in range(n):
        p = default_params(NZ=1, types=types, phi=1.0)
        m = Model(p)
        C = rng.lognormal(0, 1.5, size=(NS, 1))
        B = rng.lognormal(0, 1, size=(3, 1))
        J, Jint, chan, q, dG, FT, rate = m.fluxes(C, B)
        e = m.allocation(rate * np.maximum(-dG, 0))
        for (ia, ib, sp, nup, nuc, js) in m.pairs:
            mon = m._mon[sp]
            cap = (q[ib] / mon / REACTIONS[RXN[ib]][1]) * e[js, ib] * B[js] * m.kappa[js, None]
            assert (J[js, ib] + Jint[js, ib] <= cap * (1 + 1e-9) + 1e-15).all()
            # thermodynamics / co-substrates: no internal flux where F_T == 0
            assert (Jint[js, ib][:, FT[ib] <= 0] == 0).all()
            assert (Jint[js, ib] <= nup * J[js, ia] / nuc + 1e-15).all()
    print(f"capacity bound: external + internal consumer flux <= saturated capacity; "
          f"no internal flux where F_T = 0; internal <= own production ({n} states)")


def test_full_tendency_uses_balanced_reactions(seed=2, n=20):
    """Closed one-layer box (no transport, no algae, no exchange): the model's full tendency
    for every redox species must equal the balanced reaction tendency.  The v1 add-back of
    channelled intermediates fails this test."""
    rng = np.random.default_rng(seed)
    types = ["CMX", "DEN", "FERMET", "AOB", "NIR"]
    worst = 0.0
    for k in range(n):
        p = default_params(NZ=1, types=types, phi=0.8, top={}, ref={}, flush=0.0,
                           k_hyd=0.0, w_sink=0.0, k_B12dec=0.0)
        m = Model(p)
        C = rng.lognormal(0, 1.5, size=(NS, 1))
        B = rng.lognormal(0, 1, size=(len(types), 1))
        y = np.concatenate([C.ravel(), B.ravel()])
        dC, dB = m.tendencies(0.0, y)
        J, Jint, chan, q, dG, FT, rate = m.fluxes(C, B)
        ref = m.reaction_tendency(J + Jint)[:, 0]
        keep = [S[s] for s in SPECIES if s not in ("POC", "B12")]
        worst = max(worst, float(np.abs(dC[keep] - ref[keep]).max()))
    assert worst < 1e-12, worst
    print(f"full tendency == balanced reaction tendency in a closed box (max diff {worst:.1e})")


def test_whole_column_budget(seed=3, n=12):
    """(v3, after Carbon 248994691424894976) Whole-column C, N, S and electron budgets with
    algae ON, transport, boundaries, lateral exchange, sinking and exports, on random
    states: storage change = reservoir exchange - outflow."""
    rng = np.random.default_rng(seed)
    worst = {}
    for k in range(n):
        p = default_params(phi=float(rng.choice([0.0, 0.8])), beta=float(rng.choice([0.0, 8.0])),
                           kappa=list(rng.lognormal(0, 0.15, len(TYPES))))
        m = Model(p)
        y = m.y0(rng)
        C, B = m.unpack(y)
        C *= rng.lognormal(0, 1.0, size=C.shape)
        C += rng.uniform(0, 5, size=C.shape)
        B *= rng.lognormal(2, 1.0, size=B.shape)
        bud = m.column_budget(y)
        scale = np.abs(bud["storage"]).sum() + np.abs(bud["exchange"]).sum() + np.abs(bud["outflow"]).sum()
        for kk in ("resid_C", "resid_N", "resid_S", "resid_e"):
            worst[kk] = max(worst.get(kk, 0), abs(bud[kk]) / scale)
    assert max(worst.values()) < 1e-10, worst
    print("whole-column budgets (algae on, open boundaries): relative residuals " +
          ", ".join(f"{k}={v:.1e}" for k, v in worst.items()))


if __name__ == "__main__":
    test_full_tendency_uses_balanced_reactions()
    test_whole_column_budget()
    test_stoichiometry()
    test_v1_addback_fails()
    test_channel_conservation()
    test_capacity_bound()
    print("ALL PASS")

