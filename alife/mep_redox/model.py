"""mep_redox model: 1-D stratified water column with thermodynamically explicit redox
metabolism carried by "metabolic types" (genomes = sets of redox steps), plus phototrophic
eukaryotes (algae, type ALG) as the explicit source of organic C and O2, as metabolite donors
(exudate DOC, DMSP, necromass POC) and as cobalamin (B12) consumers.

Units: depth m, time d, concentrations mmol m^-3 (= uM; B12 in nM), biomass uM C, energies
kJ mol^-1.  Entropy production sigma = sum_r J_r (-dG_r) / T in J K^-1 m^-3 d^-1
(J in mmol m^-3 d^-1, dG in kJ mol^-1).  Column values are integrated over depth (m^-2).

MEP is NOT imposed in the population model: types grow on captured energy minus costs, and
the community is whatever grows.  MEP enters only as an external yardstick (mep.py).
"""
import numpy as np
from scipy.integrate import solve_ivp
from scipy.sparse import lil_matrix

R = 8.314e-3          # kJ mol^-1 K^-1
T = 288.15            # K
RT = R * T

# ---------------------------------------------------------------- chemistry
SPECIES = ["DOC", "POC", "O2", "NO3", "NO2", "NH4", "SO4", "H2S", "H2", "CH4", "DMSP", "B12"]
S = {s: i for i, s in enumerate(SPECIES)}
NS = len(SPECIES)
FIXED = {"CO2": 1e-3, "N2": 5e-4, "H+": 10 ** -7.5, "H2O": 1.0}   # activities
ALL_SPECIES = SPECIES + list(FIXED)
# standard Gibbs energies of formation (kJ/mol, aqueous, 298 K; used at 288 K: approximation)
GF = {"DOC": -152.65,  # CH2O as glucose/6
      "POC": -152.65, "O2": 16.32, "NO3": -111.34, "NO2": -32.22, "NH4": -79.37,
      "SO4": -744.63, "H2S": -27.87, "H2": 17.72, "CH4": -34.39,
      "DMSP": -152.65 + 0.2 * -27.87,   # per C; chosen so that DMP ~ AER (approximation)
      "B12": 0.0,
      "CO2": -386.0, "N2": 18.19, "H2O": -237.18, "H+": 0.0}
RN = 16.0 / 106.0      # NH4 released per organic C oxidised (Redfield remineralisation)

# name: (stoichiometry incl. fixed species, electrons, O2-inhibited?)
REACTIONS = {
    "AER": ({"DOC": -1, "O2": -1, "CO2": 1, "H2O": 1}, 4, False),
    # DMSP per C (C5H10O2S / 5 = C H2 O0.4 S0.2), balanced in C, H, O, S, charge
    "DMP": ({"DMSP": -1, "O2": -1.2, "CO2": 1, "H2O": 0.8, "H2S": 0.2}, 4.8, False),
    "NAR": ({"DOC": -1, "NO3": -2, "NO2": 2, "CO2": 1, "H2O": 1}, 4, True),
    "NIR": ({"DOC": -3, "NO2": -4, "H+": -4, "CO2": 3, "N2": 2, "H2O": 5}, 12, True),
    "SRB": ({"DOC": -2, "SO4": -1, "H+": -2, "H2S": 1, "CO2": 2, "H2O": 2}, 8, True),
    "FER": ({"DOC": -1, "H2O": -1, "CO2": 1, "H2": 2}, 4, True),
    "MET": ({"H2": -4, "CO2": -1, "CH4": 1, "H2O": 2}, 8, True),
    "HSR": ({"H2": -4, "SO4": -1, "H+": -2, "H2S": 1, "H2O": 4}, 8, True),
    "AMO": ({"NH4": -1, "O2": -1.5, "NO2": 1, "H+": 2, "H2O": 1}, 6, False),
    "NXR": ({"NO2": -1, "O2": -0.5, "NO3": 1}, 2, False),
    "AMX": ({"NH4": -1, "NO2": -1, "N2": 1, "H2O": 2}, 3, True),
    "SOX": ({"H2S": -1, "O2": -2, "SO4": 1, "H+": 2}, 8, False),
    "SNR": ({"H2S": -1, "NO3": -4, "SO4": 1, "NO2": 4, "H+": 2}, 8, True),
    "MOX": ({"CH4": -1, "O2": -2, "CO2": 1, "H2O": 2}, 8, False),
}
RXN = list(REACTIONS)
RI = {r: i for i, r in enumerate(RXN)}
NR = len(RXN)
NU_ALL = np.zeros((NR, len(ALL_SPECIES)))
for r, (st, ne, inh) in REACTIONS.items():
    for sp, v in st.items():
        NU_ALL[RI[r], ALL_SPECIES.index(sp)] = v
NU = NU_ALL[:, :NS].copy()                      # tracked-species stoichiometry
for r in RXN:                                   # organic-N remineralisation
    oc = -REACTIONS[r][0].get("DOC", 0)        # organic N only in DOC/POC (DMSP has none)
    NU[RI[r], S["NH4"]] += RN * max(0.0, oc)
NE = np.array([REACTIONS[r][1] for r in RXN], float)
O2INH = np.array([REACTIONS[r][2] for r in RXN])
DG0 = NU_ALL @ np.array([GF[s] for s in ALL_SPECIES])
LN_FIXED = np.log(np.array([FIXED[s] for s in FIXED]))
KS = {"DOC": 10.0, "O2": 1.0, "NO3": 1.0, "NO2": 1.0, "NH4": 1.0, "SO4": 100.0,
      "H2S": 1.0, "H2": 0.05, "CH4": 1.0, "DMSP": 0.5, "B12": 0.01}
SUBSTR = [[S[s] for s in SPECIES if NU_ALL[i, S[s]] < 0] for i in range(NR)]

# photosynthesis bookkeeping (algae).  Photon free energy following Vallino & Huber 2018's
# construction (photon Gibbs energy times a thermodynamic efficiency); the numbers are ours.
N_PHOTON = 10.0                 # photons per C fixed
MU_PHOTON = 0.75 * 176.0        # kJ per mol photons (680 nm, eta_I = 0.75)

# ------------------------------------------------------------- metabolic types
TYPES = {
    "AER": ["AER"], "ROS": ["AER", "DMP"],            # ROS: Roseobacter-like DMSP degrader
    "NAR": ["NAR"], "NIR": ["NIR"], "DEN": ["NAR", "NIR"],
    "AERNAR": ["AER", "NAR"], "AERDEN": ["AER", "NAR", "NIR"],
    "SRB": ["SRB"], "SRBH": ["HSR"], "SRBX": ["SRB", "HSR"],
    "FER": ["FER"], "MET": ["MET"], "FERMET": ["FER", "MET"],
    "AOB": ["AMO"], "NOB": ["NXR"], "CMX": ["AMO", "NXR"], "AMX": ["AMX"],
    "SOB": ["SOX"], "SNR": ["SNR"], "SOBX": ["SOX", "SNR"], "MOB": ["MOX"],
    "ALG": [],                                        # phototrophic eukaryote (special)
}
B12_PRODUCERS = {"ROS"}   # cobamide maker (leaky public good); one supplier type so that
# pool restrictions on split/stacked genomes do not change who makes B12
CHANNEL_PAIRS = [("NAR", "NIR", "NO2", 2.0, 4.0), ("AMO", "NXR", "NO2", 1.0, 1.0),
                 ("FER", "MET", "H2", 2.0, 4.0), ("FER", "HSR", "H2", 2.0, 4.0)]

# KO mapping for synthetic observations (realistic ambiguities kept: narGH = nxrAB,
# amoABC = pmoABC; sat/aprA/dsrAB shared by sulfate reducers on DOC and on H2).
STEP_KOS = {
    "AER": ["K02274", "K02276", "K00404"],            # coxA coxC ccoN
    "DMP": ["K17486"],                                # dmdA
    "NAR": ["K00370", "K00371", "K00374"],            # narG narH narI
    "NIR": ["K00368", "K04561", "K00376"],            # nirK norB nosZ
    "SRB": ["K00958", "K00394", "K11180", "K11181"],  # sat aprA dsrA dsrB
    "FER": ["K00532", "K00656"],                      # hydA pflD
    "MET": ["K00399", "K00401", "K00402"],            # mcrA mcrB mcrG
    "HSR": ["K00958", "K00394", "K11180", "K11181", "K06282"],   # + hyaA
    "AMO": ["K10944", "K10945", "K10946", "K10535"],  # amoABC hao
    "NXR": ["K00370", "K00371"],                      # nxrAB (same KOs as narGH)
    "AMX": ["K20932", "K20935"],                      # hzsA hdh
    "SOX": ["K17224", "K17226", "K17227", "K17218"],  # soxB soxY soxZ sqr
    "SNR": ["K17218", "K00370", "K00371"],            # sqr + narGH
    "MOX": ["K10944", "K10945", "K10946", "K16157"],  # pmoABC mmoX
    "PHO": ["K02703", "K02706", "K01601"],            # psbA psbD rbcL
    "B12syn": ["K02232", "K02230", "K02233", "K00798"],  # cobQ cobN cobS cobA/btuR
    "B12dep": ["K00548"],                             # metH (B12-dependent)
}


def default_params(**kw):
    p = dict(
        NZ=30, dz=1.0, mld=5.0, Kz_mix=20.0, Kz_deep=0.3,   # m, m^2/d
        cell_mix=1.0,                  # multiplier on Kz for cells (dispersal)
        # algae (phototrophic eukaryotes)
        I0=1.0, k_w=0.3, k_chl=0.01,   # relative surface light; water and self-shading 1/m
        Ik=0.2, mu_alg=2.0, eps_alg=0.5, f_ex=0.2, f_dmsp=0.05,
        alg_sink=0.1, alg_mort=0.1, alg_d2=0.002, K_B12=0.01, q_B12=1e-3, K_Nalg=0.5,
        force_amp=0.0, force_period=30.0,
        # B12 supply by bacteria (leaky public good)
        y_B12=2e-2, c_B12=0.02, k_B12dec=0.05,
        # physics
        w_sink=1.0, k_hyd=0.1, f_necro=0.5,
        flush=0.002,
        top={"O2": 280.0, "NO3": 5.0, "NH4": 0.0, "SO4": 250.0, "H2S": 0.0, "CH4": 0.0,
             "H2": 0.0, "NO2": 0.0, "B12": 0.002},
        ref={"SO4": 250.0},
        k0=10.0, Gmin=1.0, chi=1.0, K_O2inh=1.0,
        Y_E=1e-3, yield_mode="energy", Y_e=0.06,
        maint=0.01, mort=0.02, d2=0.002, immig=1e-4,
        # ---- switchable kinetic symmetry breakers ----
        alpha=1.0, c_step=0.0, beta=8.0, phi=0.0, K_I_NO2=np.inf,
        types=list(TYPES), kappa=None,
    )
    p.update(kw)
    return p


class Model:
    fixed = False

    def __init__(self, p):
        self.p = p
        self.types = list(p["types"])
        self.NT = len(self.types)
        NZ = p["NZ"]
        self.z = (np.arange(NZ) + 0.5) * p["dz"]
        self.M = np.zeros((self.NT, NR))
        for j, t in enumerate(self.types):
            for r in (TYPES[t] if t in TYPES else [t]):
                self.M[j, RI[r]] = 1.0
        self.n = self.M.sum(1)
        self.nsafe = np.maximum(self.n, 1.0)
        self.ialg = self.types.index("ALG") if "ALG" in self.types else None
        self.b12p = np.array([t in B12_PRODUCERS for t in self.types], float)
        self.kappa = np.ones(self.NT) if p["kappa"] is None else np.asarray(p["kappa"], float)
        Kz = np.where(self.z < p["mld"], p["Kz_mix"], p["Kz_deep"])
        self.Kf = np.minimum(Kz[:-1], Kz[1:])
        self.Ktop = Kz[0]
        self.top_mask = np.zeros(NS, bool)
        self.top_val = np.zeros(NS)
        for s, v in p["top"].items():
            self.top_mask[S[s]] = True
            self.top_val[S[s]] = v
        self.ref = np.zeros(NS)
        for s, v in p["ref"].items():
            self.ref[S[s]] = v
        self.ks = np.array([KS.get(s, 1.0) for s in SPECIES])
        self.pairs = []
        for (a, b, sp, nup, nuc) in CHANNEL_PAIRS:
            js = [j for j in range(self.NT) if self.M[j, RI[a]] and self.M[j, RI[b]]]
            if js:
                self.pairs.append((RI[a], RI[b], S[sp], nup, nuc, np.array(js)))
        self.NY = (NS + self.NT) * NZ

    # ------------------------------------------------------------ helpers
    def unpack(self, y):
        NZ = self.p["NZ"]
        return y[:NS * NZ].reshape(NS, NZ), y[NS * NZ:].reshape(self.NT, NZ)

    def diffuse(self, X, Kmul=1.0, top=None):
        dz = self.p["dz"]
        F = -Kmul * self.Kf * (X[:, 1:] - X[:, :-1]) / dz
        d = np.zeros_like(X)
        d[:, :-1] -= F / dz
        d[:, 1:] += F / dz
        if top is not None:
            mask, val = top
            d[mask, 0] += Kmul * self.Ktop * (val[mask] - X[mask, 0]) / (dz * dz / 2)
        return d

    def thermo(self, C):
        a = np.maximum(C, 1e-9) * 1e-6
        lnA = np.vstack([np.log(a), np.repeat(LN_FIXED[:, None], C.shape[1], 1)])
        dG = DG0[:, None] + RT * (NU_ALL @ lnA)
        FT = np.clip(1.0 - np.exp((dG / NE[:, None] + self.p["Gmin"]) / (self.p["chi"] * RT)),
                     0.0, 1.0)
        return dG, FT

    def specific_rates(self, C):
        p = self.p
        Cp = np.maximum(C, 0.0)
        dG, FT = self.thermo(C)
        mon = Cp / (Cp + self.ks[:, None])
        FK = np.ones((NR, C.shape[1]))
        for i in range(NR):
            for s in SUBSTR[i]:
                FK[i] *= mon[s]
        inh = np.where(O2INH[:, None], p["K_O2inh"] / (p["K_O2inh"] + Cp[S["O2"]][None, :]), 1.0)
        self._mon = mon                                  # kept for channelling bounds
        return p["k0"] * FK * FT * inh, dG, FT

    def light(self, t, B):
        p = self.p
        f = max(1.0 + p["force_amp"] * np.sin(2 * np.pi * t / p["force_period"]), 0.0)
        if self.ialg is None:
            shade = np.zeros(p["NZ"])
        else:
            ba = np.maximum(B[self.ialg], 0) * p["dz"]
            shade = p["k_chl"] * (np.cumsum(ba) - 0.5 * ba)
        return p["I0"] * f * np.exp(-p["k_w"] * self.z - shade)

    def allocation(self, power):
        """Effective enzyme per unit biomass, e[type, reaction, z]."""
        p = self.p
        if p["beta"] > 0:
            w = (power[None] + 1e-6) * self.M[:, :, None]
            w = w / (w.max(1, keepdims=True) + 1e-30)
            w = np.where(self.M[:, :, None] > 0, w ** p["beta"], 0.0)
            Om = w / (w.sum(1, keepdims=True) + 1e-30)
        else:
            Om = np.broadcast_to((self.M / self.nsafe[:, None])[:, :, None],
                                 (self.NT, NR, power.shape[1]))
        return Om * (self.nsafe ** (1.0 - p["alpha"]))[:, None, None]

    def fluxes(self, C, B, t=0.0):
        """Per type x reaction x depth fluxes (mmol rxn m^-3 d^-1) and pieces."""
        p = self.p
        q, dG, FT = self.specific_rates(C)
        rate = q / NE[:, None]
        e = self.allocation(rate * np.maximum(-dG, 0.0))
        Bp = np.maximum(B, 0.0)
        J = rate[None] * e * (Bp * self.kappa[:, None])[:, None, :]
        Jint = np.zeros_like(J)
        chan = []
        if p["phi"] > 0:
            # Channelling (revised after Astra 248992826624753664): a stacked genome's
            # consumer step may take a fraction phi of the intermediate made by its own
            # producer step directly, without the intermediate's Monod limit.  It still obeys
            # (i) the co-substrate Monod terms, F_T (incl. G_min) and O2 inhibition of the
            # consumer step, evaluated at bulk activities, and (ii) the same enzyme
            # capacity: external + internal consumer flux <= capacity at saturating
            # intermediate.
            mon = self._mon
            for (ia, ib, sp, nup, nuc, js) in self.pairs:
                prod = nup * J[js, ia]                      # intermediate made (js, NZ)
                q_sat = q[ib] / np.maximum(mon[sp], 1e-12)  # q with intermediate saturated
                q_sat = np.where(mon[sp] > 1e-12, q_sat, 0.0)
                cap = (q_sat / NE[ib])[None, :] * e[js, ib] * Bp[js] * self.kappa[js, None]
                room = np.maximum(cap - J[js, ib], 0.0)
                want = p["phi"] * prod / nuc
                ji = want * room / (want + room + 1e-30)
                Jint[js, ib] += ji
                chan.append((sp, js, ji * nuc))
        return J, Jint, chan, q, dG, FT, rate

    def algae(self, t, C, B):
        """Gross C fixation by algae (mmol C m^-3 d^-1) and B12 limitation factor."""
        p = self.p
        NZ = p["NZ"]
        if self.ialg is None:
            return np.zeros(NZ), np.ones(NZ)
        I = self.light(t, B)
        b12 = np.maximum(C[S["B12"]], 0)
        fb = b12 / (b12 + p["K_B12"])
        din = np.maximum(C[S["NH4"]], 0) + np.maximum(C[S["NO3"]], 0)
        fb = fb * din / (din + p["K_Nalg"])          # (v3) N limitation of algae
        Ba = np.maximum(B[self.ialg], 0)
        G = p["mu_alg"] / p["eps_alg"] * I / (I + p["Ik"]) * fb * Ba
        return G, fb

    @staticmethod
    def reaction_tendency(Jtot):
        """Tracked-species tendency from catabolic reactions (incl. channelled internal
        turnover), mmol m^-3 d^-1.  Tested in tests/test_conservation.py."""
        return NU.T @ Jtot.sum(0)

    def tendencies(self, t, y, parts=False):
        p = self.p
        C, B = self.unpack(y)
        J, Jint, chan, q, dG, FT, rate = self.fluxes(C, B, t)
        Jtot = J + Jint
        dC = self.reaction_tendency(Jtot)
        # NU.T @ Jtot already nets the channelled intermediate: the producer step adds it,
        # the internal consumer step removes it.  (v1 added it back here: a conservation
        # error found by Astra, relay 248992826624753664.)
        Bp = np.maximum(B, 0)
        # ---- algae: gross fixation G -> eps to biomass, f_ex exudate DOC, f_dmsp DMSP,
        #      the remainder respired by the alga itself (uses its own O2)
        G, fb = self.algae(t, C, B)
        grow_alg = p["eps_alg"] * G
        ex = p["f_ex"] * G
        dm = p["f_dmsp"] * G
        resp = G - grow_alg - ex - dm
        # (v3, whole-column budget after Carbon 248994691424894976) organic N at Redfield in
        # algal biomass and exudate DOC, taken up as NH4 first then NO3 (NO3 assimilation
        # releases 2 O2 per N); DMSP (C H2 O0.4 S0.2 per C) is made from CO2 + 0.2 SO4 and
        # releases 1.6 O2 per C.
        nupt = RN * (grow_alg + ex)
        a4 = np.maximum(C[S["NH4"]], 0)
        a3 = np.maximum(C[S["NO3"]], 0)
        u4 = nupt * a4 / (a4 + a3 + 1e-30)
        u3 = nupt - u4
        dC[S["DOC"]] += ex
        dC[S["DMSP"]] += dm
        dC[S["NH4"]] -= u4
        dC[S["NO3"]] -= u3
        dC[S["SO4"]] -= 0.2 * dm
        dC[S["O2"]] += (G - dm) - resp + 1.6 * dm + 2.0 * u3
        dC[S["B12"]] -= p["q_B12"] * grow_alg
        # ---- POC dynamics
        hyd = p["k_hyd"] * np.maximum(C[S["POC"]], 0)
        dC[S["POC"]] -= hyd
        dC[S["DOC"]] += hyd
        sink = p["w_sink"] * np.maximum(C[S["POC"]], 0) / p["dz"]
        dC[S["POC"]] -= sink
        dC[S["POC"], 1:] += sink[:-1]
        dC[S["DOC"], -1] += sink[-1]                 # sediment return as DOC
        dC += self.diffuse(C, 1.0, (self.top_mask, self.top_val))
        dC += p["flush"] * (self.ref[:, None] - C)
        dC[S["B12"]] -= p["k_B12dec"] * np.maximum(C[S["B12"]], 0)
        # ---- biomass (needed for necromass even when biomass is frozen)
        if p["yield_mode"] == "energy":
            gain = p["Y_E"] * (Jtot * np.maximum(-dG, 0)[None]).sum(1)
        else:
            gain = p["Y_e"] * (Jtot * NE[None, :, None]).sum(1)
        if np.isfinite(p["K_I_NO2"]):
            gain = gain * (p["K_I_NO2"] / (p["K_I_NO2"] + np.maximum(C[S["NO2"]], 0)))[None]
        b12rel = p["y_B12"] * self.b12p[:, None] * gain
        gain = gain * (1.0 - p["c_B12"] * self.b12p[:, None])
        mortB = (p["mort"] + p["c_step"] * (self.n - 1))[:, None] * Bp + p["d2"] * Bp * Bp
        loss = p["maint"] * Bp + mortB
        if self.ialg is not None:
            ia = self.ialg
            gain[ia] = grow_alg
            b12rel[ia] = 0.0
            mortB[ia] = p["alg_mort"] * Bp[ia] + p["alg_d2"] * Bp[ia] ** 2
            loss[ia] = mortB[ia]
        if not self.fixed:
            dC[S["B12"]] += b12rel.sum(0)
            # (v3) chemotroph biomass is catalyst only (grows on captured energy without
            # drawing C, N or electrons), so its necromass is NOT returned to POC; only
            # algal necromass (a real C/N/electron flow) is.
            if self.ialg is not None:
                dC[S["POC"]] += p["f_necro"] * mortB[self.ialg]
        else:
            # frozen biomass: necromass flux prescribed (from the population state it copies)
            dC[S["POC"]] += self.necro
        if parts:
            return dict(J=J, Jint=Jint, chan=chan, dG=dG, FT=FT, rate=rate, G=G, fb=fb,
                        ex=ex, dm=dm, resp=resp, u4=u4, u3=u3, hyd=hyd, gain=gain, mortB=mortB,
                        b12rel=b12rel, loss=loss)
        if self.fixed:
            return dC.ravel(), None
        mix = self.diffuse(B, p["cell_mix"])
        if self.ialg is not None and p["cell_mix"] != 1.0:
            mix[self.ialg] = self.diffuse(B[self.ialg:self.ialg + 1], 1.0)[0]   # algae: water Kz
        dB = gain - loss + p["immig"] + mix - p["flush"] * B
        if self.ialg is not None:
            sa = p["alg_sink"] * Bp[self.ialg] / p["dz"]
            dB[self.ialg] -= sa
            dB[self.ialg, 1:] += sa[:-1]
        if getattr(self, "knockout", None) is not None:
            dB[self.knockout] = -B[self.knockout]      # removed types decay to zero, no immigration
        return dC.ravel(), dB.ravel()

    def f(self, t, y):
        dC, dB = self.tendencies(t, y)
        return dC if dB is None else np.concatenate([dC, dB])

    def sparsity(self):
        NZ = self.p["NZ"]
        nv = NS + (0 if self.fixed else self.NT)
        Sp = lil_matrix((nv * NZ, nv * NZ), dtype=np.int8)
        for z in range(NZ):
            rows = [v * NZ + z for v in range(nv)]
            for r_ in rows:
                for c_ in rows:
                    Sp[r_, c_] = 1
            for v in range(nv):
                for dzz in (-1, 1):
                    if 0 <= z + dzz < NZ:
                        Sp[v * NZ + z, v * NZ + z + dzz] = 1
        if self.ialg is not None and not self.fixed:   # self-shading couples algae downward
            a0 = (NS + self.ialg) * NZ
            for z in range(NZ):
                for z2 in range(z):
                    for v in range(nv):
                        Sp[v * NZ + z, a0 + z2] = 1
        return Sp.tocsr()

    def y0(self, rng=None, B0=1e-2):
        NZ = self.p["NZ"]
        C0 = np.zeros((NS, NZ))
        C0[S["O2"]] = 280.0 * np.exp(-self.z / 5)
        C0[S["SO4"]] = 250.0
        C0[S["NO3"]] = 5.0
        C0[S["NH4"]] = 1.0
        C0[S["B12"]] = 0.01
        B = np.full((self.NT, NZ), B0)
        if rng is not None:
            B = B * rng.lognormal(0, 0.5, size=B.shape)
        if self.ialg is not None:
            B[self.ialg] = 1.0
        return np.concatenate([C0.ravel(), B.ravel()])

    def run(self, y0, t_end, t_eval=None, rtol=1e-5, atol=1e-7, max_step=np.inf):
        if not hasattr(self, "_sp"):
            self._sp = self.sparsity()
        sol = solve_ivp(self.f, (0, t_end), y0, method="BDF", t_eval=t_eval,
                        jac_sparsity=self._sp, rtol=rtol, atol=atol, max_step=max_step)
        if not sol.success:
            raise RuntimeError(sol.message)
        return sol

    # ------------------------------------------------------------ diagnostics
    def sigma_parts(self, y, t=0.0):
        """Column entropy production (J K^-1 m^-2 d^-1): chemotrophic reactions and algal
        photon dissipation (photon free energy captured minus chemical potential stored).
        Light absorbed by water is not counted.  Also returns per-depth chemotroph profile
        and per-type chemotrophic sigma."""
        C, B = self.unpack(y)
        J, Jint, chan, q, dG, FT, rate = self.fluxes(C, B, t)
        dGp = np.maximum(-dG, 0)
        per_type = ((J + Jint) * dGp[None]).sum((1, 2)) * self.p["dz"] / T
        prof = ((J + Jint).sum(0) * dGp).sum(0) / T
        G, fb = self.algae(t, C, B)
        stored = dGp[RI["AER"]]
        alg = G * np.maximum(N_PHOTON * MU_PHOTON - stored, 0) / T
        dz = self.p["dz"]
        return float(prof.sum() * dz), float(alg.sum() * dz), prof, per_type

    def npp(self, y, t=0.0):
        """Column net primary production of algae (gross fixation minus algal respiration),
        mmol C m^-2 d^-1: the organic-C (and O2) supply that chemotrophs can dissipate."""
        C, B = self.unpack(y)
        G, fb = self.algae(t, C, B)
        p = self.p
        return float((G * (p["eps_alg"] + p["f_ex"] + p["f_dmsp"])).sum() * p["dz"])

    def column_budget(self, y, t=0.0):
        """Whole-column budgets (v3): storage change vs boundary exchange and exchange with
        the fixed reservoirs (CO2, N2), for C, N, S and electrons (degree of reduction
        relative to CO2, H2O, H+, NO3-, SO4--).  Returns dict of residuals (mmol m^-2 d^-1)
        and the terms.  Algal biomass and DOC/POC carry Redfield organic N (gamma 8 per N).
        Chemotroph biomass is catalyst only and is excluded."""
        p = self.p
        dz = p["dz"]
        C, B = self.unpack(y)
        dC, dB = self.tendencies(t, y)
        dC = dC.reshape(NS, -1)
        P = self.tendencies(t, y, parts=True)
        Jt = (P["J"] + P["Jint"]).sum(0)
        dBa = dB.reshape(self.NT, -1)[self.ialg] if (dB is not None and self.ialg is not None) \
            else np.zeros(p["NZ"])
        content = {  # per unit of each pool: (C, N, S, gamma)
            "DOC": (1, RN, 0, 4 + 8 * RN), "POC": (1, RN, 0, 4 + 8 * RN), "ALG": (1, RN, 0, 4 + 8 * RN),
            "O2": (0, 0, 0, -4), "NO3": (0, 1, 0, 0), "NO2": (0, 1, 0, 2), "NH4": (0, 1, 0, 8),
            "SO4": (0, 0, 1, 0), "H2S": (0, 0, 1, 8), "H2": (0, 0, 0, 2), "CH4": (1, 0, 0, 8),
            "DMSP": (1, 0, 0.2, 6.4)}
        stor = np.zeros(4)
        for sp, v in content.items():
            rate = dBa if sp == "ALG" else dC[S[sp]]
            stor += np.array(v) * rate.sum() * dz
        # exchanges with fixed reservoirs
        co2 = (NU_ALL[:, ALL_SPECIES.index("CO2")] @ Jt).sum() * dz - (P["G"] - P["resp"]).sum() * dz
        n2 = (NU_ALL[:, ALL_SPECIES.index("N2")] @ Jt).sum() * dz
        # boundary and loss terms (all expressed as content leaving the column)
        out = np.zeros(4)
        Ktop = self.Ktop
        for sp, v in content.items():
            if sp == "ALG":
                continue
            s_ = S[sp]
            if self.top_mask[s_]:
                out -= np.array(v) * Ktop * (self.top_val[s_] - C[s_, 0]) / (dz / 2)
            out -= np.array(v) * p["flush"] * (self.ref[s_] - C[s_]).sum() * dz
        Ba = np.maximum(B[self.ialg], 0) if self.ialg is not None else np.zeros(p["NZ"])
        alg_out = 0.0
        if self.ialg is not None:
            alg_out = ((1 - p["f_necro"]) * P["mortB"][self.ialg].sum() * dz
                       + p["alg_sink"] * Ba[-1] + p["flush"] * B[self.ialg].sum() * dz
                       - p["immig"] * p["NZ"] * dz)
        out += np.array(content["ALG"]) * alg_out
        # B12 and chemotroph catalyst are outside the budget
        src = np.array([-co2, -2 * n2, 0.0, -10 * n2])   # CO2 fixed in, N2 lost
        resid = stor - (src - out)
        return dict(resid_C=resid[0], resid_N=resid[1], resid_S=resid[2], resid_e=resid[3],
                    storage=stor, exchange=src, outflow=out)

    def sigma_column(self, y, t=0.0):
        return self.sigma_parts(y, t)[0]


class FixedEnzyme(Model):
    """Chemistry-only model with a prescribed enzyme field E[r, z] and prescribed algal
    biomass: the 'song' without the 'singers'.  greedy_beta: per-depth catalyst total is
    re-allocated at every instant to the reactions of highest instantaneous power
    (Vallino's short-horizon limit).  B12 is held at its supplied profile (not a control)."""
    fixed = True

    def __init__(self, p, E, Balg, necro, greedy_beta=None):
        q = dict(p)
        q.update(types=list(RXN) + ["ALG"], kappa=None, alpha=1.0, phi=0.0,
                 beta=(greedy_beta if greedy_beta else 0.0), c_step=0.0)
        super().__init__(q)
        E = np.asarray(E, float)
        if greedy_beta is not None:
            self.types = ["OMNI", "ALG"]
            self.NT = 2
            self.M = np.vstack([np.ones(NR), np.zeros(NR)])
            self.n = self.M.sum(1)
            self.nsafe = np.ones(2)
            self.kappa = np.ones(2)
            self.pairs = []
            self.ialg = 1
            self.b12p = np.zeros(2)
            self.Bfix = np.vstack([E.sum(0), Balg])
        else:
            self.Bfix = np.vstack([E, Balg])
        self.necro = np.asarray(necro, float)
        self.NY = NS * self.p["NZ"]

    def unpack(self, y):
        return y.reshape(NS, self.p["NZ"]), self.Bfix

    def f(self, t, y):
        dC, _ = self.tendencies(t, y)
        dC = dC.reshape(NS, -1)
        dC[S["B12"]] = 0.0
        return dC.ravel()


def enzyme_field(model, C, B):
    """Effective enzyme per reaction and depth implied by a population state."""
    q, dG, FT = model.specific_rates(C)
    rate = q / NE[:, None]
    e = model.allocation(rate * np.maximum(-dG, 0.0))
    return (e * (np.maximum(B, 0) * model.kappa[:, None])[:, None, :]).sum(0)
