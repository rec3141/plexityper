"""bq_lattice world: individual-based 2-D lattice with diffusing metabolites (numpy only).

Mechanisms (nothing about who depends on whom is coded; only chemistry, cost and loss):
- One cell per site, periodic L x L lattice.  Water sites mix cells (pairwise swaps) and have
  diffusion D_w; particle sites (discs, fraction phi) do not mix cells and have D_p = D_w * pdiff.
- Genome = bitmask over K functions.  Every cell needs every product each step (demand theta=1).
  * leaky compatible  L0..L4 : producer exports a fraction `leak` of yield Y into a diffusing field
  * OX  (leaky)        : O2-producing route; product also exported; emits O2 into an O2 field
  * ANX (leaky)        : O2-sensitive route; efficiency (1 - iota*has_OX) / (1 + iota*O2/K_O)
  * private P0..P3     : product is never exported; a non-carrier can grow only inside a
                         static abiotic patch of that product (coverage priv_cov)
  Ambient abiotic supply amb per site-step for every leaky product (the "medium").
- Growth prob per step g = rmax * (1-c)^n_carried * min_j supply_j/theta (Liebig).
  Birth into a random empty Moore neighbour; death prob mu.  Gene loss per gene per division
  u_loss (Black Queen: losing saves cost).  HGT: per cell-step prob h_hgt, gain one function a
  random neighbour carries.  iota ramps from 0 to its target over the first `ramp` steps
  (an oxygenation event), so the ancestor (a full prototroph) is viable at t=0.
- Provenance: in tagged windows each leaky field is split into source tags (producer genotype,
  ENV=ambient, PRE=content at window start, OTHER=minor producers).  Uptake removes tags in
  proportion, so every unit a recipient takes is attributed to a source genotype.
"""
import numpy as np

KL = 5                          # leaky compatible
OX, ANX = 5, 6                  # incompatible pair (both leaky)
PRIV = [7, 8, 9, 10]            # private
K = 11
FL = 7                          # leaky fields: 0..6
LEAKY = list(range(FL))
NAMES = [f"L{i}" for i in range(KL)] + ["OX", "ANX"] + [f"P{i}" for i in range(4)]
KAR = np.arange(K)
POW = (1 << KAR).astype(np.int64)
FULL = (1 << K) - 1
LEAKYMASK = (1 << FL) - 1
EUK = 1 << K                     # flag bit: phototrophic eukaryote (fixed genome, never evolves)
EUKGENO = EUK | (1 << 0) | (1 << 1) | (1 << OX)   # exports L0, L1 and OX (organic C + O2)
NG = 1 << (K + 1)                # genotype index space
B12 = 4                          # L4 plays the cobalamin role: the eukaryote requires it

DEFAULTS = dict(L=64, T=5000, ramp=1500, leak=0.5, cost=0.03, D=0.3, pdiff=0.1, iota=1.0,
                phi=0.2, Y=4.0, theta=1.0, rmax=0.5, mu=0.05, u_loss=1e-3, h_hgt=2e-4,
                amb=0.02, decay=0.03, priv_cov=0.95, o_y=1.0, o_resp=0.2, o_decay=0.02,
                o_bg=0.5, o_relax=0.1, K_O=0.5, o_dmult=2.0, p_mix=0.5, init_occ=0.5,
                priv_bg=0.0, euk_frac=0.1, Y_e=8.0, rmax_e=0.4, mu_e=0.04, shade_n=3.0,
                light_p=0.3, resp_p_mult=3.0, stack_discount=0.0)


def diff_ops(Dmap, maxc=0.2):
    n = int(np.ceil(Dmap.max() / maxc)) if Dmap.max() > 0 else 0
    if n == 0:
        return None
    d = Dmap / n
    Ex = np.minimum(d, np.roll(d, -1, 0)).astype(np.float32)
    Ey = np.minimum(d, np.roll(d, -1, 1)).astype(np.float32)
    return Ex, Ey, n


def diffuse(A, ops):
    if ops is None:
        return
    Ex, Ey, n = ops
    for _ in range(n):
        fx = Ex * (np.roll(A, -1, -2) - A)
        fy = Ey * (np.roll(A, -1, -1) - A)
        A += fx - np.roll(fx, 1, -2) + fy - np.roll(fy, 1, -1)


def smooth_patches(rng, L, cov, scale=6):
    z = rng.standard_normal((L, L))
    k = np.fft.fftfreq(L)
    kx, ky = np.meshgrid(k, k, indexing="ij")
    f = np.real(np.fft.ifft2(np.fft.fft2(z) * np.exp(-(kx ** 2 + ky ** 2) * (scale * np.pi) ** 2)))
    return (f > np.quantile(f, 1 - cov)).astype(np.float32)


class World:
    def __init__(self, p, rng):
        self.p = p = {**DEFAULTS, **p}
        self.rng = rng
        L = p["L"]
        self.L, self.N = L, L * L
        # particles: discs until fraction phi
        part = np.zeros((L, L), bool)
        xx, yy = np.meshgrid(np.arange(L), np.arange(L), indexing="ij")
        while part.mean() < p["phi"]:
            cx, cy, r = rng.integers(L), rng.integers(L), rng.uniform(3.5, 7.0)
            dx = np.minimum(abs(xx - cx), L - abs(xx - cx))
            dy = np.minimum(abs(yy - cy), L - abs(yy - cy))
            part |= dx ** 2 + dy ** 2 <= r * r
        self.part2 = part
        self.part = part.ravel()
        Dmap = np.where(part, p["D"] * p["pdiff"], p["D"]).astype(np.float32)
        self.dops = diff_ops(Dmap)
        self.oops = diff_ops(Dmap * p["o_dmult"])
        # private patch maps (abiotic supply of private products)
        self.priv_map = np.stack([np.maximum(smooth_patches(rng, L, p["priv_cov"], 8).ravel(),
                                             p["priv_bg"]) for _ in PRIV], 1)
        # neighbours (Moore, periodic)
        idx = np.arange(self.N)
        x, y = idx // L, idx % L
        self.nbr = np.stack([((x + dx) % L) * L + (y + dy) % L
                             for dx in (-1, 0, 1) for dy in (-1, 0, 1) if dx or dy], 1)
        # checkerboard swap pairs for water mixing
        self.pairs = []
        for ax in (0, 1):
            for par in (0, 1):
                a = idx[(x if ax == 0 else y) % 2 == par]
                ax_, ay_ = a // L, a % L
                b = (((ax_ + 1) % L) * L + ay_) if ax == 0 else (ax_ * L + (ay_ + 1) % L)
                keep = ~self.part[a] & ~self.part[b]
                self.pairs.append((a[keep], b[keep]))
        # state
        self.geno = np.full(self.N, -1, np.int32)
        occ = rng.random(self.N) < p["init_occ"]
        self.geno[occ] = FULL
        self.geno[occ & (rng.random(self.N) < p["euk_frac"])] = EUKGENO
        self.M = np.zeros((FL, L, L), np.float32)
        self.M2 = self.M.reshape(FL, -1)
        self.O2 = np.full((L, L), p["o_bg"], np.float32)
        self.O22 = self.O2.ravel()
        self.t = 0
        self.tags = None
        self.evolve = True
        self.track_growth = False
        self.ghost = None   # counterfactual: frozen non-exporting consumers (bool per site)

    # ---------------- provenance tagging ----------------
    def start_tags(self, ntop=12, only=None):
        """Split each leaky field into tags: 0=PRE, 1=ENV, 2=OTHER, 3.. = top producer genotypes
        (restricted to genotypes in `only` if given, e.g. the kept set of a counterfactual)."""
        NT = ntop + 3
        L = self.L
        alive = self.geno >= 0
        cnt = np.bincount(self.geno[alive], minlength=NG)
        self.tag_geno = np.full((FL, NT), -1, np.int64)
        self.tag_lut = np.full((FL, NG), 2, np.int64)
        gs = np.nonzero(cnt)[0]
        if only is not None:
            gs = np.array([g for g in gs if int(g) in set(only)], dtype=np.int64)
        for j in range(FL):
            prod = gs[(gs >> j) & 1 == 1]
            prod = prod[np.argsort(-cnt[prod])][:ntop]
            self.tag_geno[j, 3:3 + len(prod)] = prod
            self.tag_lut[j, prod] = 3 + np.arange(len(prod))
        self.T = np.zeros((FL, NT, L, L), np.float32)
        self.T[:, 0] = self.M
        self.T2 = self.T.reshape(FL, NT, -1)
        self.NT = NT
        # flux[j, tag, recipient genotype*2 + habitat]
        self.flux = np.zeros((FL, NT, NG * 2))
        self.tags = True
        self.accum = False

    def stop_tags(self):
        self.tags = None
        self.T = self.T2 = None

    # ---------------- one step ----------------
    def step(self):
        p, rng = self.p, self.rng
        iota = p["iota"] * min(1.0, self.t / p["ramp"]) if p["ramp"] > 0 else p["iota"]
        alive = self.geno >= 0
        idx = np.nonzero(alive)[0]
        n = len(idx)
        g = self.geno[idx].astype(np.int64)
        B = ((g[:, None] >> KAR) & 1).astype(np.float32)
        ncar = B.sum(1)
        isE = (g & EUK) > 0
        Yv = np.where(isE, p["Y_e"], p["Y"]).astype(np.float32)[:, None]
        eff = np.ones((n, FL), np.float32)
        o2 = self.O22[idx]
        eff[:, ANX] = (1 - iota * B[:, OX]) / (1 + iota * o2 / p["K_O"])
        # light attenuation inside particles (revision after Carbon's review): OX (oxygenic
        # photosynthesis analogue) runs at light_p there, for both product and O2
        lightf = np.where(self.part[idx], p["light_p"], 1.0).astype(np.float32)
        eff[:, OX] *= lightf
        prod = B[:, :FL] * eff * Yv
        internal = np.minimum(prod * (1 - p["leak"]), p["theta"])
        deficit = p["theta"] - internal
        isE_ = (g & EUK) > 0
        if isE_.any():   # eukaryote takes up only what it requires (B12); revision after review
            deficit[isE_] = np.where(np.arange(FL) == B12, deficit[isE_], 0.0)
        here = self.M2[:, idx].T
        take = np.minimum(here, deficit)
        supply = (internal + take) / p["theta"]
        if isE.any():   # eukaryote needs only B12 (L4); other needs met internally / from medium
            supply[isE] = np.where(np.arange(FL) == B12, supply[isE], 1.0)
        # tag bookkeeping for uptake
        if self.tags:
            frac = np.where(here > 0, take / np.maximum(here, 1e-30), 0).astype(np.float32)
            Tsub = self.T2[:, :, idx]                      # (FL, NT, n)
            if self.accum:
                rk = g * 2 + self.part[idx]
                key = (np.arange(self.NT)[:, None] * (NG * 2) + rk[None, :]).ravel()
                for j in range(FL):
                    w = Tsub[j] * frac[:, j]
                    self.flux[j] += np.bincount(key, w.ravel(), minlength=self.NT * NG * 2
                                                ).reshape(self.NT, -1)
            self.T2[:, :, idx] = Tsub * (1 - frac.T)[:, None, :]
        self.M2[:, idx] -= take.T
        # private supply: carrier -> 1, else abiotic patch
        priv = np.maximum(B[:, PRIV], self.priv_map[idx])
        priv[isE] = 1.0
        smin = np.minimum(supply.min(1), priv.min(1))
        if isE.any():   # self-shading: eukaryote growth falls with eukaryote neighbours
            gn = self.geno[self.nbr[idx[isE]]]
            ne = ((gn >= 0) & ((gn & EUK) > 0)).sum(1)
            shade = np.ones(n, np.float32)
            shade[isE] = np.maximum(0.0, 1 - ne / p["shade_n"])
        else:
            shade = 1.0
        # planted positive control (off by default): a genome carrying all five leaky compatible
        # functions pays for them as one (economy of scope) -> selects stacked helpers
        n_eff = ncar - p["stack_discount"] * 4 * (B[:, :5].sum(1) == 5)
        grow = np.where(isE, p["rmax_e"] * smin * shade, p["rmax"] * (1 - p["cost"]) ** n_eff * smin)
        gh_ = self.ghost[idx] if self.ghost is not None else None
        if gh_ is not None:
            grow = np.where(gh_, 0.0, grow)
        if self.track_growth:
            gk = (g * 2 + self.part[idx])
            live = ~gh_ if gh_ is not None else slice(None)
            self.gsum += np.bincount(gk[live], grow[live], minlength=NG * 2)
            self.nsum += np.bincount(gk[live], minlength=NG * 2)
        # export + ambient
        exp_ = (prod * p["leak"]).T
        if gh_ is not None:
            exp_ = exp_ * (~gh_)[None, :]
        self.M2[:, idx] += exp_
        self.M2 += p["amb"]
        if self.tags:
            tg = self.tag_lut[:, g]                         # (FL, n)
            for j in range(FL):
                self.T2[j] += np.bincount(tg[j] * self.N + idx, exp_[j],
                                          minlength=self.NT * self.N).reshape(self.NT, -1)
            self.T2[:, 1] += p["amb"]
        # O2
        act = np.ones(n, np.float32) if gh_ is None else (~gh_).astype(np.float32)
        # ghosts are redox-inert: no O2 production, no respiration (revision after review)
        self.O22[idx] += p["o_y"] * B[:, OX] * np.where(isE, 2.0, 1.0) * act * lightf
        # particulate-C respiration: cells on particles respire resp_p_mult times faster
        rmult = np.where(self.part[idx], p["resp_p_mult"], 1.0)
        self.O22[idx] -= np.minimum(self.O22[idx], p["o_resp"] * act * rmult)
        w = ~self.part
        self.O22[w] += p["o_relax"] * (p["o_bg"] - self.O22[w])
        self.O2 *= (1 - p["o_decay"])
        diffuse(self.O2, self.oops)
        # decay + diffusion of products
        self.M *= (1 - p["decay"])
        diffuse(self.M, self.dops)
        if self.tags:
            self.T *= (1 - p["decay"])
            diffuse(self.T, self.dops)
        # births
        div = rng.random(n) < grow
        par = idx[div]
        tgt = self.nbr[par, rng.integers(0, 8, len(par))]
        ok = self.geno[tgt] < 0
        par, tgt = par[ok], tgt[ok]
        perm = rng.permutation(len(tgt))
        par, tgt = par[perm], tgt[perm]
        _, first = np.unique(tgt, return_index=True)
        par, tgt = par[first], tgt[first]
        child = self.geno[par].astype(np.int64)
        if self.evolve and len(child):
            loss = (rng.random((len(child), K)) < p["u_loss"]) @ POW
            child = np.where(child & EUK, child, child & ~loss)
        # deaths (among pre-existing cells)
        dmask = rng.random(n) < np.where(isE, p["mu_e"], p["mu"])
        if gh_ is not None:
            dmask &= ~gh_
        dead = idx[dmask]
        self.geno[dead] = -1
        self.geno[tgt] = child
        # HGT
        if self.evolve and p["h_hgt"] > 0:
            alive = np.nonzero(self.geno >= 0)[0]
            hc = alive[rng.random(len(alive)) < p["h_hgt"]]
            if len(hc):
                donor = self.geno[self.nbr[hc, rng.integers(0, 8, len(hc))]].astype(np.int64)
                gh = self.geno[hc].astype(np.int64)
                diff = np.where((donor >= 0) & ((donor & EUK) == 0) & ((gh & EUK) == 0), donor & ~gh, 0)
                diff &= FULL
                bits = ((diff[:, None] >> KAR) & 1) * rng.random((len(hc), K))
                pick = bits.argmax(1)
                gain = np.where(diff > 0, 1 << pick, 0)
                self.geno[hc] = (gh | gain).astype(np.int32)
        # water mixing
        a, b = self.pairs[rng.integers(4)]
        sw = rng.random(len(a)) < p["p_mix"]
        a, b = a[sw], b[sw]
        self.geno[a], self.geno[b] = self.geno[b].copy(), self.geno[a].copy()
        if self.ghost is not None:
            self.ghost[a], self.ghost[b] = self.ghost[b].copy(), self.ghost[a].copy()
        self.t += 1

    def snapshot_state(self):
        return self.geno.copy(), self.M.copy(), self.O2.copy()

    def restore(self, st):
        self.geno[:] = st[0]
        self.M[:] = st[1]
        self.O2[:] = st[2]
