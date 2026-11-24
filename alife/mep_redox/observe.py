"""Synthetic observations in the shape of the team's real data: depth-resolved chemistry,
MAG abundance (DNA), per-MAG CDS/KO records with contig and CDS IDs, feature order, blank-KO
CDSs and completeness, and per-CDS RNA counts.  DNA, RNA and chemistry share event_id and
sample_id.  Observers get MAG ids only; the MAG -> metabolic type key is ground truth."""
import os
import numpy as np
from model import TYPES, STEP_KOS, B12_PRODUCERS, SPECIES, S, RXN, RI, T

HOUSEKEEPING = ["K02886", "K02906", "K02874", "K02946", "K02961", "K02967", "K03043",
                "K02470", "K03553", "K04043", "K01872", "K02355", "K03046", "K01889"]
SAMPLE_DEPTHS = [1, 3, 5, 7, 9, 11, 14, 17, 20, 24, 28]
# reaction -> specialist type used when emitting observations from enzyme fields
REACTION_CARRIER = {"AER": "AER", "DMP": "ROS", "NAR": "NAR", "NIR": "NIR", "SRB": "SRB",
                    "FER": "FER", "MET": "MET", "HSR": "SRBH", "AMO": "AOB", "NXR": "NOB",
                    "AMX": "AMX", "SOX": "SOB", "SNR": "SNR", "MOX": "MOB"}


def build_catalogue(seed=12345):
    """One genome per metabolic type: contigs, CDS order, KO calls (blank for many)."""
    rng = np.random.default_rng(seed)
    names = list(TYPES)
    ids = [f"MAG{k:03d}" for k in rng.permutation(len(names)) + 1]
    cat = {}
    for t, mid in zip(names, ids):
        ops = []
        steps = TYPES[t] if t != "ALG" else ["PHO"]
        for st in steps:
            ops.append((st, list(STEP_KOS[st])))
        if t in B12_PRODUCERS:
            ops.append(("B12syn", list(STEP_KOS["B12syn"])))
        if t == "ALG":
            ops.append(("B12dep", ["K00548"]))          # metH, no metE: B12 auxotroph
        else:
            ops.append(("metE" if rng.random() < 0.6 else "metH",
                        ["K00549"] if rng.random() < 0.6 else ["K00548"]))
        n_hk = int(rng.integers(50, 90))
        hk = [(HOUSEKEEPING[rng.integers(len(HOUSEKEEPING))] if rng.random() < 0.55 else "")
              for _ in range(n_hk)]
        n_ctg = int(rng.integers(4, 10))
        contigs = [[] for _ in range(n_ctg)]
        for g in hk:
            contigs[rng.integers(n_ctg)].append(("hk", g))
        for lab, kos in ops:                                  # operons kept contiguous
            c = contigs[rng.integers(n_ctg)]
            pos = int(rng.integers(len(c) + 1))
            c[pos:pos] = [(lab, k) for k in kos]
        recs = []
        for ci, c in enumerate(contigs):
            x = 1
            for fi, (lab, ko) in enumerate(c):
                ln = int(rng.integers(600, 2400))
                recs.append(dict(contig_id=f"{mid}_c{ci + 1:02d}", cds_id=f"{mid}_c{ci + 1:02d}_{fi + 1:04d}",
                                 feature_index=fi + 1, start=x, end=x + ln - 1,
                                 strand="+" if rng.random() < 0.5 else "-", ko=ko, role=lab))
                x += ln + int(rng.integers(20, 300))
        cat[t] = dict(mag_id=mid, cds=recs, length=sum(r["end"] - r["start"] + 1 for r in recs))
    return cat


def _write(path, header, rows):
    with open(path, "w") as fh:
        fh.write("\t".join(header) + "\n")
        for r in rows:
            fh.write("\t".join(str(r[h]) for h in header) + "\n")


def emit(outdir, truthdir, run_id, casts, cat, seed, dz=1.0):
    """casts: list of dict(t, C (NS,NZ), B {type: profile}, flux {(type, rxn): profile},
    G profile).  Writes observer files to outdir and the key to truthdir."""
    rng = np.random.default_rng(seed)
    os.makedirs(outdir, exist_ok=True)
    os.makedirs(truthdir, exist_ok=True)
    zi = [int(d / dz - 0.5 + 0.5) for d in SAMPLE_DEPTHS]
    # assembly: per-run contig loss -> completeness
    kept, comp_rows, cds_rows = {}, [], []
    for t, g in cat.items():
        target = rng.uniform(0.6, 1.0)
        ctgs = sorted({r["contig_id"] for r in g["cds"]})
        order = list(rng.permutation(ctgs))
        keep = set(ctgs)
        n_tot = len(g["cds"])
        for c in order:
            n_keep = sum(1 for r in g["cds"] if r["contig_id"] in keep and r["contig_id"] != c)
            if n_keep / n_tot >= target and len(keep) > 1:
                keep.discard(c)
        rec = [r for r in g["cds"] if r["contig_id"] in keep]
        kept[t] = rec
        comp = len(rec) / n_tot
        clen = {}
        for r in rec:
            clen[r["contig_id"]] = max(clen.get(r["contig_id"], 0), r["end"] + 150)
        comp_rows.append(dict(mag_id=g["mag_id"], domain="Eukaryota" if t == "ALG" else "Bacteria",
                              completeness=round(min(1, comp + rng.normal(0, 0.03)) * 100, 1),
                              contamination=round(abs(rng.normal(0, 1.0)), 1), n_contigs=len(keep),
                              n_cds=len(rec), genome_size_bp=sum(clen.values()),
                              contig_ids=",".join(sorted(clen)),
                              contig_lengths=",".join(str(clen[c]) for c in sorted(clen))))
        for r in rec:
            cds_rows.append(dict(mag_id=g["mag_id"], contig_id=r["contig_id"], protein_id=r["cds_id"],
                                 feature_index=r["feature_index"], start=r["start"], end=r["end"],
                                 strand=r["strand"], KO=r["ko"]))
    _write(os.path.join(outdir, "mag_completeness.tsv"),
           ["mag_id", "domain", "completeness", "contamination", "n_contigs", "n_cds",
            "genome_size_bp", "contig_ids", "contig_lengths"], comp_rows)
    _write(os.path.join(outdir, "per_mag_cds.tsv"),
           ["mag_id", "contig_id", "protein_id", "feature_index", "start", "end", "strand", "KO"],
           cds_rows)
    _write(os.path.join(truthdir, "mag_key.tsv"), ["mag_id", "metabolic_type", "steps"],
           [dict(mag_id=g["mag_id"], metabolic_type=t, steps="+".join(TYPES[t]) or "PHO")
            for t, g in cat.items()])
    ev_rows, chem_rows, ab_rows, rna_rows = [], [], [], []
    for k, cs in enumerate(casts):
        ev = f"{run_id}_e{k + 1}"
        ev_rows.append(dict(event_id=ev, run_id=run_id, day=round(cs["t"], 2), assays="chem;dna;rna"))
        for d, z in zip(SAMPLE_DEPTHS, zi):
            sid = f"{ev}_z{d:02d}"
            row = dict(event_id=ev, sample_id=sid, depth_m=d, T_K=T, pH=7.5)
            for sp in SPECIES:
                if sp == "POC":
                    continue
                v = cs["C"][S[sp], z] * rng.lognormal(0, 0.05)
                dl = 5e-4 if sp == "B12" else 0.01
                row[sp] = f"{v:.4g}" if v >= dl else "0"
            chem_rows.append(row)
            # DNA
            names = list(cat)
            bio = np.array([max(cs["B"].get(t, np.zeros(1))[z] if t in cs["B"] else 0, 0) for t in names])
            L = np.array([cat[t]["length"] for t in names], float)
            w = bio * L
            unb = 0.25 * w.sum() + 1e-12
            pr = np.append(w, unb) / (w.sum() + unb)
            cnt = rng.multinomial(1_000_000, pr)
            for t, c in zip(names, cnt[:-1]):
                ab_rows.append(dict(sample_id=sid, event_id=ev, mag_id=cat[t]["mag_id"], reads=int(c)))
            ab_rows.append(dict(sample_id=sid, event_id=ev, mag_id="unbinned", reads=int(cnt[-1])))
            # RNA: pathway CDS ~ step flux, housekeeping ~ biomass, algal PHO ~ fixation
            keys, wts = [], []
            for t in names:
                b = bio[names.index(t)]
                if b <= 0:
                    continue
                for r in kept[t]:
                    role = r["role"]
                    if role in RI:
                        fl = cs["flux"].get((t, role))
                        x = 0.0 if fl is None else max(fl[z], 0)
                        x *= 20.0
                    elif role == "PHO":
                        x = 5.0 * max(cs["G"][z], 0)
                    elif role == "B12syn":
                        x = 0.3 * b
                    else:
                        x = 0.05 * b
                    if x > 0:
                        keys.append((t, r))
                        wts.append(x * rng.lognormal(0, 0.3))
            if wts:
                wts = np.array(wts)
                cnt = rng.multinomial(1_000_000, wts / wts.sum())
                for (t, r), c in zip(keys, cnt):
                    if c > 0:
                        rna_rows.append(dict(sample_id=sid, event_id=ev, mag_id=cat[t]["mag_id"],
                                             protein_id=r["cds_id"], KO=r["ko"], reads=int(c)))
    _write(os.path.join(outdir, "events.tsv"), ["event_id", "run_id", "day", "assays"], ev_rows)
    _write(os.path.join(outdir, "chemistry.tsv"),
           ["event_id", "sample_id", "depth_m", "T_K", "pH"] + [s for s in SPECIES if s != "POC"],
           chem_rows)
    _write(os.path.join(outdir, "abundance.tsv"), ["sample_id", "event_id", "mag_id", "reads"], ab_rows)
    _write(os.path.join(outdir, "rna_counts.tsv"),
           ["sample_id", "event_id", "mag_id", "protein_id", "KO", "reads"], rna_rows)


def cast_from_population(m, y, t):
    P = m.tendencies(t, y, parts=True)
    C, B = m.unpack(y)
    Jt = P["J"] + P["Jint"]
    flux = {}
    for j, tn in enumerate(m.types):
        for r in RXN:
            if m.M[j, RI[r]]:
                flux[(tn, r)] = Jt[j, RI[r]]
    return dict(t=t, C=C.copy(), B={tn: B[j].copy() for j, tn in enumerate(m.types)},
                flux=flux, G=P["G"])


def cast_from_fixed(f, y, t):
    """Observations from an enzyme field: each reaction's enzyme is assigned to its
    specialist carrier type (REACTION_CARRIER)."""
    P = f.tendencies(t, y, parts=True)
    C, B = f.unpack(y)
    Jt = (P["J"] + P["Jint"]).sum(0)
    Bt, flux = {}, {}
    if f.types[0] == "OMNI":
        E = None
        # greedy: enzyme follows the instantaneous allocation
        e = f.allocation(P["rate"] * np.maximum(-P["dG"], 0))[0] * B[0][None]
    else:
        e = B[:len(RXN)]
    for r in RXN:
        tn = REACTION_CARRIER[r]
        Bt[tn] = Bt.get(tn, 0) + e[RI[r]]
        flux[(tn, r)] = Jt[RI[r]]
    Bt["ALG"] = B[f.ialg]
    return dict(t=t, C=C.copy(), B=Bt, flux=flux, G=P["G"])
