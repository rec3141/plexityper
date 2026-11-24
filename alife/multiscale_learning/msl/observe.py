"""Synthetic observations: sampling designs, sequencing noise, compositional closure, and
catalogue (per-MAG CDS records with contig/CDS IDs, order, blank KOs, completeness).

All observation tables are keyed by a shared event_id (station x arm x sampling step), so the
DNA, RNA, SNV and amplicon assays of one collection carry the same key and cannot be split across
train/test folds.
"""
import numpy as np
from .model import M, R, L, DT

DAY_STEPS = int(round(1 / DT))   # 10
HOURS = np.arange(DAY_STEPS) * 24 * DT   # 0, 2.4, ..., 21.6


def step_offset(day, hour_slot):
    return int(day) * DAY_STEPS + int(hour_slot)


def make_designs(obs_years, seed=11):
    """Returns {design: list of (day, hour_slot)}. hour_slot 5 = 12:00."""
    rng = np.random.default_rng(seed)
    D = {}
    ndays = obs_years * 365
    D['daily_ideal'] = [(d, 5) for d in range(ndays)]
    D['weekly_2y_noon'] = [(d, 5) for d in range(0, 2 * 365, 7)]
    D['weekly_2y'] = [(d, int(rng.choice([3, 4, 5, 6]))) for d in range(0, 2 * 365, 7)]
    # weekly with the clock time drifting 4.8 h later each visit: aliases the diel cycle
    D['weekly_2y_drift'] = [(d, (k * 2) % DAY_STEPS) for k, d in enumerate(range(0, 2 * 365, 7))]
    # Lake-Winnipeg-like: monthly open-water (May-Oct) for 4 years (years 2-5), 15% missed
    mw = []
    for y in range(1, 5):
        for mo_start in (121, 152, 182, 213, 244, 274):
            if rng.random() < 0.15:
                continue
            mw.append((y * 365 + mo_start + int(rng.integers(0, 28)), int(rng.choice([4, 5, 6]))))
    D['monthly_4y_gaps'] = mw
    # CMO-like: two winter campaigns (years 3 and 4, days 20-95), every 3-4 days
    camp = []
    for y in (2, 3):
        d = y * 365 + 20
        while d < y * 365 + 95:
            camp.append((d, 5)); d += int(rng.choice([3, 4]))
    D['two_campaigns'] = camp
    # sub-daily burst: every 4.8 h for 14 days in summer of year 2
    D['subdaily_burst'] = [(365 + 190 + k // 5, (k % 5) * 2) for k in range(14 * 5)]
    return {k: [(d, h) for d, h in v if d < ndays] for k, v in D.items()}


def union_offsets(designs):
    return np.unique([step_offset(d, h) for v in designs.values() for d, h in v])


# ------------------------------------------------------------------ catalogue
GENE_NAMES = (['hk%02d' % k for k in range(10)] + ['bio%d' % m for m in range(M)] +
              ['upt%d' % m for m in range(M)] + ['cat%d' % r for r in range(R)])


def ko_of(gene):
    """Synthetic KO identifiers in the unused K90001+ range, so p2_scorer.parse_ko_field (which
    keeps only K + 5 digits) accepts them; they are NOT KEGG orthologs (see synthetic_ko_map.tsv)."""
    return 'K%05d' % (90001 + GENE_NAMES.index(gene))


def build_catalogue(rec, pool, rng, rel_thresh=5e-4, recovery=0.8, n_hyp=40):
    """MAG catalogue: lineages whose max relative abundance over the recorded window passes a
    threshold, each recovered with probability `recovery`. Gene content = consensus (F > 0.5,
    averaged over the window). Lost biosynthesis genes leave a blank-KO remnant CDS 30% of the time.
    """
    B = rec['B'].astype(float)
    rel = B / B.sum(1, keepdims=True)
    elig = np.nonzero(rel.max(0) > rel_thresh)[0]
    keep = elig[rng.random(len(elig)) < recovery]
    Fm = rec['F'].mean(0)
    mags, cds = [], []
    for n, i in enumerate(keep):
        mag = 'MAG%04d' % (n + 1)
        genes = []
        for gname in GENE_NAMES:
            if gname.startswith('bio'):
                m = int(gname[3:])
                if Fm[i, m] > 0.5:
                    genes.append((gname, ko_of(gname)))
                elif rng.random() < 0.3:
                    genes.append((gname + '_remnant', ''))
            else:
                genes.append((gname, ko_of(gname)))
        genes += [('hyp%02d' % k, '') for k in range(n_hyp)]
        order = rng.permutation(len(genes))
        genes = [genes[o] for o in order]
        n_contig = int(rng.integers(8, 21))
        cuts = np.sort(rng.choice(np.arange(1, len(genes)), n_contig - 1, replace=False))
        contig_of = np.searchsorted(cuts, np.arange(len(genes)), side='right')
        target = rng.uniform(0.6, 1.0)
        contigs = rng.permutation(n_contig)
        dropped = set()
        for c in contigs:     # drop contigs until completeness reaches target
            kept = np.isin(contig_of, list(set(range(n_contig)) - dropped - {c})).mean()
            if kept >= target:
                dropped.add(c)
        present_gene = {}
        pos = {}
        clen = {}
        for gidx, (gname, ko) in enumerate(genes):
            c = contig_of[gidx]
            if c in dropped:
                continue
            k = pos.get(c, 0); pos[c] = k + 1
            start = 100 + k * 1100 + int(rng.integers(0, 200))
            end = start + int(rng.integers(600, 1000))
            cid = '%s_c%02d' % (mag, c)
            cds.append(dict(mag_id=mag, contig_id=cid, protein_id='%s_%d' % (cid, k + 1),
                            feature_index=k + 1, start=start, end=end,
                            strand=int(rng.choice([1, -1])), KO=ko, EC=''))
            clen[cid] = end + 200
            present_gene[gname] = True
        hk_present = np.mean([('hk%02d' % k) in present_gene for k in range(10)])
        comp_true = 1 - len(dropped) / n_contig if n_contig else 1.0
        kept_frac = np.mean([contig_of[g] not in dropped for g in range(len(genes))])
        mags.append(dict(mag_id=mag, local=int(i), lineage=int(rec['idx'][i]),
                         completeness_true=round(float(kept_frac), 3),
                         completeness_est=round(float(np.clip(hk_present + rng.normal(0, 0.05), 0, 1)), 3),
                         contamination_est=round(float(abs(rng.normal(0, 0.01))), 3),
                         n_contigs=n_contig - len(dropped), present=present_gene, domain='Bacteria',
                         contig_lengths=[clen[c_] for c_ in sorted(clen)]))
    # one eukaryotic (algal) bin: Eukaryota domain, B12-uptake but no B12 synthesis, releases 1-3
    egenes = [('hk%02d' % k, ko_of('hk%02d' % k)) for k in range(10)] + \
             [('bio%d' % m, ko_of('bio%d' % m)) for m in (1, 2, 3)] + [('upt0', ko_of('upt0'))] + \
             [('hyp%03d' % k, '') for k in range(200)]
    for k, (gname, ko) in enumerate(egenes):
        c = k // 20; cid = 'EUK0001_c%02d' % c; kk = k % 20
        cds.append(dict(mag_id='EUK0001', contig_id=cid, protein_id='%s_%d' % (cid, kk + 1),
                        feature_index=kk + 1, start=100 + kk * 1500, end=100 + kk * 1500 + 1200,
                        strand=1, KO=ko, EC=''))
    return mags, cds


# ------------------------------------------------------------------ per-event observation
def observe_events(rec, rows, pool, cfg, mags, rng, dna_depth=2e6, rna_depth=2e6, amp_depth=4e4,
                   dm_alpha_dna=400., dm_alpha_rna=150., dm_alpha_amp=300.):
    """Observation tables for record rows `rows`. Returns dict of arrays (events x features)."""
    idx = rec['idx']; gs = pool['gsize'][idx]
    loc = np.array([m['local'] for m in mags], int)
    nE = len(rows)
    A = rec['A'][rows].astype(float) if cfg['euk'] == 'dynamic' else np.full(nE, cfg['A_flat'] if cfg['euk'] == 'flat' else 0.0)
    out = dict(rows=np.asarray(rows))
    # ---- amplicon (all lineages + chloroplast reads from algae): ASV-like table
    B = rec['B'][rows].astype(float)
    w = np.hstack([B, (0.5 * A)[:, None]])
    out['amp'] = _dm_counts(w, amp_depth, dm_alpha_amp, rng)
    # ---- metagenome reads: per lineage (gsize-weighted) + eukaryote reads (large genomes)
    wd = np.hstack([B * gs[None, :], (8.0 * A)[:, None]])
    dna = _dm_counts(wd, dna_depth * np.exp(rng.normal(0, 0.3, nE)), dm_alpha_dna, rng)
    out['dna_lineage_reads'] = dna
    mag_reads = dna[:, loc]
    unbinned = dna[:, :-1].sum(1) - mag_reads.sum(1)
    out['mag_reads'] = mag_reads; out['unbinned_reads'] = unbinned; out['euk_reads'] = dna[:, -1]
    cov = mag_reads * 150.0 / (gs[loc][None, :] * 1e6)                    # mean depth (x)
    out['mag_cov'] = cov
    # ---- gene-relative coverage (genome layer observable: carrier fraction F)
    F = rec['F'][rows][:, loc, :].astype(float)
    lam = cov[:, :, None] * 1.0
    gene_reads = rng.poisson(lam * F)
    with np.errstate(invalid='ignore', divide='ignore'):
        out['gene_cov_ratio'] = np.where(lam > 0.5, gene_reads / lam, np.nan)   # (E, MAG, M)
    # ---- RNA: transcripts per (lineage, gene); eukaryote transcripts dominate in blooms
    Fa = rec['F'][rows].astype(float); Es = rec['Esyn'][rows].astype(float)
    Eu = rec['Eup'][rows].astype(float); Ec = rec['Ecat'][rows].astype(float)
    lev = np.concatenate([np.ones(B.shape + (1,)), 2 * Fa * Es, 1.5 * Eu, 3 * Ec], axis=2)  # hk,bio,upt,cat
    lev = lev * np.exp(rng.normal(0, 0.3, lev.shape))
    wr = (B[:, :, None] * lev).reshape(nE, -1)
    # algal transcription is strongly diel (photosynthesis genes peak around midday); this
    # modulates the eukaryote share of the RNA library and, by closure, every prokaryotic share
    hour = (rec['t'][rows] % 1.0) * 24.0
    diel_euk = 1.0 + 0.8 * np.cos(2 * np.pi * (hour - 12.0) / 24.0)
    wr = np.hstack([wr, (20.0 * A * diel_euk)[:, None]])
    rna = _dm_counts(wr, rna_depth * np.exp(rng.normal(0, 0.3, nE)), dm_alpha_rna, rng)
    rna_lin = rna[:, :-1].reshape(nE, B.shape[1], 1 + 2 * M + R)
    out['rna_mag'] = rna_lin[:, loc, :]                                    # (E, MAG, genes)
    out['rna_unbinned'] = rna_lin.sum((1, 2)) - out['rna_mag'].sum((1, 2))
    out['rna_euk'] = rna[:, -1]
    # ---- SNVs on MAGs
    p = rec['Psnv'][rows][:, loc, :].astype(float)
    depth = rng.poisson(np.repeat(cov[:, :, None], L, axis=2))
    e = 0.005
    alt = rng.binomial(depth, p * (1 - e) + (1 - p) * e)
    with np.errstate(invalid='ignore', divide='ignore'):
        out['snv_freq'] = np.where(depth >= 5, alt / np.maximum(depth, 1), np.nan)
    out['snv_depth'] = depth
    # ---- measured environment (with sensor noise); storms and pools are NOT measured
    out['env_T'] = rec['T'][rows].astype(float) + rng.normal(0, 0.2, nE)
    out['env_S'] = rec['S'][rows].astype(float) * np.exp(rng.normal(0, 0.15, (nE, R)))
    out['chl'] = A * np.exp(rng.normal(0, 0.2, nE))                         # chlorophyll proxy
    return out


def _dm_counts(w, depth, alpha, rng):
    """Dirichlet-multinomial counts with per-row depth; rows of w are unnormalised weights."""
    w = np.maximum(w, 0) + 1e-15
    p = w / w.sum(1, keepdims=True)
    depth = np.broadcast_to(np.asarray(depth, float), (w.shape[0],))
    out = np.empty(w.shape, np.int64)
    for k in range(w.shape[0]):
        q = rng.gamma(alpha * p[k] + 1e-12)
        q = q / q.sum()
        out[k] = rng.multinomial(int(depth[k]), q)
    return out


def write_tables(prefix, design, events, obs, mags, cds, gene_names=None):
    """Write TSVs in the shapes of the team's real tables. events: list of dicts."""
    import os
    os.makedirs(prefix, exist_ok=True)
    with open(f'{prefix}/events.tsv', 'w') as fh:
        keys = list(events[0].keys())
        fh.write('\t'.join(keys) + '\n')
        for e in events:
            fh.write('\t'.join(str(e[k]) for k in keys) + '\n')
    ev = [e['event_id'] for e in events]
    magids = [m['mag_id'] for m in mags]
    # abundance (MAG read counts + unbinned + eukaryote)
    with open(f'{prefix}/dna_mag_counts.tsv', 'w') as fh:
        fh.write('mag_id\t' + '\t'.join('DNA_' + e for e in ev) + '\n')
        for j, mid in enumerate(magids):
            fh.write(mid + '\t' + '\t'.join(str(int(x)) for x in obs['mag_reads'][:, j]) + '\n')
        fh.write('unbinned\t' + '\t'.join(str(int(x)) for x in obs['unbinned_reads']) + '\n')
        fh.write('EUK0001\t' + '\t'.join(str(int(x)) for x in obs['euk_reads']) + '\n')
    with open(f'{prefix}/amplicon_counts.tsv', 'w') as fh:
        fh.write('asv_id\t' + '\t'.join('AMP_' + e for e in ev) + '\n')
        n = obs['amp'].shape[1] - 1
        for j in range(n):
            fh.write('ASV%04d\t' % (j + 1) + '\t'.join(str(int(x)) for x in obs['amp'][:, j]) + '\n')
        fh.write('chloroplast\t' + '\t'.join(str(int(x)) for x in obs['amp'][:, -1]) + '\n')
    gn = ['hk%02d' % k for k in range(1)] + ['bio%d' % m for m in range(M)] + \
         ['upt%d' % m for m in range(M)] + ['cat%d' % r for r in range(R)]
    with open(f'{prefix}/rna_counts_long.tsv', 'w') as fh:
        fh.write('event_id\tmag_id\tko\tcount\n')
        for k, e in enumerate(ev):
            for j, mid in enumerate(magids):
                for gi, g in enumerate(gn):
                    if g.startswith('bio') and not mags[j]['present'].get(g, False):
                        continue
                    c = int(obs['rna_mag'][k, j, gi])
                    if c:
                        fh.write(f'{e}\t{mid}\t{ko_of(g)}\t{c}\n')
            fh.write(f'{e}\tunbinned\t\t{int(obs["rna_unbinned"][k])}\n')
            fh.write(f'{e}\tEUK0001\t\t{int(obs["rna_euk"][k])}\n')
    with open(f'{prefix}/gene_cov_ratio_long.tsv', 'w') as fh:
        fh.write('event_id\tmag_id\tko\tcoverage_ratio\n')
        for k, e in enumerate(ev):
            for j, mid in enumerate(magids):
                for m in range(M):
                    v = obs['gene_cov_ratio'][k, j, m]
                    if np.isfinite(v):
                        fh.write(f'{e}\t{mid}\t{ko_of("bio%d" % m)}\t{v:.3f}\n')
    with open(f'{prefix}/snv_freq_long.tsv', 'w') as fh:
        fh.write('event_id\tmag_id\tsite\tdepth\talt_freq\n')
        for k, e in enumerate(ev):
            for j, mid in enumerate(magids):
                for s in range(L):
                    v = obs['snv_freq'][k, j, s]
                    if np.isfinite(v):
                        fh.write(f'{e}\t{mid}\t{mid}_snv{s}\t{int(obs["snv_depth"][k, j, s])}\t{v:.3f}\n')
    with open(f'{prefix}/env.tsv', 'w') as fh:
        fh.write('event_id\ttemp\tS0\tS1\tS2\tchl\n')
        for k, e in enumerate(ev):
            fh.write(f'{e}\t{obs["env_T"][k]:.2f}\t' + '\t'.join(f'{x:.3f}' for x in obs['env_S'][k])
                     + f'\t{obs["chl"][k]:.3f}\n')
    with open(f'{prefix}/genomes.tsv', 'w') as fh:
        fh.write('mag_id\tdomain\tcompleteness_est\tcontamination_est\tn_contigs\tgenome_size_bp\tcontig_lengths\n')
        for m in mags:
            cl = m['contig_lengths']
            fh.write(f"{m['mag_id']}\t{m['domain']}\t{m['completeness_est']}\t{m['contamination_est']}\t"
                     f"{m['n_contigs']}\t{sum(cl)}\t{','.join(str(x) for x in cl)}\n")
        ecl = {}
        for c in cds:
            if c['mag_id'] == 'EUK0001':
                ecl[c['contig_id']] = max(ecl.get(c['contig_id'], 0), c['end'] + 200)
        cl = [ecl[k] for k in sorted(ecl)]
        fh.write(f"EUK0001\tEukaryota\tnan\tnan\t{len(cl)}\t{sum(cl)}\t{','.join(str(x) for x in cl)}\n")
    with open(f'{prefix}/synthetic_ko_map.tsv', 'w') as fh:
        fh.write('KO\tsynthetic_role\n')
        for g in GENE_NAMES:
            fh.write(f'{ko_of(g)}\t{g}\n')
    with open(f'{prefix}/per_mag_cds.tsv', 'w') as fh:
        keys = ['mag_id', 'contig_id', 'protein_id', 'feature_index', 'start', 'end', 'strand', 'KO', 'EC']
        fh.write('\t'.join(keys) + '\n')
        for c in cds:
            fh.write('\t'.join(str(c[k]) for k in keys) + '\n')
    with open(f'{prefix}/per_mag.tsv', 'w') as fh:     # the legacy two-column view
        fh.write('mag_id\tKO\n')
        for c in cds:
            if c['KO']:
                fh.write(f"{c['mag_id']}\t{c['KO']}\n")
