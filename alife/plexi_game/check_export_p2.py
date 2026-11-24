"""Check Plexi Pond's TSV export against the team's P2 parser (scratch/p2_scorer.py).

Usage: python3 -I check_export_p2.py <export.tsv> [path/to/p2_scorer.py]
Parses each MAG's per_mag section with parse_per_mag_rows(complete_cds_table=True) and compares
the parser's edge_status for every pathway CDS with the simulated boundary truth.
"""
import importlib.util, sys
from collections import defaultdict

path = sys.argv[1]
scorer = sys.argv[2] if len(sys.argv) > 2 else '/data/dev/plexityper/scratch/p2_scorer.py'
import os; sys.path.insert(0, os.path.dirname(os.path.abspath(scorer)))
spec = importlib.util.spec_from_file_location('p2_scorer', scorer)
p2 = importlib.util.module_from_spec(spec); sys.modules['p2_scorer'] = p2; spec.loader.exec_module(p2)

sections, cur = defaultdict(list), None
for line in open(path):
    line = line.rstrip('\n')
    if line.startswith('## '):
        cur = line[3:].split(' ')[0] if not line.startswith('## KO') else cur
        continue
    if line and cur: sections[cur].append(line.split('\t'))

def table(name):
    rows = sections[name]; hdr = rows[0]
    return [dict(zip(hdr, r)) for r in rows[1:]]

per = table('per_mag'); ctg = table('contigs'); truth = table('boundary'); comp = table('completeness')
by_mag = defaultdict(list); [by_mag[r['mag_id']].append(r) for r in per]
lens = defaultdict(list); [lens[r['mag_id']].append(int(r['contig_length'])) for r in ctg]
dom = {r['mag_id']: r['domain'] for r in comp}
tmap = {(r['mag_id'], r['protein_id']): r['boundary_true'] for r in truth}
agree = n = 0; status_counts = defaultdict(int)
for mag, rows in by_mag.items():
    _, hits = p2.parse_per_mag_rows(rows, mag, domain=dom[mag], contig_lengths=lens[mag], complete_cds_table=True)
    for h in hits:
        if not h.kos: continue
        t = tmap[(mag, h.protein_id)]; n += 1; status_counts[(t, h.edge_status)] += 1
        agree += (t == h.edge_status)
print('MAGs', len(by_mag), '| eukaryotic', sum(1 for m in dom.values() if m == 'Eukaryota'),
      '| pathway CDS', n, '| parser edge_status == simulated boundary:', agree, '/', n)
for k, v in sorted(status_counts.items()): print('  truth %-8s parser %-8s %d' % (k[0], k[1], v))
sys.exit(0 if n and agree == n else 1)
