"""P6 field counting: one assignment rule for DNA (capacity) and RNA (expression).

Protocol: P6 (frozen) + P6-A1 (approved) + run spec; table notes §5.2.
Input: SAM text, grouped so that all records of one read name are contiguous. Production input is
strobealign output streamed in read order (grouped by construction). Grouping is checked over a
bounded window of the last 200,000 names (constant memory); a violation raises ValueError.
CDS table: TSV with columns cds_id, contig, start (1-based), end (inclusive).

Rules (identical for DNA and RNA):
- Drop records with flags 0x4 (unmapped), 0x100 (secondary), 0x800 (supplementary).
- Optional MAPQ filter (sensitivity: --min-mapq 20).
- A fragment is one SE read, or one PE pair (both mates' primary records combined).
- Overlap = aligned reference bases (CIGAR M/=/X; D, N, I, S, H not counted) intersecting a CDS.
  Overlapping-mate convention: for PE, the aligned blocks of both mates are merged into their
  union first, so a reference base covered by both mates counts once per fragment.
- Assign the fragment to the CDS with strictly the largest overlap; exact ties at the maximum
  -> unassigned ('tie'); no overlap -> unassigned ('none').
Outputs per CDS: fragments (count) and overlap_bp (sum of the assigned fragments' overlap with it).
DNA capacity uses overlap_bp (depth x length over uniquely assigned fragments); RNA uses fragments.
"""
import bisect
import collections
import csv
import re
import sys

CIG = re.compile(r'(\d+)([MIDNSHP=X])')
DROP = 0x4 | 0x100 | 0x800


def load_cds(path):
    by = collections.defaultdict(list)
    with open(path) as fh:
        for r in csv.DictReader(fh, delimiter='\t'):
            by[r['contig']].append((int(r['start']), int(r['end']), r['cds_id']))
    idx = {}
    for c, lst in by.items():
        lst.sort()
        idx[c] = ([s for s, _, _ in lst], lst, max(e - s for s, e, _ in lst) + 1)
    return idx


def aligned_blocks(pos, cigar):
    """1-based inclusive reference blocks covered by aligned bases (M/=/X)."""
    out, ref = [], pos
    for n, op in CIG.findall(cigar):
        n = int(n)
        if op in 'M=X':
            out.append((ref, ref + n - 1)); ref += n
        elif op in 'DN':
            ref += n
    return out


def overlaps(idx, contig, blocks, acc):
    if contig not in idx:
        return
    starts, lst, maxlen = idx[contig]
    for bs, be in blocks:
        i = bisect.bisect_left(starts, bs - maxlen)
        while i < len(lst) and lst[i][0] <= be:
            s, e, cid = lst[i]
            ov = min(e, be) - max(s, bs) + 1
            if ov > 0:
                acc[cid] += ov
            i += 1


def merge(bl):
    out = []
    for s, e in sorted(bl):
        if out and s <= out[-1][1] + 1:
            out[-1] = (out[-1][0], max(out[-1][1], e))
        else:
            out.append((s, e))
    return out


def assign(acc):
    if not acc:
        return None, 'none'
    best = max(acc.values())
    top = [c for c, v in acc.items() if v == best]
    return (top[0], 'ok') if len(top) == 1 else (None, 'tie')


def count(sam_lines, idx, min_mapq=0):
    frag = collections.Counter(); bp = collections.Counter(); fate = collections.Counter()
    cur, blocks, seen_primary = None, collections.defaultdict(list), False

    def flush():
        if cur is None:
            return
        if not seen_primary:
            fate['dropped'] += 1; return
        acc = collections.Counter()
        for contig, bl in blocks.items():
            overlaps(idx, contig, merge(bl), acc)
        cid, why = assign(acc)
        fate[why] += 1
        if cid:
            frag[cid] += 1; bp[cid] += acc[cid]

    recent = collections.deque(maxlen=200000); recent_set = set()   # bounded grouping check
    for line in sam_lines:
        if line.startswith('@'):
            continue
        f = line.rstrip('\n').split('\t')
        qn, flag, rname, pos, mapq, cigar = f[0], int(f[1]), f[2], int(f[3]), int(f[4]), f[5]
        if qn != cur:
            flush()
            if qn in recent_set:
                raise ValueError(f'read {qn} not contiguous: input must be name-grouped')
            if cur is not None:
                if len(recent) == recent.maxlen:
                    recent_set.discard(recent[0])
                recent.append(cur); recent_set.add(cur)
            cur, blocks, seen_primary = qn, collections.defaultdict(list), False
        if flag & DROP or mapq < min_mapq or cigar == '*':
            continue
        seen_primary = True
        blocks[rname].extend(aligned_blocks(pos, cigar))
    flush()
    return frag, bp, fate


def main():
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument('--cds', required=True)
    ap.add_argument('--min-mapq', type=int, default=0)
    ap.add_argument('--out', required=True)
    a = ap.parse_args()
    idx = load_cds(a.cds)
    frag, bp, fate = count(sys.stdin, idx, a.min_mapq)
    with open(a.out, 'w') as fh:
        fh.write('cds_id\tfragments\toverlap_bp\n')
        for c in sorted(frag):
            fh.write(f'{c}\t{frag[c]}\t{bp[c]}\n')
    with open(a.out + '.fate', 'w') as fh:
        for k, v in sorted(fate.items()):
            fh.write(f'{k}\t{v}\n')


if __name__ == '__main__':
    main()
