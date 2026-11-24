"""Aggregate per-BAM covered intervals into MAG breadth (spec scratch/p2_breadth_job_spec.md).

breadth(MAG, unit) = |union of covered intervals over the unit's BAMs, on the MAG's contigs|
                     / full MAG length (all contigs, including zero-coverage ones).

Units: each BAM alone (always written), and events when an event map is given
(TSV with columns library, event; library = flowcell_barcode as in the BAM names; per P2 l.51
an event is campaign x station, with same-cast depths and replicates grouped).

  python3 -I scratch/p2_breadth_aggregate.py <dir with *.covered.bed.gz> [event_map.tsv]
Writes results/p2/breadth/breadth_per_bam.tsv and, with a map, breadth_per_event.tsv.
"""
import csv
import glob
import gzip
import os
import sys
from collections import defaultdict

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(ROOT, "results", "p2", "breadth")
C2M = os.path.join(OUT, "contig_to_mag.tsv")
TARGETS = os.path.join(OUT, "targets_v2.bed")
TARGETS_SHA12 = "f48d55128984"


def load_c2m():
    c2m, mag_len = {}, defaultdict(int)
    with open(C2M) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            c2m[r["contig"]] = r["mag"]
            mag_len[r["mag"]] += int(r["length"])
    return c2m, mag_len


def load_extra(c2m):
    """Non-representative target contigs (unbinned, eukaryotic, non-rep bins): contig -> length.
    Kept out of representative MAG lengths; reported per contig (Astra 249754118818021376)."""
    import hashlib
    with open(TARGETS, "rb") as fh:
        assert hashlib.sha256(fh.read()).hexdigest()[:12] == TARGETS_SHA12, "targets BED hash mismatch"
    extra = {}
    with open(TARGETS) as fh:
        for line in fh:
            c, s, e = line.split("\t")[:3]
            if c not in c2m:
                extra[c] = int(e) - int(s)
    return extra


def contig_breadth(ivs_list, extra):
    per_contig = defaultdict(list)
    for iv in ivs_list:
        for c, lst in iv.items():
            if c in extra:
                per_contig[c].extend(lst)
    return {c: union_len(l) / extra[c] for c, l in per_contig.items()}


def read_intervals(path):
    iv = defaultdict(list)
    with gzip.open(path, "rt") as fh:
        for line in fh:
            c, s, e = line.split("\t")[:3]
            iv[c].append((int(s), int(e)))
    return iv


def union_len(intervals):
    tot, cur_s, cur_e = 0, None, None
    for s, e in sorted(intervals):
        if cur_e is None or s > cur_e:
            if cur_e is not None:
                tot += cur_e - cur_s
            cur_s, cur_e = s, e
        else:
            cur_e = max(cur_e, e)
    if cur_e is not None:
        tot += cur_e - cur_s
    return tot


def breadth(ivs_list, c2m, mag_len):
    per_contig = defaultdict(list)
    for iv in ivs_list:
        for c, lst in iv.items():
            per_contig[c].extend(lst)
    cov = defaultdict(int)
    for c, lst in per_contig.items():
        m = c2m.get(c)
        if m is not None:
            cov[m] += union_len(lst)
    return {m: cov.get(m, 0) / L for m, L in mag_len.items()}


def main():
    d = sys.argv[1]
    emap = sys.argv[2] if len(sys.argv) > 2 else None
    c2m, mag_len = load_c2m()
    extra = load_extra(c2m)
    mags = sorted(mag_len)
    files = sorted(glob.glob(os.path.join(d, "*.covered.bed.gz")))
    libs = [os.path.basename(f).replace(".covered.bed.gz", "") for f in files]
    per_lib = {}
    with open(os.path.join(OUT, "breadth_per_bam.tsv"), "w") as fh:
        fh.write("mag\t" + "\t".join(libs) + "\n")
        cols, extra_cols = [], []
        for f, lib in zip(files, libs):
            iv = read_intervals(f)
            per_lib[lib] = iv if emap else None
            cols.append(breadth([iv], c2m, mag_len))
            extra_cols.append(contig_breadth([iv], extra))
        for m in mags:
            fh.write(m + "\t" + "\t".join(f"{c[m]:.5f}" for c in cols) + "\n")
    with open(os.path.join(OUT, "breadth_extra_contigs_per_bam.tsv"), "w") as fh:
        # long format, nonzero only: library, contig, breadth
        fh.write("library\tcontig\tbreadth\n")
        for lib, ec in zip(libs, extra_cols):
            for c in sorted(ec):
                fh.write(f"{lib}\t{c}\t{ec[c]:.5f}\n")
    print(f"per-BAM: {len(libs)} libraries x {len(mags)} MAGs; {len(extra)} extra contigs (per contig)")
    if emap:
        ev = defaultdict(list)
        with open(emap) as fh:
            for r in csv.DictReader(fh, delimiter="\t"):
                ev[r["event"]].append(r["library"])
        missing = sorted({l for ls in ev.values() for l in ls} - set(libs))
        unmapped = sorted(set(libs) - {l for ls in ev.values() for l in ls})
        events = sorted(ev)
        with open(os.path.join(OUT, "breadth_per_event.tsv"), "w") as fh:
            fh.write("mag\t" + "\t".join(events) + "\n")
            cols = [breadth([per_lib[l] for l in ev[e] if l in per_lib], c2m, mag_len) for e in events]
            for m in mags:
                fh.write(m + "\t" + "\t".join(f"{c[m]:.5f}" for c in cols) + "\n")
        with open(os.path.join(OUT, "breadth_extra_contigs_per_event.tsv"), "w") as fh:
            fh.write("event\tcontig\tbreadth\n")
            for e in events:
                ec = contig_breadth([per_lib[l] for l in ev[e] if l in per_lib], extra)
                for c in sorted(ec):
                    fh.write(f"{e}\t{c}\t{ec[c]:.5f}\n")
        print(f"per-event: {len(events)} events; libraries missing from outputs: {len(missing)}; "
              f"libraries without an event: {len(unmapped)}")


if __name__ == "__main__":
    main()
