"""Write the P6 synthetic-benchmark block map (protocol P6 0a9b44648826, 'Synthetic benchmark').

Reads alife/mep_redox/truth/unblinding_key.tsv and writes results/p6/synthetic_block_map.tsv
with ONLY two columns: dataset_id, block (random label B01..). An emergent run and the
fixed-enzyme yardstick variants derived from it share a block. Nothing else from the key
(scenario names, seeds, emergent/yardstick status) is written or printed.
"""
import csv
import os
import random
import sys

YARDSTICK_SUFFIXES = ("best_found_ref", "mep_opt", "greedy", "random1", "random2")
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
KEY = os.path.join(ROOT, "alife", "mep_redox", "truth", "unblinding_key.tsv")
OUT = os.path.join(ROOT, "results", "p6", "synthetic_block_map.tsv")


def base_run(run_id):
    parts = run_id.split("__")
    if parts[-1] in YARDSTICK_SUFFIXES:
        parts = parts[:-1]
    return "__".join(parts)


def main():
    with open(KEY) as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    bases = sorted({base_run(r["run_id"]) for r in rows})
    rng = random.SystemRandom()
    labels = [f"B{i + 1:02d}" for i in range(len(bases))]
    rng.shuffle(labels)
    base_to_label = dict(zip(bases, labels))
    os.makedirs(os.path.dirname(OUT), exist_ok=True)
    with open(OUT, "w") as fh:
        fh.write("dataset_id\tblock\n")
        for r in sorted(rows, key=lambda r: r["dataset_id"]):
            fh.write(f"{r['dataset_id']}\t{base_to_label[base_run(r['run_id'])]}\n")
    sizes = sorted(sum(1 for r in rows if base_run(r["run_id"]) == b) for b in bases)
    # Only aggregate counts are printed.
    print(f"datasets={len(rows)} blocks={len(bases)} block_sizes={sizes}", file=sys.stderr)


if __name__ == "__main__":
    main()
