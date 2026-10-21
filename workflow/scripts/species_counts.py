#!/usr/bin/env python3
"""Per-read species assignment from a SAM stream on stdin (no header, supplementary removed).

Records of one read must be consecutive (minimap2 output order, or `samtools collate`).
Each read is assigned to the species of its primary alignment. Reads whose primary and
secondary alignments hit both species are additionally counted as cross-species ambiguous.
"""
import argparse
import sys

ap = argparse.ArgumentParser()
ap.add_argument("--map", required=True, help="TSV: reference name -> species tag")
ap.add_argument("--species", nargs=2, required=True, metavar=("HOST", "XENO"))
ap.add_argument("--sample", required=True)
ap.add_argument("--space", required=True, help="tx or gn")
ap.add_argument("--out", required=True)
args = ap.parse_args()

a, b = args.species
ref2sp = {}
with open(args.map) as fh:
    for line in fh:
        ref, sp = line.rstrip("\n").split("\t")[:2]
        ref2sp[ref] = sp

counts = {"total": 0, "unmapped": 0, a: 0, b: 0, "amb_" + a: 0, "amb_" + b: 0}
unknown = set()
seen_done = 0

cur, prim, hits = None, None, set()


def flush():
    counts["total"] += 1
    if prim is None:
        counts["unmapped"] += 1
        return
    counts[prim] += 1
    if len(hits) > 1:
        counts["amb_" + prim] += 1


for line in sys.stdin:
    qname, flag, rname = line.split("\t", 3)[:3]
    if qname != cur:
        if cur is not None:
            flush()
        cur, prim, hits = qname, None, set()
    flag = int(flag)
    if flag & 4:
        continue
    sp = ref2sp.get(rname)
    if sp is None:
        unknown.add(rname)
        continue
    hits.add(sp)
    if not flag & 256:
        prim = sp
if cur is not None:
    flush()

if unknown:
    sys.exit(f"{len(unknown)} reference(s) missing from species map, e.g. {sorted(unknown)[:3]}")

cols = ["sample", "space", "total", "unmapped", a, b, "amb_" + a, "amb_" + b]
with open(args.out, "w") as fh:
    fh.write("\t".join(cols) + "\n")
    fh.write("\t".join([args.sample, args.space] + [str(counts[c]) for c in cols[2:]]) + "\n")
