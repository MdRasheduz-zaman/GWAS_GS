#!/usr/bin/env python3
"""Demultiplex a streamed GBS subset by the study's real barcodes.
Read layout (discovered from the data): [inline barcode][CWGC ApeKI remnant + genomic].
We strip ONLY the barcode (the CWGC is genomic and kept for alignment).
Writes the top-N most-abundant samples to per-sample FASTQs."""
import sys, collections

TOPN = int(sys.argv[1]) if len(sys.argv) > 1 else 8
bcfile = "../repo/GBS_barcode_plate_info_blb.txt"
infq   = "../bioinfo/subset.fastq"
outdir = "../bioinfo/samples"
import os; os.makedirs(outdir, exist_ok=True)

# barcode -> genotype, and N-tolerant lookup keyed by base1+base3..L (ignore failed cycle-2)
bc = {}
for ln in open(bcfile):
    f = ln.rstrip("\n").split("\t")
    if f[0] == "Flowcell" or len(f) < 3: continue
    bc[f[1]] = f[2]
keys = collections.defaultdict(dict)
for b in bc:
    keys[len(b)][b[0] + b[2:]] = b
lengths = sorted(keys)

# pass 1: count reads per barcode
cnt = collections.Counter()
def barcode_of(seq):
    for L in lengths:
        if len(seq) < L + 4: continue
        k = seq[0] + seq[2:L]
        b = keys[L].get(k)
        if b: return b
    return None

n = matched = 0
with open(infq) as fh:
    for h in fh:
        s = fh.readline().rstrip("\n"); fh.readline(); fh.readline()
        n += 1
        b = barcode_of(s)
        if b: cnt[b] += 1; matched += 1
top = [b for b, _ in cnt.most_common(TOPN)]
print(f"reads={n} matched={matched} ({100*matched/n:.1f}%) total_samples={len(cnt)}")
print("Selected top", TOPN, "samples:")
for b in top:
    print(f"  {bc[b]:20s} barcode={b:9s} reads={cnt[b]}")

# pass 2: write per-sample FASTQ, barcode stripped (CWGC remnant kept = genomic)
outfh = {b: open(f"{outdir}/{bc[b]}.fastq", "w") for b in top}
topset = set(top)
written = collections.Counter()
with open(infq) as fh:
    while True:
        h = fh.readline()
        if not h: break
        s = fh.readline(); p = fh.readline(); q = fh.readline()
        seq = s.rstrip("\n")
        b = barcode_of(seq)
        if b in topset:
            L = len(b)
            outfh[b].write(h + seq[L:] + "\n" + p + q[L:])
            written[b] += 1
for f in outfh.values(): f.close()
print("Wrote per-sample FASTQs to", outdir)
# emit sample list for the shell pipeline
with open(f"{outdir}/samples.txt", "w") as f:
    for b in top: f.write(bc[b] + "\n")
