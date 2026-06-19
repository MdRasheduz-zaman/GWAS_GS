#!/usr/bin/env python3
"""Transparent genotype caller: pileup -> 0/1/2 dosage matrix.
This is the heart of "how 0/1/2 is made". Real callers (NGSEP, bcftools) use
genotype-likelihood models; we use the same underlying logic in its simplest,
readable form: at each position, per sample, count REF vs ALT supporting reads,
then call the dosage from the allele fraction.
  hom-ref (0/0)=0   het (0/1)=1   hom-alt (1/1)=2   (low depth = missing)
"""
import re, sys

MINDP   = 5      # need >=5 reads to call a sample
HET_LO, HET_HI = 0.15, 0.85   # alt-fraction thresholds (inbred beans -> mostly 0 or 2)
MIN_CALLED = 6   # keep SNP only if >=6 of 8 samples genotyped
MIN_MAF = 0.05

pile = "../bioinfo/pileup.txt"
samples = [l.strip() for l in open("../bioinfo/samples/samples.txt")]
nS = len(samples)

indel = re.compile(r'[+-](\d+)')
def parse_bases(s, ref):
    """return (ref_count, {alt_base: count})"""
    i = 0; refc = 0; alt = {}
    s = s.upper()
    while i < len(s):
        c = s[i]
        if c == '^': i += 2; continue        # read start + mapQ
        if c == '$': i += 1; continue         # read end
        if c in '+-':                          # indel: skip its bases
            m = indel.match(s, i); n = int(m.group(1))
            i = m.end() + n; continue
        if c in '.,': refc += 1
        elif c in 'ACGT': alt[c] = alt.get(c, 0) + 1
        # '*' deletion / 'N' ignored
        i += 1
    return refc, alt

rows = []   # (snp_id, ref, alt, [dosages])
for ln in open(pile):
    f = ln.rstrip("\n").split("\t")
    chrom, pos, ref = f[0], f[1], f[2].upper()
    if ref not in "ACGT": continue
    per = []          # (refc, altdict) per sample
    altpool = {}
    for k in range(nS):
        dp = int(f[3 + 3*k]); bases = f[4 + 3*k]
        rc, ad = parse_bases(bases, ref) if dp else (0, {})
        per.append((rc, ad))
        for b, n in ad.items(): altpool[b] = altpool.get(b, 0) + n
    if not altpool: continue                       # invariant site
    altb = max(altpool, key=altpool.get)            # main ALT allele
    dosages = []; called = 0; altalleles = 0; totalleles = 0
    for rc, ad in per:
        ac = ad.get(altb, 0); dp = rc + ac
        if dp < MINDP: dosages.append("NA"); continue
        frac = ac / dp
        g = 0 if frac <= HET_LO else (2 if frac >= HET_HI else 1)
        dosages.append(str(g)); called += 1
        altalleles += g; totalleles += 2
    if called < MIN_CALLED: continue
    maf = min(altalleles, totalleles - altalleles) / totalleles
    if maf < MIN_MAF: continue
    rows.append((f"S01_{pos}", ref, altb, dosages))

print(f"samples: {nS}   variable biallelic SNPs called: {len(rows)}")
# write 0/1/2 matrix (samples x SNPs) in the SAME shape as GB_BLB$geno
with open("../bioinfo/genotypes_012.tsv", "w") as o:
    o.write("sample\t" + "\t".join(r[0] for r in rows) + "\n")
    for si, s in enumerate(samples):
        o.write(s + "\t" + "\t".join(r[3][si] for r in rows) + "\n")

# quick stats + a preview identical in spirit to GB_BLB$geno[1:8,1:6]
flat = [g for r in rows for g in r[3] if g != "NA"]
from collections import Counter
c = Counter(flat)
tot = sum(c.values())
print("genotype composition:  0(homref)=%d  1(het)=%d  2(homalt)=%d   het-rate=%.1f%%"
      % (c['0'], c['1'], c['2'], 100*c['1']/tot))
print("\nPreview (first 8 SNPs x first 4 samples) -- this IS the 0/1/2 matrix:")
hdr = [r[0] for r in rows[:8]]
print("%-12s " % "sample" + " ".join("%-11s" % h for h in hdr))
for si in range(min(4, nS)):
    print("%-12s " % samples[si] + " ".join("%-11s" % rows[j][3][si] for j in range(min(8,len(rows)))))
print("\nWrote ../bioinfo/genotypes_012.tsv  (samples x SNPs, 0/1/2 + NA)")
