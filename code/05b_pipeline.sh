#!/usr/bin/env bash
# ============================================================================
# REAL GBS -> 0/1/2 pipeline on demultiplexed samples (subset, chr01 only).
# Mirrors the paper's logic (they used Cutadapt+Bowtie2+NGSEP); we use the
# equivalent, already-installed tools: fastp + bwa + bcftools.
#   reads -> trim -> align -> sort -> pileup -> genotype call -> 0/1/2 matrix
# ============================================================================
set -euo pipefail
cd "$(dirname "$0")/../bioinfo"

BWA=/Users/md.rasheduzzaman/miniconda3/envs/metagx-bio/bin/bwa
SAMTOOLS=/Users/md.rasheduzzaman/miniconda3/envs/metagx-bio/bin/samtools
FASTP=/Users/md.rasheduzzaman/miniconda3/envs/metagx-bio/bin/fastp
BCFTOOLS=/Users/md.rasheduzzaman/miniconda3/envs/metagx-amr/bin/bcftools
REF=chr01.fa

mkdir -p bam
: > align_stats.tsv
echo -e "sample\traw_reads\ttrimmed\tmapped_chr01\tmap_rate" >> align_stats.tsv

while read S; do
  raw=$(( $(wc -l < samples/$S.fastq) / 4 ))
  # 1) TRIM: quality + adapter (ApeKI 50bp reads); quiet
  $FASTP -i samples/$S.fastq -o samples/$S.trim.fastq \
         -q 20 -l 30 --json /dev/null --html /dev/null 2>/dev/null
  trimmed=$(( $(wc -l < samples/$S.trim.fastq) / 4 ))
  # 2) ALIGN to chr01 with read-group = sample name (so VCF columns are named)
  $BWA mem -R "@RG\tID:$S\tSM:$S" -t 4 $REF samples/$S.trim.fastq 2>/dev/null \
    | $SAMTOOLS sort -o bam/$S.bam -
  $SAMTOOLS index bam/$S.bam
  mapped=$($SAMTOOLS view -c -F 4 bam/$S.bam)
  rate=$(awk -v m=$mapped -v t=$trimmed 'BEGIN{printf "%.1f%%", t? 100*m/t:0}')
  echo -e "$S\t$raw\t$trimmed\t$mapped\t$rate" >> align_stats.tsv
  rm -f samples/$S.trim.fastq
done < samples/samples.txt

echo "=== alignment summary ==="; column -t align_stats.tsv

# 3) JOINT VARIANT CALL across all samples -> multi-sample VCF
$SAMTOOLS faidx $REF
$BCFTOOLS mpileup -f $REF -q 20 -Q 20 -a FORMAT/DP bam/*.bam 2>/dev/null \
  | $BCFTOOLS call -mv -Oz -o calls.vcf.gz 2>/dev/null
$BCFTOOLS index calls.vcf.gz

# 4) FILTER like the paper (biallelic SNPs, drop very low quality)
$BCFTOOLS view -m2 -M2 -v snps -e 'QUAL<30' calls.vcf.gz -Oz -o snps.vcf.gz
nsnp=$($BCFTOOLS view -H snps.vcf.gz | wc -l | tr -d ' ')
echo "=== biallelic SNPs (QUAL>=30) on chr01 subset: $nsnp ==="
echo "=== example genotype calls (CHROM POS REF ALT then per-sample GT) ==="
$BCFTOOLS query -f '%CHROM\t%POS\t%REF\t%ALT[\t%GT]\n' snps.vcf.gz | head -8

du -sh . | awk '{print "bioinfo dir size:",$1}'
echo "PIPELINE_DONE"
