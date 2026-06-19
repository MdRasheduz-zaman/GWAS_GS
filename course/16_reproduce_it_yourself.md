# Lesson 16 — Reproduce It Yourself

> **The question:** You've understood the study. Now *run* it. This lesson is a hands-on guide to
> the scripts in `code/` — what each does, how to run it, and how the original authors' scripts in
> `repo/` map onto the paper. Reproduction is the final test of understanding.

---

## 16.1 Setup (once)

```r
# In R (>= 4.0):
install.packages(c("BGLR", "rrBLUP", "SpATS", "ggplot2"))   # prediction, spatial, figures
# These need a compiler toolchain (clang + gfortran + make), which ships with Xcode
# Command Line Tools on macOS: run  xcode-select --install  once if missing.
```
**Bioinformatics tools** for Lesson 4 (`bwa`, `samtools`, `fastp`) install via conda/bioconda:
`conda install -c bioconda bwa samtools fastp`.

> **Why two packages didn't install the first time (and weren't Mac-specific):** the current CRAN
> build of **`SFSI`** requires a newer R than was on this machine (R 4.3) — a version-availability
> issue, not a macOS issue. So we implement its **regularized selection index from scratch**
> (penalized regression — Lessons 7 & 11) instead. **`GAPIT3`/FarmCPU** needs extra deps; we
> reproduce the GWAS *concept* with a transparent PC-corrected scan (Lesson 9). `SpATS` and `BGLR`
> compiled fine once the toolchain was confirmed.

Most data is already in `repo/GB_BLB.RData` (phenotypes, genotypes, NIRS, SNP positions — Lesson
2). For the **real genotyping pipeline (Lesson 4)** the raw reads are on NCBI BioProject
**PRJNA1138671** / ENA. The full smallest lane is ~10.7 GB (361M reads); we **stream just a
subset** straight from ENA without downloading the whole file:
```bash
# stream ~6M reads (downloads only a fraction, not the full 10.7 GB)
curl -s "https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR299/016/SRR29913416/SRR29913416.fastq.gz" \
  | gunzip | head -24000000 > bioinfo/subset.fastq
```
and align to **one chromosome** to keep disk/time small. The logic is identical at full scale.

Run scripts from inside the `code/` folder (paths are relative):
```bash
cd code
Rscript 01_explore_data.R
```

---

## 16.2 Our reproduction scripts (`code/`)

| Script | Lessons | What it does | Key outputs |
|--------|---------|--------------|-------------|
| `01_explore_data.R` | 1–4 | Loads data; trait distributions, correlations, MAF, SNPs/chromosome | `figures/01_*`, `02_*`, `03_*`; prints yield↔SW = +0.41, SW↔texture = −0.36 |
| `02_gblup_from_scratch.R` | 5–7 | Builds **G**; PCA + kinship; GBLUP **3 ways** (scratch / rrBLUP / BGLR) | `figures/04_pca`, `05_kinship`; **accuracy 0.64** ×3 |
| `03_gwas.R` | 9 | PC-corrected single-marker GWAS; Bonferroni; Manhattan | `figures/06_manhattan`; threshold **2.16×10⁻⁵**; 4 (yield) vs 105 (color) hits |
| `04_across_cycle.R` | 12–14 | Across-cycle ST vs MT; accuracy vs % new-cycle lines added | `figures/07_across_cycle`; rising curves, MT > ST |
| `05a_demux.py` | 4 | Demultiplexes the study's real streamed GBS reads by barcode | per-sample FASTQs; 97.9% barcode match |
| `05b_pipeline.sh` | 4 | Trim → align (chr01) → sort (`fastp`/`bwa`/`samtools`) | `bioinfo/bam/*.bam`, `align_stats.tsv` |
| `05c_call_genotypes.py` | 4 | Pileup → **0/1/2** via allele fraction; het diagnostics | `bioinfo/genotypes_012.tsv`, `figures/08_real_genotype_matrix.png` |
| `06_spatial_demo.R` | 3 | Runs **SpATS**; recovers genotype BLUPs from a noisy field | `figures/09_spats_demo.png`; r 0.46→0.84 |
| `07_heritability_accuracy.R` | 5,13,15 | h²_g vs GBLUP accuracy for all 2018 traits | `figures/10_*`; **r(h², acc)=0.80** |
| `toy_04_reads_to_012.R` | 4 | toy reads→dosage at one SNP | `figures/12_*` |
| `toy_05_breeding_value.R` | 5 | toy 5×10 g = M·α, fully worked | `figures/11_*` |
| `toy_05b_quantgen.R` | 5 | variance components & H²; additive vs dominance (σ²_A/σ²_D); breeder's equation | `figures/23_*`, `24_*` |
| `toy_06_pedigree_kinship.R` | 6 | pedigree A (tabular) vs marker G; Mendelian sampling & cryptic relatedness | `figures/22_*` |
| `toy_08_kernels.R` | 8 | toy distance→Gaussian kernel, bandwidth effect | `figures/15_*` |
| `toy_09_gwas.R` | 9 | toy per-SNP tests → mini Manhattan + Bonferroni | `figures/16_*` |
| `toy_11_rsi.R` | 11 | toy index: OLS overfit vs ridge (p>n) | `figures/17_*` |
| `toy_10_winners_curse.R` | 10 | toy: GWAS hits inflate (winner's curse); fixed < shrink | `figures/20_*` |
| `toy_12_multitrait.R` | 12 | toy: correlated secondary rescues hidden primary | `figures/19_*` |
| `toy_14_across_cycle.R` | 14 | toy: two cycles; accuracy rises as new-cycle lines added | `figures/21_*` |
| `toy_13_cv.R` | 13 | toy hide/predict/correlate + 100-split spread | `figures/18_*` |

(Also `toy` figures 13 = G matrix and 14 = BLUP shrinkage are produced by inline snippets.)

Each script is heavily commented and re-implements the method *from first principles* before (or
instead of) calling a package — so you can read every line. Run them in order; later scripts assume
you understand earlier figures.

🔬 **What "reproduced" means here.** We confirmed our pipeline matches the paper on independently
checkable anchors: GBLUP yield accuracy ≈ **0.64** (paper ~0.60); Bonferroni = **2.16×10⁻⁵**
(paper's exact value); heritability→accuracy and heritability→GWAS-power patterns; and the
across-cycle accuracy-rises-with-updating trend with MT ahead of ST. We deliberately used lighter
settings (fewer CV reps, shorter MCMC) so scripts finish in minutes on a laptop, so absolute
numbers differ slightly from the paper's 100-rep / 12,000-iteration runs — the *patterns* match.

---

## 16.3 The authors' original scripts (`repo/`) — map to the paper

The repo's scripts are written for a computing **cluster** (note the `#SBATCH` headers and the
`SLURM_ARRAY_TASK_ID` job index — each of 1,000 array jobs does one trait × one of 100
repetitions). To run one locally, set the job id manually:

```r
Sys.setenv(SLURM_ARRAY_TASK_ID = "1")   # then source the script; it picks trait 1, rep 1
```

| Repo script | Paper section | What it computes |
|-------------|---------------|------------------|
| `0.Spatial_Analysis.r` | Phenotypic data (L3) | SpATS spatial model → **BLUPs** (the `pheno` table). Needs the raw Excel, not shipped. |
| `1.BLB_GBLUP_RKHS_traits20XX.R` | GBLUP vs kernels (L7–8) | GBLUP + 3 Gaussian kernels (K1/K2/K3) + **KA**, per trait |
| `2.BLB_GBLUP_KA_GWAS_RSI_20XX.R` | RSI, GWAS, MT (L9–13) | builds NIRS **RSI** (SFSI), runs **FarmCPU** GWAS (GAPIT), ST/MT GBLUP & KA, ± GWAS fixed effects |
| `3.BLB_..._YD/app_20XX.R` | Across cycles, MT (L12,L14) | the cross-cycle, multi-trait runs with correlated secondaries (SW for yield, Text for App) and varying % of cycle 2 in training |

🧠 **Reading tip.** Every repo script follows the same skeleton you now recognize: *load data →
build G and distance D → subset to a trait → 70/30 split → fit models → `cor(pred, truth)` → save*.
Once you see that skeleton, the 320-line scripts become readable.

---

## 16.4 A guided mini-exercise (≈ 15 min)

Cement understanding by *changing* things and predicting the result first:

1. **Trait swap.** In `02_gblup_from_scratch.R`, change `trait <- "yd18"` to `"col18"` (color,
   high $h^2$). *Predict:* accuracy should be **much higher** (≈0.9). Run it. Were you right? (This
   is heritability→accuracy from Lesson 5/7, live.)
2. **Break the relatedness.** In `04_across_cycle.R`, set `props <- c(0)` only. *Predict:* the
   hardest case (no new-cycle lines) → lowest accuracy, biggest MT advantage. Confirm.
3. **Kill the secondary correlation.** In `04`, replace `sw <- Y$sw19` with a random vector
   `sw <- rnorm(nrow(Y))`. *Predict:* MT should collapse to ≈ ST (no $\sigma_{g_{12}}$ to exploit —
   Lesson 12). Confirm.

If your predictions matched the runs, you understand the study. If they didn't, the surprise is
exactly where to re-read.

---

## 16.5 Troubleshooting
- **`SFSI`/`GAPIT3` won't install** → they need compilers/extra deps. You don't need them: our
  `03_gwas.R` reproduces the GWAS *concept*, and the RSI is penalized regression you understand
  from Lesson 7/11.
- **BGLR prints pages of iterations** → normal MCMC logging; we wrap it in `sink()` in our scripts.
- **Numbers differ slightly between runs** → BGLR is stochastic and splits are random; set the seed
  (we do) and/or raise repetitions for stability.

---

## 16.6 Where to go next
- Re-read **[Lesson 0's mental map](00_mental_map.md)** — every box should now feel concrete.
- Read the paper (`ref_paper.pdf`) end to end; you have the vocabulary now.
- Try a fourth script: multi-trait for **appearance + texture** (mirror `04` with
  `y <- Y$app19; sw <- Y$text19`) and see if you reproduce the appearance gains.

**You can now read a quantitative-genetics / genomic-prediction paper, reproduce its core, and
explain every step to someone else. That's the whole goal.** 🌱
