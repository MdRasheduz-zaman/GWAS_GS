# Lesson 2 — The Data: Anatomy of `GB_BLB`

> **The question:** before any modeling, *what exactly do we have?* A model is only as
> good as your understanding of its inputs. We open the real file and name every part.

---

## 2.1 One file to rule them all

The authors ship a single R data file, `repo/GB_BLB.RData`. Loading it gives one object,
`GB_BLB`, which is a **list** (think: a labeled drawer with 5 compartments):

```r
load("repo/GB_BLB.RData")
names(GB_BLB)
#> "pheno"  "geno"  "nirs18"  "nirs19"  "SNP_Position"
```

🔬 **In the data** (from `code/01_explore_data.R`), the five compartments are:

| Compartment | Shape | What it holds |
|-------------|-------|---------------|
| `pheno` | **415 × 21** | one row per line; the cleaned trait values (BLUPs) for 2018 & 2019 |
| `geno` | **415 × 2315** | one row per line; **2,315 SNP markers** coded 0/1/2 |
| `nirs18` | **415 × 551** | near-infrared **spectrum** per line, 2018 (one column per wavelength) |
| `nirs19` | **415 × 551** | same, 2019 |
| `SNP_Position` | **2315 × 3** | for each SNP: its name, chromosome, and base-pair position |

The magic that makes everything work: **all matrices share the same 415 rows, in the same
order.** Row 7 of `pheno`, row 7 of `geno`, row 7 of `nirs18` are *the same bean line*. This
alignment is what lets a model connect "this DNA" ↔ "this phenotype."

⚠️ **Common confusion — "415 rows but only 272 have 2018 yield?"** Yes. The matrix has a row
for *every* line, but cycle-2 lines (rows 273–415) were not grown in 2018, so their 2018
cells are `NA` (missing). Missingness is *information*, not an error — it encodes the
breeding-cycle design (Lesson 1).

---

## 2.2 `pheno` — the trait table

21 columns: a `taxa` (line name) column plus 20 trait columns. The naming convention is
`trait` + `year`:

```
taxa | yd18 yd19 | sw18 sw19 | dm18 dm19 | df18 df19 | app18 app19 |
       text18 text19 | col18 col19 | l18 l19 | a18 a19 | b18 b19
```

So `yd18` = yield in 2018, `app19` = canning appearance in 2019, etc. These are **not** raw
field measurements — they are **BLUPs**: one cleaned, field-noise-adjusted value per line per
year. *How* raw plot data became BLUPs is the entire subject of **Lesson 3**.

🧠 **Why store BLUPs, not raw data?** Because every downstream model wants *one number per
line* ("the genetic value we're trying to predict"), not a tangle of replicate plots
contaminated by which corner of the field they sat in.

---

## 2.3 `geno` — the marker matrix

- 2,315 columns, each a **SNP** (single-nucleotide polymorphism — a single DNA position where
  lines differ, e.g. some have an **A**, others a **G**).
- Column names like `S01_79737` mean *chromosome 01, position 79,737 bp*.
- Values are **0, 1, 2** = the count of one of the two alleles a line carries at that SNP.
  - For a homozygous self-pollinating bean, you mostly see **0** (zero copies) or **2** (two
    copies); **1** (heterozygous) is rarer.

🧮 **What the numbers mean.** Pick a SNP with alleles A/G. Code the "G" allele:
- `AA` → 0 copies of G → **0**
- `AG` → 1 copy → **1**
- `GG` → 2 copies → **2**

This 0/1/2 encoding turns DNA into something we can do arithmetic on — average it, correlate
it, multiply it by effect sizes. That arithmetic *is* genomic prediction.

> 🔬 **In the data.** Allele coding range is exactly 0–2; SNPs are spread across all **11
> chromosomes** (166–273 SNPs each). Minor-allele frequencies (Lesson 4) range from ~0.007 up
> to 0.5. See `figures/03_maf_histogram.png`.

---

## 2.4 `nirs18` / `nirs19` — the light fingerprint

Near-infrared spectroscopy (NIRS) shines infrared light (1100–2200 nm) at dry beans and
records how much is absorbed at each wavelength. The result is a **spectrum**: 551 numbers
per line (one every 2 nm). Different chemistry (water, starch, protein, seed-coat compounds)
absorbs differently, so the spectrum is an indirect *chemical fingerprint* of the bean.

🌱 **Breeding logic.** A NIRS scan is **fast, cheap, and non-destructive** — you keep the
seed. *If* the spectrum carries information about canning quality, it could be a cheap stand-in
for the expensive taste panel. Lesson 11 tests exactly this idea (spoiler: it didn't help
much here, but it could help in other studies and the *method* is important to understand).

⚠️ The column names differ slightly between years (`nm18_1100` vs `nm19_nm19_1100`) — a tiny
real-world data-cleaning wrinkle you'd fix before combining them.

---

## 2.5 `SNP_Position` — the map

A 3-column table — `SNP`, `Chromosome`, `Position` — giving each marker's physical location
on the *P. vulgaris* v2.1 reference genome. You need this for two things:
1. **GWAS** (Lesson 9): to draw a **Manhattan plot** you must know where each SNP sits.
2. **Interpreting hits**: a significant SNP at chr03:5,000,000 can be looked up against known
   genes.

---

## 2.6 What's *not* in the file (and where it lives)

- **Raw plot data** (individual field plots with row/column positions) — used by the spatial
  analysis (Lesson 3). The repo's `0.Spatial_Analysis.r` reads a separate Excel file
  (`BBL_phenotype_2020.xlsx`) not included here; the *output* of that step (the BLUPs) is what
  lives in `pheno`. So we **explain** Lesson 3 fully but reproduce it conceptually.
- **Raw DNA sequence reads** — these are huge; they live on NCBI (BioProject **PRJNA1138671**)
  and were processed into `geno` by the bioinformatics pipeline in Lesson 4.

---

## 2.7 A 30-second hands-on

```r
load("repo/GB_BLB.RData")
str(GB_BLB, max.level = 1)            # see the 5 compartments
GB_BLB$geno[1:5, 1:6]                 # a corner of the marker matrix (0/1/2)
summary(GB_BLB$pheno$yd18)            # distribution of 2018 yield BLUPs
sum(!is.na(GB_BLB$pheno$yd18))        # 272  -> cycle-1 lines
```

You now know what every input is. Next we earn the right to call `pheno` "clean".

👉 Next: **[Lesson 3 — Phenotyping & Spatial BLUPs](03_phenotyping_spatial_BLUPs.md).**
