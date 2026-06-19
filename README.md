# Genomic Prediction in Black Beans — Reproduce & Learn

A complete, beginner-friendly **course + reproduction** built around the study:

> **Izquierdo, Wright & Cichy (2025).** *GWAS-assisted and multitrait genomic prediction for
> improvement of seed yield and canning quality traits in a black bean breeding panel.*
> **G3** 15(3): jkaf007. [doi:10.1093/g3journal/jkaf007](https://doi.org/10.1093/g3journal/jkaf007)

The goal: teach every **biological, breeding, statistical, and mathematical** idea in this paper
to a BSc-level reader, *while reproducing the analysis from the real data.*

## Start here
- 📖 **Read online:** **[`course/README.md`](course/README.md)** — the syllabus and Lesson 0 (the mental map).
- 📄 **Single PDF:** **`course.pdf`** — the whole course (75 pages: cover, clickable TOC, all 17
  lessons + glossary, with math, diagrams, and figures rendered).

The course is 16 lessons + a glossary, each grounded in numbers/figures we reproduced from the
authors' own data.

> **Rebuild the PDF:** `node code/build_pdf.js` (creates `course.html`), then print it with headless
> Chrome to `course.pdf` (see `code/build_pdf.js` header). Requires `markdown-it` (installed in `code/`).

## Layout
```
GWAS_GS/
├── ref_paper.pdf      # the study
├── repo/              # authors' original scripts + data (GB_BLB.RData), cloned from GitHub
├── course/            # 16 lessons (00–16) + GLOSSARY.md   ← the teaching material
├── code/              # our clean, from-scratch reproduction scripts (01–04)
└── figures/           # figures generated from the real data
```

## What we reproduced (independently verified against the paper)
| Result | Ours | Paper | Lesson |
|--------|------|-------|--------|
| **Raw reads → 0/1/2 matrix** (real GBS pipeline on the study's own SRA data) | **built it** | NGSEP pipeline | 4 |
| GBLUP yield accuracy (3 implementations agree) | **0.64** | ~0.60 | 7 |
| Bonferroni GWAS threshold | **2.16×10⁻⁵** | 2.16×10⁻⁵ | 9 |
| Heritability → GWAS power (color vs yield hits) | **105 vs 4** | high-vs-low $h^2$ | 9 |
| Heritability caps accuracy across traits | **r(h²,acc)=0.80** | core claim | 5 |
| Trait antagonism (seed weight ↔ texture) | **−0.36** | negative | 1 |
| SpATS removes field noise → recovers genetics (BLUP) | **r 0.46→0.84** | method | 3 |
| Across-cycle accuracy rises as new lines added; MT > ST | **reproduced** | +63% yield | 14 |

## Run it
```bash
cd code
Rscript 01_explore_data.R        # data, distributions, correlations, MAF
Rscript 02_gblup_from_scratch.R  # G matrix, PCA, GBLUP three ways
Rscript 03_gwas.R                # GWAS Manhattan, Bonferroni
Rscript 04_across_cycle.R        # across-cycle ST vs MT (the headline finding)
Rscript 06_spatial_demo.R        # SpATS field-noise removal (real run)
# Real GBS reads -> 0/1/2 (needs bwa, samtools, fastp; ~1GB disk):
python3 05a_demux.py 8 && bash 05b_pipeline.sh && python3 05c_call_genotypes.py
```
Requires R ≥ 4.0 with `BGLR`, `rrBLUP`, `SpATS`, `ggplot2`. The bioinformatics demo uses
`bwa`/`samtools`/`fastp`. See **[Lesson 16](course/16_reproduce_it_yourself.md)** for the full
guide, including how to run the authors' original cluster scripts.

## The one big idea
> Extra information (GWAS hits, NIRS, correlated traits) helps a genomic-prediction model **only
> when it is both stable and genuinely correlated** with the target. Correlated traits (multi-trait
> models) delivered; GWAS-as-fixed-effects and NIRS did not. For polygenic traits, GBLUP's humble
> "shrink everything equally" is hard to beat.
