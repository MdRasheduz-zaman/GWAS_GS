# Genomic Prediction in Black Beans — A Course for Beginners

> A step-by-step, beginner-friendly course built around one real study:
> **Izquierdo, Wright & Cichy (2025), *G3* — "GWAS-assisted and multitrait genomic
> prediction for improvement of seed yield and canning quality traits in a black
> bean breeding panel."** ([DOI: 10.1093/g3journal/jkaf007](https://doi.org/10.1093/g3journal/jkaf007))

We are going to **reproduce this study** and, along the way, explain *every* moving
part — the **biology**, the **breeding**, the **statistics**, and the **math** — in
language a 2nd/3rd-year BSc student can follow. By the end you will not just "know what
GWAS and genomic selection are"; you will have a **mental map** of the whole research
question and be able to run the analysis yourself.

---

## Who this is for
- BSc students in plant breeding, genetics, agronomy, biology, or data science.
- Anyone who has heard "GWAS", "GBLUP", "genomic selection", "BLUP", "kinship matrix"
  and wants to *really* understand them, not just nod along.

**Assumed background:** high-school algebra, the idea of a mean and a correlation, and a
willingness to meet a little matrix notation (we explain it as we go). No prior genetics
or R required — we build it up.

---

## How to read this course
Each lesson follows the same rhythm:
1. **The question** — what are we trying to find out, and *why now* in the pipeline?
2. **Intuition first** — a picture/analogy before any symbol.
3. **The math, gently** — every symbol defined, small worked numbers.
4. **In the data** — real numbers/figures reproduced from *this* study.
5. **Why they did it this way** — the breeding logic behind the choice.

Boxes you'll see:
- 🧠 **Intuition** — the plain-language mental model.
- 🧮 **Math** — the equation, with every term defined.
- 🧸 **Toy example** — a tiny made-up dataset (a handful of lines/SNPs) where we compute *every
  number by hand and visualize it*, **before** touching the real 415-line data.
- 🔭 **Zoom out** — the bridge: "now picture that toy stretched to the real 415 × 2,315 scale."
- 🌱 **Breeding logic** — why a breeder cares.
- 🔬 **In the data** — a real, reproduced result from the black bean dataset.
- ⚠️ **Common confusion** — a trap to avoid.

> **The core teaching move:** every hard idea appears first on a 🧸 **toy** you can hold in your
> head (e.g. 5 lines × 10 SNPs), then 🔭 **zooms out** to the real data. Toy → real, every time.

---

## Syllabus

| # | Lesson | What you'll be able to do |
|---|--------|---------------------------|
| 0 | [The Mental Map](00_mental_map.md) | Draw the whole study as one flowchart and state the research question |
| 1 | [Beans, Breeding & the Traits](01_biology_and_breeding.md) | Explain the crop, the traits, and *why* yield vs. quality fight each other |
| 2 | [The Data: anatomy of `GB_BLB`](02_the_data.md) | Open the real dataset and know what every piece is |
| 3 | [Phenotyping & Spatial BLUPs](03_phenotyping_spatial_BLUPs.md) | Turn messy field plots into one clean number per line; understand heritability |
| 4 | [Genotyping by Sequencing & SNPs](04_genotyping_GBS_SNPs.md) | Explain how DNA becomes a 0/1/2 marker matrix, and every filter applied |
| 5 | [Quantitative Genetics Foundations](05_quant_genetics_foundations.md) | Define breeding value, additive effect, and heritability with equations |
| 6 | [The Genomic Relationship Matrix **G**](06_genomic_relationship_matrix.md) | Build VanRaden's **G** by hand and read a kinship heatmap |
| 7 | [GBLUP — genomic prediction's workhorse](07_GBLUP.md) | Derive GBLUP, see it = ridge regression, and reproduce accuracy 0.64 |
| 8 | [RKHS & Gaussian Kernels](08_RKHS_kernels.md) | Explain kernels, bandwidth, and kernel averaging (KA) |
| 9 | [GWAS with FarmCPU](09_GWAS_FarmCPU.md) | Run/understand a GWAS, PCs, Bonferroni, Manhattan plots |
| 10 | [GWAS-assisted Genomic Prediction](10_GWAS_assisted_GP.md) | Add markers as fixed effects — and explain why it *hurt* here |
| 11 | [NIRS & Regularized Selection Indices](11_NIRS_RSI.md) | Use spectra as a cheap secondary trait via penalized regression |
| 12 | [Multi-Trait Models](12_multitrait_models.md) | Borrow strength from correlated traits; why MT beat ST across cycles |
| 13 | [Cross-Validation & Prediction Accuracy](13_cross_validation_accuracy.md) | Measure "how good is the prediction?" honestly; read Table 1 |
| 14 | [Across Breeding Cycles & Updating](14_across_cycles_updating.md) | Predict a *new* cycle and explain why the training set must keep growing |
| 15 | [Results & Take-home Messages](15_results_and_takehome.md) | Summarize what the study found and what it means for breeders |
| 16 | [Reproduce It Yourself](16_reproduce_it_yourself.md) | Run the scripts in `code/` and regenerate the figures |

---

## What's in this folder

```
GWAS_GS/
├── ref_paper.pdf          # the study (G3, 2025)
├── repo/                  # the authors' original code + data (cloned from GitHub)
│   ├── GB_BLB.RData       # phenotypes, genotypes, NIRS, SNP positions
│   └── *.R                # their analysis scripts
├── course/                # ← you are here: the lessons
├── code/                  # our clean, runnable reproduction scripts
└── figures/               # figures we generate from the real data
```

> **Reproducibility note.** Every number quoted in 🔬 *In the data* boxes is produced by
> the scripts in `code/`. We re-implement the core methods *from scratch* in base R first
> (so you can see the machinery), then confirm against the same packages the authors used
> (`BGLR`, `rrBLUP`). Where we ran the real reproduction, you'll see the exact figure name.

Start with **[Lesson 0: The Mental Map](00_mental_map.md).**
