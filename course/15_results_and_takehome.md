# Lesson 15 — Results & Take-home Messages

> **The question:** Zoom all the way back out. In one page: *what did the study find, what does it
> mean for a bean breeder, and what general principles should you carry to any genomic-prediction
> problem?*

---

## 15.1 The four objectives, answered

| # | Objective | Answer | Evidence |
|---|-----------|--------|----------|
| 1 | ST vs MT genomic prediction | **MT ≈ ST within a cycle; MT ≫ ST across cycles** (+63% yield, +41% appearance) | Figs. 3–4 (L12, L14) |
| 2 | NIRS index as a secondary trait | **No improvement** (low correlation, intact-seed scan) | Fig. 3, Supp. T3 (L11) |
| 3 | GWAS SNPs as fixed effects | **Lowered accuracy for all traits/years** (unstable hits, polygenic) | Figs. 3–4 (L10) |
| 4 | New-cycle lines needed for updating | **Accuracy rises steadily as new-cycle lines are added** (yield up to ~1.8-fold by 50%) | Fig. 4 (L14) |

---

## 15.2 The numbers worth remembering

- **Prediction is feasible.** Within related material, accuracies ran **0.39–0.93** depending on
  trait; **color/SW** (high $h^2$) reached **~0.93**, **yield/appearance** (low $h^2$) sat
  **~0.4–0.6** — *moderate-to-high, and good enough to select.*
- **KA ≈ GBLUP**, with KA consistently **~+0.03** better → the authors carried both forward.
- **Heritability ran the show:** color (high $h^2$) easy; days-to-maturity's $h^2$ fell 0.77→0.49
  across years and its accuracy fell with it.
- **MT across cycles:** **+63% (yield)**, **+41% (appearance)** — the headline gains.

🔬 Our independent reproductions agree with the paper: GBLUP yield accuracy **0.64** (three ways,
L7); Bonferroni threshold **2.16×10⁻⁵** (L9); high $h^2$ color far more GWAS-detectable than
low $h^2$ yield (105 vs 4 hits, L9); and the across-cycle accuracy curve rising with added new-cycle
lines (L14, `figures/07_across_cycle.png`).

---

## 15.3 The ONE big idea (if you remember nothing else)

> **Extra information helps a prediction model only when it is both *stable* and *genuinely
> correlated* with the target — and you can't fake either by being confident.**

The study is, at heart, three tests of this principle:
- **Correlated traits (SW, Text):** stable ✓ + strongly correlated ✓ → **helped a lot** (MT).
- **NIRS index:** stable ✓ but weakly correlated ✗ → **didn't help**.
- **GWAS hits as fixed effects:** strongly "relevant"-looking but **unstable** ✗ → **hurt**
  (winner's curse, false certainty).

And the corollary that makes GBLUP's humility a virtue:

> **For polygenic traits, treating all markers as small, equally-shrunk, uncertain effects
> (GBLUP) is close to optimal. Injecting false certainty (fixing GWAS hits) degrades it.**

---

## 15.4 The breeder's playbook from this study

1. **Use genomic prediction** — GBS is cheap, accuracies are selection-worthy.
2. **Default to GBLUP or KA**; KA if you want a consistent small edge.
3. **Don't bolt GWAS hits on as fixed effects** for polygenic traits — it backfires.
4. **Do add a cheap, heritable, correlated trait** (seed weight, texture) via a **multi-trait
   model** — especially when predicting **new** material.
5. **Continually update the training set** with a sample of each new cycle's phenotyped lines;
   models decay as germplasm advances.
6. **Set expectations by heritability** — low $h^2$ traits will never predict like high $h^2$ ones.

---

## 15.5 Why the *sequence* of the study makes sense

Trace the logic — each step *had* to come before the next:
1. **Clean phenotypes** (BLUPs, L3) — else you predict field noise.
2. **Quality markers** (L4) → **relatedness G** (L6) — the substrate of prediction.
3. **Establish the GP baseline** (GBLUP/KA, L7–8) — you need a yardstick before testing add-ons.
4. **Test each add-on against that baseline** — GWAS (L9–10), NIRS (L11), multi-trait (L12) — one
   variable at a time, honestly scored (L13).
5. **Stress-test in the real deployment regime** (across cycles, L14) — where the winner (MT) and
   the operational rule (update training) finally emerge.

You can't judge "does GWAS help?" without a baseline; you can't trust any comparison without
honest cross-validation; and you can't claim real-world value without the across-cycle test. The
ordering *is* the scientific method applied to breeding.

---

## 15.6 Limitations the authors are honest about
- Only **one preprocessing** (SNV) for NIRS; derivatives/wavelets might do better.
- **Intact dry seeds** scanned — may miss cooked-quality signal.
- **Conservative Bonferroni** may have missed small-effect loci; a function-aware analysis of loci
  is needed before GWAS-assisted GP could help.
- Modest panel (**415 lines**); larger panels show bigger updating gains.

---

## 15.7 What you should now be able to say
- State all **four objectives and their answers**, and the headline numbers (0.4–0.93 accuracy;
  +63%/+41% MT across cycles; KA +0.03; GWAS/NIRS no help).
- Articulate the **one big idea** (stable **and** correlated information helps; false certainty
  hurts) and the **breeder's playbook**.
- Explain **why the study's steps are ordered** as they are.

👉 Final: **[Lesson 16 — Reproduce It Yourself](16_reproduce_it_yourself.md)** — run the code and
regenerate every figure.
