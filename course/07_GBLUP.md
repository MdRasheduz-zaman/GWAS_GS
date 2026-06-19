# Lesson 7 — GBLUP: Genomic Prediction's Workhorse

> **The question:** We have clean phenotypes (BLUPs, L3), a marker matrix (L4), a relatedness
> matrix **G** (L6), and the idea that prediction = exploiting relatedness (L5). **GBLUP** is
> the model that turns all of that into an actual predicted number for every line — including
> lines we never phenotyped. This is the core of the entire paper.

---

## 7.1 The model in one line

🧮 **GBLUP model.**
$$
\mathbf{y} = \mathbf{1}\mu + \mathbf{g} + \mathbf{e}, \qquad
\mathbf{g} \sim N(\mathbf 0,\ \mathbf G\,\sigma_g^2), \qquad
\mathbf{e} \sim N(\mathbf 0,\ \mathbf I\,\sigma_e^2)
$$
- $\mathbf{y}$ — vector of phenotypes (BLUPs) for the lines we *did* measure.
- $\mu$ — overall mean; $\mathbf 1$ — a column of ones.
- $\mathbf{g}$ — the breeding values we want (Lesson 5), one per line.
- **The key line:** $\mathbf{g} \sim N(\mathbf 0, \mathbf G \sigma_g^2)$ says *breeding values
  are correlated exactly according to G* (Lesson 6). Relatives have correlated genetic values.
- $\mathbf{e}$ — independent environmental noise.

🧠 **Intuition.** We're telling the model: "Genetic merit is a smooth quantity over the
relatedness map G. If line A (unmeasured) sits genetically near lines B, C, D (measured, all
high-yielding), then A is probably high-yielding too." GBLUP formalizes "you resemble your
relatives" into optimal arithmetic.

---

## 7.2 How it predicts an *unseen* line (the part that feels magic)

Split lines into **training** (phenotyped) and **test** (genotyped only). Because all breeding
values are jointly normal with covariance G, the prediction for the test lines is the textbook
conditional-mean formula:

🧮
$$
\hat{\mathbf{g}}_{\text{test}} \;=\; \mathbf{G}_{\text{test,train}}\,\big(\mathbf{G}_{\text{train,train}} + \lambda \mathbf{I}\big)^{-1}\,(\mathbf{y}_{\text{train}}-\mu), \qquad \lambda = \frac{\sigma_e^2}{\sigma_g^2}
$$

Decode it:
- $\mathbf{G}_{\text{test,train}}$ — **how related each test line is to each training line** (the
  bridge across the gap).
- $(\mathbf{G}_{\text{train,train}}+\lambda\mathbf I)^{-1}(\mathbf y_{\text{train}}-\mu)$ — turns
  training phenotypes into weights, *down-weighting* noise via $\lambda$.
- Multiply: each test line's prediction = a **relatedness-weighted blend of training
  phenotypes.** A test line gets pulled toward the phenotypes of the training lines it most
  resembles. **That is the whole mechanism.**

⚠️ **Common confusion — "but the test line has no phenotype!"** Correct — and it doesn't need
one. Its *genotype* tells us its relatedness ($\mathbf{G}_{\text{test,train}}$) to lines that
*do* have phenotypes, and relatedness is enough. This is why a single leaf sample suffices to
predict a line you never grew to maturity.

The ratio $\lambda = \sigma_e^2/\sigma_g^2$ is the inverse signal-to-noise — high noise → large
$\lambda$ → more shrinkage. (Same quantity that set BLUP shrinkage in L3 and equals
$(1-h^2)/h^2$.)

### 🧸 Toy — predict a hidden line from its relatives (5 lines)

Take 5 lines; we know the phenotype of four and **hide L3**. Using the toy G (entries =
relatedness) and the formula above with $\lambda=0.5$:

| | L1 | L2 | **L3** | L4 | L5 |
|---|---|---|---|---|---|
| phenotype $y$ | 10 | 8 | **? (hidden)** | 2 | 9 |
| relatedness to L3, $G_{\text{L3},\cdot}$ | −0.65 | −0.15 | (1.35) | **+0.35** | −0.90 |

L3 is **most related to L4** (+0.35, and L4 is the *low* performer, $y=2$) and **genetically
opposite** the high performers L1 ($y=10$) and L5 ($y=9$). So the relatedness-weighted blend pulls
L3's prediction *below* the mean ($\mu=7.25$):
$$
\hat y_{\text{L3}} = \mu + \mathbf G_{\text{L3,trn}}(\mathbf G_{\text{trn,trn}}+\lambda\mathbf I)^{-1}(\mathbf y_{\text{trn}}-\mu) = \mathbf{4.27}
$$
🔭 **Zoom out:** the real GBLUP does this for each of 143 new-cycle lines at once, blending **272**
training phenotypes by a **415×415** relatedness matrix. Same formula — we *just verified it gives
0.64 accuracy on real yield, three different ways* (§7.4).

---

## 7.3 The beautiful duality: GBLUP **=** ridge regression on markers

Recall Lesson 5's two routes. Here's why they're the same.

**Route A (RR-BLUP):** estimate every marker effect $\alpha_j$, but penalize their size:
$$
\hat{\boldsymbol\alpha} = \arg\min_{\boldsymbol\alpha} \ \|\mathbf y - \mathbf Z\boldsymbol\alpha\|^2 + \lambda\|\boldsymbol\alpha\|^2
$$
then $\hat{\mathbf g} = \mathbf Z\hat{\boldsymbol\alpha}$. The $\lambda\|\boldsymbol\alpha\|^2$
term — **ridge** — is what makes the impossible $p\gg n$ problem solvable: it forbids any single
marker from taking a wild effect, spreading signal across all of them.

**Route B (GBLUP):** the conditional-mean formula in §7.2 using $\mathbf G = \mathbf Z\mathbf
Z^\top/p$.

🧮 **They are algebraically identical predictions.** Substituting $\mathbf G=\mathbf Z\mathbf
Z^\top/p$ into Route B reproduces Route A exactly (a standard linear-algebra identity, the
"kernel trick"). So:

> **GBLUP = ridge regression where every SNP is shrunk equally.**
> "Predict from relatedness (G)" and "shrink all marker effects" are *the same model* viewed
> from two ends.

🧠 **Why this matters pedagogically.** It dissolves the apparent mystery of G. GBLUP isn't a
strange new thing — it's penalized regression you could write in one line, repackaged so you
work with a 415×415 matrix instead of 2,315 effects (much faster when $p\gg n$). The assumption
hiding inside is **"all markers contribute small, equal-variance effects"** — a great fit for
**polygenic** traits like yield, and the reason GBLUP shines on them.

---

## 7.4 We reproduced it three ways — and they agree

To *prove* GBLUP is "just the mixed model above," `code/02_gblup_from_scratch.R` predicts 2018
yield (`yd18`, n = 272, 70/30 split) three ways:

| Implementation | What it is | Accuracy (cor(pred, truth) on test) |
|----------------|-----------|--------------------------------------|
| **(A) From scratch** | our own base-R mixed-model solver (the §7.2 formula, variance components by REML) | **0.640** |
| **(B) `rrBLUP::mixed.solve`** | a standard package, REML | **0.640** |
| **(C) `BGLR` RKHS(G)** | the **exact** Bayesian tool the authors used | **0.635** |

🔬 They match. The tiny BGLR difference is Monte-Carlo noise (it samples; the others solve in
closed form). **Take-home: the "black box" in the paper is reproducible from first principles —
you now understand every step that produced 0.64.** And 0.64 lines up with the paper's reported
single-trait yield accuracy (~0.60 in 2018), confirming our pipeline matches theirs.

Heart of the from-scratch code (the §7.2 formula, lightly annotated):
```r
# variance ratio lambda estimated by REML, then:
Ginv     <- solve(G[trn,trn] + lambda * diag(length(trn)))   # (G_train + λI)^-1
ghat_all <- G[, trn] %*% (Ginv %*% (y[trn] - mu))            # G_test,train  ·  weights
pred     <- mu + ghat_all                                     # predictions for ALL lines
accuracy <- cor(y[tst], pred[tst])                            # score on hidden test lines
```

---

## 7.5 What "prediction accuracy" means here (preview of L13)

We score a model by **the correlation between predicted and observed values on the held-out
test lines**. Accuracy 0.64 means predictions and reality line up moderately well — good enough
to *rank* lines and select the top ones, even if not perfect point-by-point.

⚠️ **Common confusion — accuracy vs. $R^2$ vs. heritability.** Accuracy here = Pearson $r$ (not
$r^2$). It is bounded above by $\sqrt{h^2}$-ish quantities: you can't predict noise, so a
trait's heritability caps how high $r$ can go. That's why high-$h^2$ color hits 0.93 and lower-
$h^2$ yield sits near 0.6.

---

## 7.6 Where GBLUP sits in the study

GBLUP is the **baseline** every other model is compared against:
- **vs. Gaussian kernels (KA, Lesson 8):** can a *non-linear* relatedness model beat additive G?
  (Slightly, yes.)
- **vs. GWAS-assisted (Lesson 10):** does adding top SNPs as fixed effects beat plain GBLUP?
  (No — it hurt.)
- **vs. multi-trait (Lesson 12):** does borrowing a correlated trait beat single-trait GBLUP?
  (Across cycles, yes — a lot.)

So master GBLUP and the rest of the paper is "GBLUP, plus one twist at a time."

---

## 7.7 What you should now be able to say
- **GBLUP** models $\mathbf y=\mathbf 1\mu+\mathbf g+\mathbf e$ with $\mathbf g\sim N(\mathbf 0,
  \mathbf G\sigma_g^2)$; it predicts unseen lines as a **relatedness-weighted blend of training
  phenotypes** via $\hat{\mathbf g}_{\text{test}}=\mathbf G_{\text{test,train}}(\mathbf
  G_{\text{train,train}}+\lambda\mathbf I)^{-1}(\mathbf y-\mu)$.
- It needs **no phenotype** on the test line — only its **genotype-based relatedness**.
- **GBLUP = ridge regression on markers** (equal shrinkage); the assumption is *many small equal
  effects*, ideal for polygenic traits.
- We reproduced **accuracy 0.64** for yield three independent ways — the method is fully
  demystified.

👉 Next: **[Lesson 8 — RKHS & Gaussian Kernels](08_RKHS_kernels.md)** — replacing the *linear*
relatedness in G with a flexible *non-linear* one.
