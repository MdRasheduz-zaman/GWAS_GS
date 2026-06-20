# Lesson 6 — The Genomic Relationship Matrix **G**

> **The question:** Lesson 5 said prediction = *exploiting relatedness*. So we need a precise,
> number-for-number answer to **"how genetically related is line $i$ to line $j$?"** for all
> pairs. That answer is the matrix **G**. We build it by hand on 4 lines, then on the real 415.

---

## 6.1 What G is

**G** is a square table, one row and one column per line (415 × 415 here). Entry $G_{ij}$
measures the **genomic relationship** (realized kinship) between lines $i$ and $j$:

- **Large positive** $G_{ij}$ → lines share many alleles → closely related.
- **Near zero** → unrelated.
- **Negative** → *less* related than two random lines (more genetically opposite than average).
- The **diagonal** $G_{ii}$ ≈ how inbred/typical line $i$ is (≈1 on average).

🧠 **Intuition.** G is a *friendship-by-DNA* map. Pedigree breeders drew family trees to guess
relatedness (e.g. "these two are cousins → ~12.5% shared"). Markers let us **measure** it directly —
and catch surprises a pedigree misses (two "unrelated" lines that happen to share a lot, or full
sibs that drifted apart). G is the *data-driven* upgrade of the pedigree relationship matrix.

---

## 6.2 The formula (VanRaden 2008) — exactly what the paper used

🧮 **Step by step.** Start from the marker matrix $\mathbf{M}$ (415 × 2,315, entries 0/1/2).

1. **Center & scale each SNP column** → matrix $\mathbf{Z}$. For SNP $j$ with mean $2p_j$ and
   SD $s_j$:

```math
Z_{ij} = \frac{M_{ij} - 2p_j}{s_j}
```

   Centering removes the "average dose" so we measure *deviation* from the population; scaling
   puts every SNP on a comparable footing.
   
2. **Cross-multiply and average over markers:**

```math
\boxed{ \mathbf{G} = \frac{\mathbf{Z}\mathbf{Z}^{\top}}{p} } \qquad (p = \text{number of markers})
```

In R this is literally the line the authors wrote (and we reproduced):
```r
Z <- scale(M)              # center + scale every SNP
G <- tcrossprod(Z) / ncol(M)   # Z Z' / p
```

🧠 **Why this *is* relatedness.** Writing $i,k$ for the **two lines** and $j$ for the **SNP**
(the same $j$ as in $Z_{ij}$, so the summation runs over markers), one entry is

```math
G_{ik} = \frac{1}{p}\sum_{j=1}^{p} Z_{ij}\,Z_{kj}
```

— just the **average product of standardized genotypes** across all SNPs, i.e. how *correlated*
lines $i$ and $k$ are across the genome. If they tend to be high together and low together at the
same SNPs, the products are positive and pile up → big $G_{ik}$. It's a correlation in disguise.

> 📝 **On the indices.** In $M_{ij}$ and $Z_{ij}$ the second subscript $j$ is the **SNP**, so a
> single relatedness entry needs a *third* letter for the second line — hence $G_{ik}$, summed over
> markers $j$. (Earlier, $G_{ij}$ in §6.1 just meant "the entry in row $i$, column $j$ of the
> table"; here we make the marker sum explicit, so the line pair is $i,k$.)

---

## 6.3 Build it by hand — 4 lines, 3 SNPs

Take a toy with 4 lines and 3 SNPs (we ran this in R):

```
M (0/1/2):              Z (centered+scaled):        G = ZZ'/p:
   SNP1 SNP2 SNP3           SNP1  SNP2  SNP3            L1    L2    L3    L4
L1   2    2    2      L1   0.87  1.22  0.87     L1   1.00  0.50 -0.50 -1.00
L2   2    1    2      L2   0.87  0.00  0.87     L2   0.50  0.50 -0.50 -0.50
L3   1    1    0      L3  -0.87  0.00 -0.87     L3  -0.50 -0.50  0.50  0.50
L4   1    0    0      L4  -0.87 -1.22 -0.87     L4  -1.00 -0.50  0.50  1.00
```

The three SNPs deliberately have **different allele frequencies**, so each gets its own mean
$2p_j$ and SD $s_j$ — and a heterozygote (dose 1) is *not* always "neutral":

| SNP $j$ | doses (col) | mean $2p_j$ | $\sum_i (M_{ij}-2p_j)^2$ | $\div(n{-}1)=\div3$ | $s_j=\sqrt{\,}$ |
|---|---|---|---|---|---|
| 1 | 2,2,1,1 | 1.50 | 1.00 | 0.333 | **0.58** |
| 2 | 2,1,1,0 | 1.00 | 2.00 | 0.667 | **0.82** |
| 3 | 2,2,0,0 | 1.00 | 4.00 | 1.333 | **1.15** |

(The same mean can give different $s_j$: SNP2 and SNP3 both center at 1.0, but SNP3's all-or-nothing
doses spread wider, so $s_3>s_2$.)

🖼️ **Watch the two equations run on this toy** — each symbol below points to a piece of the toy
above. Step 1 shows **how $s_j$ is built for every SNP**, then turns a genotype $M_{ij}$ into a
standardized $Z_{ij}$ (one cell worked out, plus the heterozygote cases). Step 2 takes **two lines'
$Z$-rows** ($i=1$ and $k=2$), multiplies them SNP-by-SNP, sums over the $p=3$ markers $j$, and
divides by $p$ — giving the single entry $G_{12}=0.50$. The extreme case ($G_{14}=-1.00$) shows what
"opposite at every SNP" looks like:

![building G from the toy, step by step](../figures/35_G_construction.png)

So $G_{ik}=\frac{1}{p}\sum_j Z_{ij}Z_{kj}$ is literally "**slide line $i$'s standardized row
against line $k$'s, multiply, average**." Now read the result like a breeder:
- **$G_{12}=0.50$** (L1 vs L2): they match at SNP1 and SNP3 (both 2); at SNP2 L2 is a
  **heterozygote** (dose 1) which sits *exactly* at the mean ($2p_2=1$) → $Z=0$ → that SNP adds
  nothing → still clearly related. ✔
- **$G_{14}=-1.00$** (L1 vs L4): opposite at *every* SNP → maximally *un*alike. ✔
- **$G_{11}=1.00$** (diagonal): L1's relationship with itself = how far it sits from the
  population average across SNPs.
- ⚠️ **A heterozygote isn't automatically $Z=0$.** At SNP1 ($2p_1=1.5$) a dose-1 line (L3, L4)
  gives $Z=(1-1.5)/0.58=-0.87$ — *below* average — because the population leans toward dose 2.

That's the whole idea — the colour heatmap on the right of the figure above *is* this toy G.

🔭 **Zoom out:** the real G is this exact picture at **415 × 415** — a giant heatmap whose red
blocks are families of related lines and blue regions are unrelated pairs (we draw the real one in
§6.4). The real G just does this over 2,315 SNPs and 415 lines.

> 🔬 **In the data.** For the real matrix (`code/02_gblup_from_scratch.R`): mean diagonal ≈
> **0.998** (≈1, as designed) and mean off-diagonal ≈ **−0.002** (≈0 — the *average* pair is
> unrelated, by construction of centering). Individual pairs depart from 0 — and *that variation*
> is the signal prediction feeds on.

---

## 6.3b 🧸 Pedigree kinship vs. marker reality — and the two glitches markers fix

Before markers, breeders measured relatedness from a **pedigree** (family tree). It's worth seeing
this older method *and why DNA improves on it* — it's the whole reason the paper uses a marker **G**.
Run `code/toy_06_pedigree_kinship.R`.

### The toy pedigree
Founders **P1, P2, C, F**; full sibs **A, B** (= P1×P2); and **D** (= A×C), **E** (= B×F):

```mermaid
flowchart TD
  P1((P1)) --> A((A)); P2((P2)) --> A
  P1 --> B((B)); P2 --> B
  A --> D((D)); C((C)) --> D
  B --> E((E)); F((F)) --> E
  X{{"hidden ancestor X<br/>(not in the records)"}} -.-> P1
  X -.-> C
```

### Pedigree expectation: the relationship matrix A (by hand)
The **numerator relationship matrix A** is built by one rule applied top-down (the *tabular
method*): for individual $i$ with parents $s,d$,

```math
A_{ij}=\tfrac12 (A_{js}+A_{jd}), \qquad A_{ii}=1+\tfrac12 A_{sd}
```

Unknown parents ⇒ treat as unrelated founders ($A=0$). This gives **expected** relationships:

| pair | meaning | pedigree $A$ | kinship $A/2$ |
|------|---------|--------------|----------------|
| A, B | full sibs | **0.500** | 0.250 |
| D, E | first cousins | **0.125** | 0.062 |
| P1, C | "unrelated" founders | **0.000** | 0.000 |

⚠️ Notice **A gives one fixed number per relationship type** — *every* full-sib pair is exactly
0.50, by assumption.

### Marker reality: realized G — and where it disagrees
Now we *gene-drop* DNA through the pedigree and compute the marker **G** (VanRaden, §6.2, centered on
founder allele frequencies). G measures what each pair **actually** inherited:

![pedigree A vs marker G](../figures/22_toy_pedigree_kinship.png)

🔬 Two disagreements pop out (points off the dashed line):

1. **Mendelian sampling — "every full sib is 0.50" is only an *average*.** The pedigree says
   A,B = 0.50; the markers say **G = 0.55**. Full sibs each draw a *random half* of each parent's
   genome, so realized sharing scatters around 0.50 (some sib pairs 0.4, some 0.6). Pedigree A
   can't see this; **G can**. Selecting on realized relationships is more accurate.
2. **Cryptic relatedness — markers catch what the pedigree missed.** The records call P1 and C
   *unrelated founders* (A = 0.00). But they secretly share an unrecorded ancestor **X** (they're
   really half-sibs). The markers expose it: **G = 0.28**. Incomplete pedigrees are the norm in
   real breeding programs — DNA doesn't rely on paperwork.

### Why this is the bridge to the paper
🌱 The black bean panel comes from **two programs (MSU, USDA-ARS)** with **different, often
incomplete pedigrees** — there's no single reliable family tree linking all 415 lines. And even
where pedigrees exist, **realized** relatedness (what GP actually needs, Lesson 7) beats pedigree
*expectations*. So the authors compute relatedness straight from **2,315 SNPs** — exactly the
marker **G** of §6.2. The toy is the study in miniature: *replace the family tree with DNA, gain
both accuracy (Mendelian sampling) and coverage (cryptic relatedness).*

> 🔧 **In practice (R).** Pedigree **A**: `AGHmatrix::Amatrix()`, `nadiv::makeA()`, or `pedigreemm`.
> Marker **G**: `rrBLUP::A.mat()`, `AGHmatrix::Gmatrix()`, or `sommer::A.mat()` — all implement the
> VanRaden formula we wrote by hand in §6.2. (We built both from scratch here so the machinery is
> visible.)

---

## 6.4 Reading the real G: population structure

A heatmap of G (paper's Fig. 1a) shows **blocks** along the diagonal: lines from the same
**breeding program × cycle** light up as related clusters, while across-program pairs stay near
zero. We reproduced this two ways:

🔬 **Kinship by group** (`figures/05_kinship_dist.png`): relationships are clearly **higher
*within* a breeding cycle than *between* cycles** — because each cycle is selected from a limited
set of shared parents.

🔬 **PCA of G** (`figures/04_pca_structure.png`): the first two eigenvectors of G separate the
lines into groups. (PCA of G = principal components of genetic relatedness; the leading PCs *are*
the major axes of population structure.)

🌱 **Why a breeder cares about these blocks — two consequences that drive the whole study:**
1. **Prediction works best within a block.** If your training set is closely related to your
   target lines (same cycle/program), G has strong positive entries linking them → high accuracy.
   Predict *across* blocks (cycle 1 → cycle 2) and the linking entries shrink → accuracy drops.
   **This is the entire story of Lesson 14** (and why the paper's across-cycle accuracies are
   lower until new-cycle lines are added to training).
2. **Structure confounds GWAS.** If a whole subpopulation happens to both carry an allele *and*
   have high yield (for unrelated reasons), a naive test flags a false association. That's why
   GWAS (Lesson 9) must *correct* for structure — using the very PCs we just computed from G.

---

## 6.5 Two roles G plays (don't confuse them)

| Role | Where | What G does |
|------|-------|-------------|
| **Covariance of breeding values** | GBLUP (Lesson 7) | $\mathbf{g} \sim N(\mathbf 0,  \mathbf G \sigma_g^2)$ — relatives have correlated genetic values |
| **Diagnostic of structure** | EDA / GWAS (Lesson 9) | its PCs reveal subpopulations to correct for |

The first role is the engine of prediction: by declaring "breeding values are *correlated
according to G*", the model can infer an unseen line's value from its relatives' phenotypes.

---

## 6.6 Why G beats counting shared genes one-by-one

You might ask: why not just find the trait genes and check who shares them? Because (a) we
mostly don't *know* the genes (Lesson 5's pivot), and (b) **G aggregates tiny signals from all
2,315 SNPs at once.** Complex traits are built from many small effects; G captures their *net*
relatedness without needing to identify any single one. This is precisely why genomic prediction
outperforms marker-by-marker (QTL) selection for polygenic traits like yield.

---

## 6.7 What you should now be able to say
- **G** is the 415×415 genomic relationship matrix; $G_{ij}$ = average product of standardized
  genotypes = realized relatedness.
- Build it by **center+scale markers → $\mathbf{Z}$**, then $\mathbf{G}=\mathbf{Z}\mathbf{Z}^\top/p$.
- Its **diagonal ≈ 1**, **mean off-diagonal ≈ 0**; **blocks** reveal program/cycle structure
  (within-cycle relatedness > between-cycle).
- G plays two roles: **covariance of breeding values** (powers GBLUP) and **structure
  diagnostic** (powers GWAS correction). The block structure foreshadows *within-cycle accuracy >
  across-cycle accuracy*.

👉 Next: **[Lesson 7 — GBLUP](07_GBLUP.md)** — we feed G into a mixed model and actually predict
(and reproduce accuracy 0.64).
