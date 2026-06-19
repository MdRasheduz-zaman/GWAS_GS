# Glossary — Quick Reference

Plain-language definitions of every key term in the course. Each links to the lesson where it's
developed. Skim before an exam; look up while reading.

### Biology & breeding
- **Common bean / *Phaseolus vulgaris*** — the crop; diploid, 2n=22 (**11 chromosomes**),
  self-pollinating. (L1)
- **Market class** — beans grouped by seed type (navy, pinto, **black**, …). (L1)
- **Line** — a near-homozygous, reproducible genotype (because beans self-pollinate). (L1)
- **Breeding cycle** — one round of cross → fix lines → test → select parents. Here: **cycle 1 =
  272 lines** (2018–19), **cycle 2 = 143 lines** (2019). (L1)
- **Check cultivar** — a known variety (Eclipse/Zorro/Zenith) grown repeatedly to map field
  variation. (L3)
- **Genetic gain ($\Delta G$)** — improvement per cycle; $\Delta G = i\,r\,\sigma_A$. (L5)

### Traits
- **YD** seed yield; **SW** 100-seed weight; **DF** days to flowering; **DPM** days to maturity. (L1)
- **App** canning **appearance** (trained-panel score 1–5); **Col** color/darkness; **Text**
  texture (machine peak-force); **L\*a\*b\*** objective color coordinates. (L1)
- **Canning quality** — how good beans are after the standardized canning protocol. (L1)

### Statistics & quantitative genetics
- **Phenotype** — measured trait value = mean + breeding value + environment. (L5)
- **Additive effect ($\alpha_j$)** — the per-copy value of an allele at marker $j$. (L5)
- **Breeding value ($g_i$)** — sum of additive effects, $g_i=\sum_j \alpha_j x_{ij}$; what we
  predict/select on. (L5)
- **Heritability ($h^2$, $H^2$, $h^2_g$)** — fraction of trait variance that's genetic; **caps
  prediction accuracy**. (L3, L5)
- **Fixed effect** — estimated freely, no shrinkage ("I'm sure"). **Random effect** — shrunk
  toward the mean ("probably small"). (L3, L10)
- **BLUP** — Best Linear Unbiased Prediction; the shrunken genotype estimate = one clean number
  per line. (L3)
- **Mixed model** — a model with both fixed and random effects. (L3, L7)
- **Shrinkage / regularization** — pulling estimates toward 0 to avoid overfitting when
  $p\gg n$. (L7, L11)
- **Cross-validation** — hide data, predict it, score honestly; **70/30 ×100** here, with inner
  **10-fold** for tuning. (L13)
- **Prediction accuracy ($r$)** — Pearson correlation between predicted and observed test values.
  (L13)

### Genomics
- **SNP** — single-DNA-position difference; coded **0/1/2** (allele dosage). (L2, L4)
- **GBS** — Genotyping-by-Sequencing; cheap way to find thousands of shared SNPs. (L4)
- **MAF** — minor allele frequency; rare SNPs (<0.01) filtered out. (L4)
- **LD (linkage disequilibrium)** — correlation between nearby SNPs; basis of marker tagging; LD
  pruning removes redundant SNPs. (L4)
- **Imputation** — filling missing genotypes from genetic neighbors (Beagle). (L4)
- **Population structure** — subpopulations (program/cycle) that confound GWAS; corrected with
  PCs. (L6, L9)

### Models in the study
- **G (genomic relationship matrix)** — realized relatedness, $\mathbf G=\mathbf Z\mathbf Z^\top/p$.
  (L6)
- **GBLUP** — genomic prediction via G in a mixed model; **= ridge regression on markers**. (L7)
- **RKHS / Gaussian kernel** — GBLUP with a flexible similarity kernel $\exp(-\theta d^2)$;
  captures non-additive effects. (L8)
- **Bandwidth ($\theta$, `h`)** — kernel decay rate; small = smooth/global, large = local. (L8)
- **KA (kernel averaging)** — mix of bandwidths (0.02, 1, 5) weighted by data; the chosen GP
  model. (L8)
- **GWAS** — tests each SNP for trait association; p-values, Manhattan plot. (L9)
- **FarmCPU** — the iterative GWAS method (GAPIT3) used. (L9)
- **Bonferroni** — strict significance threshold $0.05/\text{#SNPs}=2.16\times10^{-5}$. (L9)
- **GWAS-assisted GP** — significant SNPs added as fixed effects; **hurt accuracy** here. (L10)
- **NIRS** — near-infrared spectroscopy; 551-wavelength chemical fingerprint. (L11)
- **RSI (regularized selection index)** — penalized-regression weighting of spectra into one
  index. (L11)
- **SNV** — Standard Normal Variate, the NIRS preprocessing used. (L11)
- **ST vs MT** — single-trait vs multi-trait prediction. (L12)
- **Secondary trait** — a cheap, correlated trait that boosts a hard target in MT models (SW for
  yield, Text for App). (L12)
- **Genetic covariance ($\sigma_{g_{12}}$)** — the channel through which a secondary trait helps;
  0 → MT = ST. (L12)

### Model-code letters (Table 1)
- **G** = GWAS hits; **C** = correlated trait; **R** = RSI (NIRS). E.g. **MT-GC** = multi-trait +
  GWAS + correlated trait. (L13)
