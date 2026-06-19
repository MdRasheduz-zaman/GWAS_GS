#!/usr/bin/env Rscript
# ============================================================================
# GWAS reproduction (teaching version of FarmCPU's core idea):
#   single-marker test of each SNP, CORRECTED for population structure (3 PCs),
#   Bonferroni threshold, Manhattan + QQ plot.  Compares yield (low h2) vs
#   color (high h2) to show why h2 drives discovery.
# ============================================================================
suppressMessages({library(ggplot2)})
load("../repo/GB_BLB.RData")
dir.create("../figures", showWarnings = FALSE)

M   <- as.matrix(GB_BLB$geno)
pos <- GB_BLB$SNP_Position
Y   <- GB_BLB$pheno

# population-structure covariates: first 3 PCs of the marker matrix (as paper does)
Z   <- scale(M)
pca <- prcomp(Z, center = FALSE)
PCs <- pca$x[, 1:3]

run_gwas <- function(trait) {
  y   <- Y[[trait]]
  ok  <- which(!is.na(y))
  yk  <- y[ok]; Mk <- M[ok, ]; Pk <- PCs[ok, ]
  pval <- numeric(ncol(Mk))
  for (j in seq_len(ncol(Mk))) {
    snp <- Mk[, j]
    if (sd(snp) == 0) { pval[j] <- NA; next }
    fit <- lm(yk ~ Pk + snp)            # y ~ PC1+PC2+PC3 + SNP  (structure-corrected)
    cf  <- summary(fit)$coefficients
    pval[j] <- if ("snp" %in% rownames(cf)) cf["snp", 4] else NA
  }
  data.frame(SNP = colnames(Mk), Chr = pos$Chromosome, Pos = pos$Position,
             P = pval, trait = trait)
}

bonf <- 0.05 / ncol(M)
cat("Bonferroni threshold = 0.05 /", ncol(M), "=", signif(bonf, 3),
    " (-log10 =", round(-log10(bonf), 2), ")\n\n")

res <- rbind(run_gwas("yd18"), run_gwas("col18"))
for (tr in c("yd18", "col18")) {
  n_hit <- sum(res$P[res$trait == tr] < bonf, na.rm = TRUE)
  cat(sprintf("%-6s : %d SNP(s) pass Bonferroni\n", tr, n_hit))
}

# Manhattan plot
res$Chr <- factor(res$Chr)
res$logp <- -log10(res$P)
# cumulative position for plotting
res <- res[order(res$Chr, res$Pos), ]
res$idx <- ave(seq_len(nrow(res)), res$trait, FUN = seq_along)
ggplot(res, aes(idx, logp, color = Chr)) +
  geom_point(size = .8, alpha = .8) +
  geom_hline(yintercept = -log10(bonf), linetype = "dashed", color = "red") +
  facet_wrap(~trait, ncol = 1, scales = "free_x",
             labeller = as_labeller(c(yd18 = "Yield 2018 (lower heritability)",
                                      col18 = "Color 2018 (higher heritability)"))) +
  scale_color_manual(values = rep(c("#2166AC", "#92C5DE"), 6)) +
  labs(title = "GWAS Manhattan plots (PC-corrected single-marker scan)",
       subtitle = "Dashed red = Bonferroni threshold", x = "SNP (ordered by chromosome)",
       y = expression(-log[10](p))) +
  theme_bw(base_size = 11) + theme(legend.position = "none")
ggsave("../figures/06_manhattan.png", width = 9, height = 6, dpi = 120)
cat("\nManhattan -> figures/06_manhattan.png\n")
