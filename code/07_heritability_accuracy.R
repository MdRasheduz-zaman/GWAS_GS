#!/usr/bin/env Rscript
# ============================================================================
# Substantiates the course's central claim: HERITABILITY CAPS ACCURACY.
# For every 2018 trait: estimate genomic heritability h2_g (from GBLUP variance
# components) AND measure GBLUP prediction accuracy (70/30, averaged). Plot one
# against the other -> the paper's core relationship, on the real data.
# ============================================================================
suppressMessages({library(rrBLUP); library(ggplot2)})
set.seed(11)
load("../repo/GB_BLB.RData")
M <- as.matrix(GB_BLB$geno); Y <- GB_BLB$pheno
G <- tcrossprod(scale(M)) / ncol(M)

traits <- c(yd18="yield", sw18="seed weight", dm18="days to maturity",
            df18="days to flowering", app18="appearance", text18="texture",
            col18="color", l18="L*", a18="a*", b18="b*")
nrep <- 25
res <- data.frame()
for (tr in names(traits)) {
  y0 <- Y[[tr]]; idx <- which(!is.na(y0))
  y <- y0[idx]; Gs <- G[idx, idx]; n <- length(y)
  # h2_g from full-data variance components
  vc <- mixed.solve(y = y, K = Gs)
  h2 <- vc$Vu / (vc$Vu + vc$Ve)
  # prediction accuracy, 70/30 x nrep
  acc <- numeric(nrep)
  for (r in 1:nrep) {
    tst <- sample(1:n, round(0.3*n)); yNA <- y; yNA[tst] <- NA
    fit <- mixed.solve(y = yNA, K = Gs)
    acc[r] <- cor(y[tst], fit$u[tst])
  }
  res <- rbind(res, data.frame(trait=tr, label=traits[tr], h2=h2, acc=mean(acc)))
  cat(sprintf("%-7s h2_g=%.2f  accuracy=%.2f\n", tr, h2, mean(acc)))
}

cat(sprintf("\ncorrelation(h2_g, accuracy) across traits = %.2f\n", cor(res$h2, res$acc)))
ggplot(res, aes(h2, acc)) +
  geom_smooth(method="lm", se=FALSE, color="grey60", linetype="dashed") +
  geom_point(size=3, color="#1B9E77") +
  ggrepel::geom_text_repel(aes(label=label), size=3.5) +
  labs(title="Heritability caps prediction accuracy (real 2018 data)",
       subtitle=sprintf("Each point = one trait. r(h2, accuracy) = %.2f", cor(res$h2,res$acc)),
       x=expression("genomic heritability  "*h[g]^2), y="GBLUP prediction accuracy (r)") +
  theme_bw(base_size=12)
ggsave("../figures/10_heritability_vs_accuracy.png", width=7.5, height=6, dpi=120)
write.csv(res, "../figures/heritability_accuracy.csv", row.names=FALSE)
cat("Wrote figures/10_heritability_vs_accuracy.png\n")
