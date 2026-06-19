#!/usr/bin/env Rscript
# ============================================================================
# Lesson 3 demo: SpATS spatial model recovers clean genotype BLUPs from noisy
# field plots. The study's RAW plot data (BBL_phenotype_2020.xlsx) was never
# released, so we SIMULATE a realistic field with KNOWN truth and show SpATS
# (a) removes a spatial gradient and (b) recovers the true genotype values.
# ============================================================================
suppressMessages(library(SpATS)); suppressMessages(library(ggplot2))
set.seed(42)

nG <- 150; reps <- 2; nrow_f <- 20; ncol_f <- 15   # 300 plots on a 20x15 grid
true_g <- setNames(rnorm(nG, 0, 1), paste0("G", 1:nG))   # TRUE genotype effects

plots <- expand.grid(R = 1:nrow_f, P = 1:ncol_f)[1:(nG*reps), ]
plots$ID <- factor(rep(names(true_g), reps))
# smooth spatial gradient (fertility/moisture) — a diagonal bump
plots$spatial <- 2.5 * sin(plots$R/4) + 1.8 * cos(plots$P/3) + 0.05*plots$R
plots$y <- 10 + true_g[as.character(plots$ID)] + plots$spatial + rnorm(nrow(plots), 0, 1)
plots$row_f <- factor(plots$R); plots$col_f <- factor(plots$P)

cat("raw plot value vs TRUE genotype effect: r =",
    round(cor(plots$y, true_g[as.character(plots$ID)]), 2),
    " (field noise hides the genetics)\n")

fit <- SpATS(response = "y",
             spatial  = ~ PSANOVA(P, R, nseg = c(14, 18)),
             genotype = "ID", genotype.as.random = TRUE,
             random   = ~ row_f + col_f, data = plots,
             control  = list(tolerance = 1e-03, monitoring = 0))

blup <- predict(fit, which = "ID")
b <- setNames(blup$predicted.values, as.character(blup$ID))
common <- names(true_g)
cat("SpATS BLUP vs TRUE genotype effect : r =", round(cor(b[common], true_g[common]), 2),
    " (spatial noise removed -> genetics recovered)\n")
cat("generalized heritability (SpATS)   :", round(getHeritability(fit), 2), "\n")

# figure: fitted spatial surface (what SpATS subtracts) + truth-vs-BLUP
png("../figures/09_spats_demo.png", width = 1000, height = 430, res = 110)
par(mfrow = c(1, 2))
sp <- tapply(fitted(fit) - (mean(plots$y)), list(plots$R, plots$P), mean)
image(1:nrow_f, 1:ncol_f, sp, col = hcl.colors(20, "RdBu", rev = TRUE),
      xlab = "field row", ylab = "field column",
      main = "Field gradient SpATS estimates\n(and removes)")
plot(true_g[common], b[common], pch = 19, col = "#1B9E77",
     xlab = "TRUE genotype effect", ylab = "SpATS BLUP",
     main = sprintf("BLUP recovers truth (r = %.2f)", cor(b[common], true_g[common])))
abline(lm(b[common] ~ true_g[common]), col = "red", lwd = 2)
dev.off()
cat("Wrote figures/09_spats_demo.png\n")
