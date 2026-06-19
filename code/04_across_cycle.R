#!/usr/bin/env Rscript
# ============================================================================
# CAPSTONE reproduction of the study's headline findings (Objectives 1 & 4):
#  (1) ACROSS-CYCLE accuracy RISES as cycle-2 lines are added to training.
#  (2) MULTI-TRAIT (yield + seed weight) beats SINGLE-TRAIT across cycles.
# Trait: yield 2019 (yd19). Secondary: seed weight 2019 (sw19).
# ST GBLUP via rrBLUP (fast, closed form); MT via BGLR::Multitrait.
# ============================================================================
suppressMessages({library(rrBLUP); library(BGLR)})
set.seed(7)
load("../repo/GB_BLB.RData")

M <- as.matrix(GB_BLB$geno)
Y <- GB_BLB$pheno
G <- tcrossprod(scale(M)) / ncol(M)

y  <- Y$yd19            # primary trait, all 415 lines
sw <- Y$sw19            # secondary trait (correlated, cheap, heritable)

cyc1 <- 1:272
cyc2 <- 273:415
c1 <- cyc1[!is.na(y[cyc1])]      # cycle-1 lines with yield measured
c2 <- cyc2[!is.na(y[cyc2])]      # cycle-2 lines with yield measured
cat("Cycle 1 measured:", length(c1), " | Cycle 2 measured:", length(c2), "\n")
cat("Genetic corr yield~seedweight (2019):", round(cor(y, sw, use="complete.obs"), 2), "\n\n")

st_gblup <- function(train, test) {        # closed-form GBLUP prediction
  yNA <- y; yNA[setdiff(seq_along(y), train)] <- NA
  fit <- mixed.solve(y = yNA[c(train,test)], K = G[c(train,test), c(train,test)])
  pred <- as.numeric(fit$beta) + fit$u
  names(pred) <- c(train, test)
  pred[as.character(test)]
}

mt_gblup <- function(train, test) {        # multi-trait: yield + seed weight
  idx <- c(train, test)
  Ymt <- cbind(y[idx], sw[idx])            # secondary KNOWN on test lines
  Ymt[match(test, idx), 1] <- NA           # hide ONLY primary on test
  sink(tempfile())
  fm <- Multitrait(y = Ymt, ETA = list(list(K = G[idx, idx], model = "RKHS")),
                   resCov = list(type = "DIAG"), nIter = 3000, burnIn = 800)
  sink()
  pred <- fm$ETAHat[, 1]; names(pred) <- idx
  pred[as.character(test)]
}

# Lighter settings so this finishes in a few minutes on a laptop (the paper used
# 100 reps and 12,000 MCMC iterations). Raise nrep/nIter for publication-grade stability.
props <- c(0, 0.10, 0.20, 0.30, 0.40)
nrep  <- 12
out <- data.frame()
for (p in props) {
  accST <- accMT <- numeric(nrep)
  for (r in 1:nrep) {
    nadd <- round(p * length(c2))
    add  <- if (nadd > 0) sample(c2, nadd) else integer(0)
    test <- setdiff(c2, add)
    train <- c(c1, add)
    accST[r] <- cor(st_gblup(train, test), y[test])
    accMT[r] <- cor(mt_gblup(train, test), y[test])
  }
  out <- rbind(out,
    data.frame(prop = p*100, model = "ST (yield only)",      acc = mean(accST)),
    data.frame(prop = p*100, model = "MT (yield + seed wt)", acc = mean(accMT)))
  cat(sprintf("prop=%2.0f%%  ST=%.3f  MT=%.3f  (MT gain %+.0f%%)\n",
              p*100, mean(accST), mean(accMT),
              100*(mean(accMT)-mean(accST))/abs(mean(accST))))
}

suppressMessages(library(ggplot2))
ggplot(out, aes(prop, acc, color = model)) +
  geom_line(linewidth = 1) + geom_point(size = 2) +
  scale_color_manual(values = c("ST (yield only)" = "#7570B3",
                                "MT (yield + seed wt)" = "#1B9E77")) +
  labs(title = "Predicting the NEW breeding cycle (yield 2019)",
       subtitle = "Accuracy rises as cycle-2 lines join training; MT (correlated trait) leads",
       x = "% of cycle-2 lines added to the training set", y = "Prediction accuracy (r)",
       color = "") +
  theme_bw(base_size = 12) + theme(legend.position = "top")
ggsave("../figures/07_across_cycle.png", width = 8, height = 6, dpi = 120)
saveRDS(out, "../figures/across_cycle.rds")
cat("\nDONE -> figures/07_across_cycle.png\n")
