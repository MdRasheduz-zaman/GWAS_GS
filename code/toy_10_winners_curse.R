#!/usr/bin/env Rscript
# Toy for Lesson 10: WHY forcing GWAS hits in as FIXED (unshrunk) effects hurts.
# Polygenic trait (many tiny effects). We "discover" top SNPs in a training set, then
# show (A) winner's curse: their training effects are inflated vs an independent estimate,
# and (B) a predictor that TRUSTS those few unshrunk effects loses to one that SHRINKS
# all markers (GBLUP-style ridge).
suppressMessages({library(ggplot2); library(patchwork)})
set.seed(19)
n <- 200; p <- 300
M  <- matrix(rbinom(n*p, 2, 0.5), n, p)           # genotypes 0/1/2
Z  <- scale(M)
alpha <- rnorm(p, 0, 0.15)                         # MANY tiny effects (polygenic, none large)
g  <- as.vector(Z %*% alpha)
y  <- g + rnorm(n, 0, sd(g))                       # h2 ~ 0.5
disc <- 1:100; test <- 101:200                     # discovery vs test lines

# --- marginal GWAS-style scan on the DISCOVERY set; pick top-k "hits" ---
est_disc <- sapply(1:p, function(j) coef(lm(y[disc] ~ Z[disc,j]))[2])
est_test <- sapply(1:p, function(j) coef(lm(y[test] ~ Z[test,j]))[2])   # independent re-estimate
k <- 10
hits <- order(abs(est_disc), decreasing=TRUE)[1:k]
cat(sprintf("Top-%d 'hits': mean|effect| in discovery = %.3f  but independent re-estimate = %.3f  (true = %.3f)\n",
            k, mean(abs(est_disc[hits])), mean(abs(est_test[hits])), mean(abs(alpha[hits]))))

# --- predict TEST two ways ---
# (1) FIXED: trust the k hits, estimate them unshrunk by OLS, predict from them
bf <- coef(lm(y[disc] ~ Z[disc,hits]))
predF <- cbind(1, Z[test,hits]) %*% bf
# (2) SHRUNK (GBLUP-style ridge over ALL markers)
ridge <- function(X,y,lam) solve(crossprod(X)+lam*diag(ncol(X)), crossprod(X,y))
bs <- ridge(Z[disc,], y[disc], p)                  # heavy shrinkage (~ equal small effects)
predS <- Z[test,] %*% bs
accF <- cor(predF, y[test]); accS <- cor(predS, y[test])
cat(sprintf("Test accuracy:  FIXED hits = %.2f   |   SHRINK all (GBLUP) = %.2f\n", accF, accS))

# Panel A: winner's curse — discovery effect vs independent re-estimate (hits)
dfA <- data.frame(disc=est_disc, test=est_test, hit=factor(ifelse(1:p %in% hits,"selected 'hit'","other SNP")))
pA <- ggplot(dfA, aes(disc, test, color=hit)) +
  geom_abline(slope=1, linetype="dashed", color="grey60") + geom_hline(yintercept=0, color="grey80") +
  geom_point(aes(size=hit, alpha=hit)) +
  scale_color_manual(values=c("selected 'hit'"="#d73027","other SNP"="grey75"), name=NULL) +
  scale_size_manual(values=c("selected 'hit'"=2.6,"other SNP"=1), guide="none") +
  scale_alpha_manual(values=c("selected 'hit'"=1,"other SNP"=.4), guide="none") +
  labs(title="Step A: winner's curse",
       subtitle="hits picked for big DISCOVERY effects shrink toward 0 when re-estimated independently",
       x="effect estimated in discovery set", y="effect re-estimated in test set") +
  theme_bw(base_size=11) + theme(legend.position="top")

# Panel B: test accuracy fixed vs shrunk
dfB <- data.frame(model=c("FIXED: trust the GWAS hits\n(unshrunk)","SHRINK all markers\n(GBLUP)"),
                  acc=c(accF, accS))
pB <- ggplot(dfB, aes(reorder(model,acc), acc, fill=model)) + geom_col(width=.6) +
  geom_text(aes(label=sprintf("%.2f",acc)), vjust=-0.3, size=4) +
  scale_fill_manual(values=c("FIXED: trust the GWAS hits\n(unshrunk)"="#d95f02","SHRINK all markers\n(GBLUP)"="#1b9e77"), guide="none") +
  scale_y_continuous(expand=expansion(mult=c(0,.18))) +
  labs(title="Step B: which predicts the test set better?",
       subtitle="for a polygenic trait, shrinking everything beats trusting a few unstable hits",
       x=NULL, y="test accuracy (r)") + theme_bw(base_size=11)
ggsave("../figures/20_toy_winners_curse.png", (pA/pB)+plot_layout(heights=c(1.2,1)), width=8, height=7.6, dpi=120)
cat("Wrote figures/20_toy_winners_curse.png\n")
