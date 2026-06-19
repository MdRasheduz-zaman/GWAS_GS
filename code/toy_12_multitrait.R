#!/usr/bin/env Rscript
# Toy for Lesson 12: a correlated SECONDARY trait (known on the test lines) rescues
# prediction of a PRIMARY trait when relatedness can't help (the across-cycle regime).
# We model the hard case: test lines are a NEW family, UNRELATED to training, so a
# relatedness-only single-trait (ST) model can only guess the mean -> accuracy ~0.
# A multi-trait (MT) model leans on each test line's own measured secondary trait.
suppressMessages({library(ggplot2); library(patchwork); library(MASS)})
set.seed(8)

sim_once <- function(rho, n = 60, h2 = 0.5) {
  Sg <- matrix(c(1, rho, rho, 1), 2)              # genetic covariance of the 2 traits
  G  <- mvrnorm(n, c(0, 0), Sg)                   # true breeding values g1, g2
  Ve <- (1 - h2) / h2                             # residual var to hit target h2
  y1 <- G[,1] + rnorm(n, 0, sqrt(Ve))             # primary phenotype
  y2 <- G[,2] + rnorm(n, 0, sqrt(Ve))             # secondary phenotype
  trn <- 1:40; tst <- 41:n                        # test = a NEW, unrelated family
  # ST: no relatives in training -> best guess is the training mean -> ~0 correlation
  predST <- rep(mean(y1[trn]), length(tst))
  # MT: learn how y2 tracks the primary on training, apply to test (y2 is KNOWN on test!)
  b <- coef(lm(y1[trn] ~ y2[trn]))                # primary ~ secondary
  predMT <- b[1] + b[2] * y2[tst]
  c(ST = suppressWarnings(cor(predST, G[tst,1])),
    MT = cor(predMT, G[tst,1]))
}

# accuracy vs genetic correlation rho (average several reps)
rhos <- seq(0, 0.9, 0.15)
acc <- t(sapply(rhos, function(r) rowMeans(replicate(40, sim_once(r)))))
df <- data.frame(rho = rep(rhos, 2),
                 acc = c(acc[,"ST"], acc[,"MT"]),
                 model = rep(c("ST (primary only)", "MT (+ correlated trait)"), each = length(rhos)))
df$acc[is.na(df$acc)] <- 0
cat("accuracy vs genetic correlation:\n"); print(round(acc, 2))

pA <- ggplot(df, aes(rho, acc, color = model)) +
  geom_line(linewidth = 1) + geom_point(size = 2.5) +
  scale_color_manual(values = c("ST (primary only)" = "#7570b3", "MT (+ correlated trait)" = "#1b9e77"), name = NULL) +
  labs(title = "MT helps more as the two traits become more correlated",
       subtitle = "across-cycle regime: test lines unrelated to training, so ST is stuck at the mean",
       x = expression("genetic correlation between traits  "*rho), y = "accuracy on hidden primary (r)") +
  theme_bw(base_size = 11) + theme(legend.position = "top")

# one detailed scatter at rho = 0.8
set.seed(8); Sg <- matrix(c(1,.8,.8,1),2); n<-60; h2<-.5; Ve<-(1-h2)/h2
G <- mvrnorm(n,c(0,0),Sg); y1<-G[,1]+rnorm(n,0,sqrt(Ve)); y2<-G[,2]+rnorm(n,0,sqrt(Ve))
trn<-1:40; tst<-41:n; b<-coef(lm(y1[trn]~y2[trn]))
sc <- rbind(
  data.frame(obs=G[tst,1], pred=mean(y1[trn]),       m="ST: predict the mean (no skill, r~0)"),
  data.frame(obs=G[tst,1], pred=b[1]+b[2]*y2[tst],   m=sprintf("MT: use each line's secondary (r=%.2f)", cor(b[1]+b[2]*y2[tst],G[tst,1]))))
pB <- ggplot(sc, aes(obs, pred)) + geom_point(color="#2166AC", alpha=.8) +
  geom_smooth(method="lm", se=FALSE, color="grey60", linetype="dashed") + facet_wrap(~m) +
  labs(title=expression("Hidden-primary prediction at "*rho*" = 0.8"),
       x="true primary breeding value", y="predicted") + theme_bw(base_size=11)

ggsave("../figures/19_toy_multitrait.png", (pA/pB)+plot_layout(heights=c(1,1)), width=8.5, height=7.5, dpi=120)
cat("Wrote figures/19_toy_multitrait.png\n")
