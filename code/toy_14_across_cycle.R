#!/usr/bin/env Rscript
# Toy for Lesson 14: predicting a NEW, weakly-related cycle is hard; accuracy RISES as
# you fold new-cycle lines into training. Two genetic clusters (cycle 1 vs cycle 2)
# created via different allele frequencies; GBLUP via rrBLUP.
suppressMessages({library(rrBLUP); library(ggplot2); library(patchwork)})
set.seed(31)
p <- 800; n1 <- 150; n2 <- 150; h2 <- 0.6
# two subpopulations with INDEPENDENT allele frequencies -> near-zero between-cluster relatedness
f1 <- runif(p, .05, .95); f2 <- runif(p, .05, .95)
M1 <- sapply(1:p, function(j) rbinom(n1, 2, f1[j]))
M2 <- sapply(1:p, function(j) rbinom(n2, 2, f2[j]))
M  <- rbind(M1, M2); rownames(M) <- c(paste0("C1_",1:n1), paste0("C2_",1:n2))
cyc <- factor(rep(c("Cycle 1","Cycle 2"), c(n1,n2)))
alpha <- rnorm(p, 0, 0.12)
g <- as.vector(scale(M) %*% alpha)
y <- g + rnorm(n1+n2, 0, sd(g)*sqrt((1-h2)/h2))   # target heritability h2
G <- tcrossprod(scale(M)) / p

c1 <- 1:n1; c2 <- (n1+1):(n1+n2)
gblup_acc <- function(train, test){
  idx <- c(train, test)                       # training first, then test
  yNA <- y[idx]; yNA[(length(train)+1):length(idx)] <- NA   # hide test phenotypes
  fit <- mixed.solve(y=yNA, K=G[idx, idx])
  predtest <- fit$u[(length(train)+1):length(idx)]          # positional, not by name
  cor(predtest, y[test])
}
# reference: within-cycle-1 accuracy (lots of relatives)
ref <- mean(replicate(20,{t<-sample(c1,round(.3*n1)); gblup_acc(setdiff(c1,t),t)}))

props <- c(0,.25,.5,.75); nrep <- 20
acc <- sapply(props, function(pr) mean(replicate(nrep,{
  add <- if(pr>0) sample(c2, round(pr*n2)) else integer(0)
  test <- setdiff(c2, add)
  gblup_acc(c(c1, add), test)
})))
names(M)<-NULL
cat("within-cycle-1 accuracy (reference):", round(ref,2), "\n")
cat("across-cycle accuracy by % of cycle-2 added to training:\n")
print(setNames(round(acc,2), paste0(props*100,"%")))

# Panel A: PCA of G -> two clusters (weak between-cycle relatedness)
ev <- eigen(G, symmetric=TRUE)
pcd <- data.frame(PC1=ev$vectors[,1], PC2=ev$vectors[,2], cyc=cyc)
pA <- ggplot(pcd, aes(PC1,PC2,color=cyc)) + geom_point(alpha=.8,size=2) +
  scale_color_manual(values=c("Cycle 1"="#0072B2","Cycle 2"="#D55E00"), name=NULL) +
  labs(title="Two cycles = two weakly-related families",
       subtitle="PCA of the toy relationship matrix G", x="PC1", y="PC2") +
  theme_bw(base_size=11) + theme(legend.position="top")

# Panel B: accuracy rises as cycle-2 lines are added
dfB <- data.frame(prop=props*100, acc=acc)
pB <- ggplot(dfB, aes(prop, acc)) +
  geom_hline(yintercept=ref, linetype="dashed", color="#0072B2") +
  annotate("text", x=8, y=ref+.03, label="within-cycle-1 accuracy", color="#0072B2", size=3, hjust=0) +
  geom_line(color="#D55E00", linewidth=1) + geom_point(color="#D55E00", size=3) +
  geom_text(aes(label=sprintf("%.2f",acc)), vjust=-0.8, size=3.5) +
  scale_y_continuous(expand=expansion(mult=c(.05,.15))) +
  labs(title="Predicting cycle 2: accuracy climbs as you add cycle-2 lines to training",
       x="% of cycle-2 lines added to training", y="accuracy on remaining cycle-2 (r)") +
  theme_bw(base_size=11)
ggsave("../figures/21_toy_across_cycle.png", (pA/pB)+plot_layout(heights=c(1,1)), width=8, height=7.4, dpi=120)
cat("Wrote figures/21_toy_across_cycle.png\n")
