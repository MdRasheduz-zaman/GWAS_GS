#!/usr/bin/env Rscript
# Toy for Lesson 11: a selection index from many predictors (p ~ n) needs REGULARIZATION.
# 12 individuals, 10 "wavelengths"; OLS overfits (perfect on train, junk on test);
# a ridge penalty shrinks weights and generalizes. Shows index = X %*% beta.
suppressMessages({library(ggplot2); library(reshape2); library(patchwork)})
set.seed(21)
n <- 26; p <- 14                       # p (14 wavelengths) > training size (12): OLS is doomed
X <- matrix(rnorm(n*p), n, p, dimnames=list(NULL, paste0("nm",1:p)))
# truth: only 2 wavelengths matter; the rest is noise
b_true <- c(2.0, 0,0, -1.5, rep(0,10))
y <- as.vector(X %*% b_true) + rnorm(n, 0, 4)
trn <- 1:12; tst <- 13:26              # 12 train, 14 test  (p > n_train -> overfit territory)

ridge <- function(X,y,lam) solve(crossprod(X)+lam*diag(ncol(X)), crossprod(X,y))
bO <- ridge(X[trn,],y[trn],1e-6)     # ~OLS (no penalty) — underdetermined, weights explode
bR <- ridge(X[trn,],y[trn],6)        # penalized
acc <- function(b) c(train=cor(X[trn,]%*%b,y[trn]), test=cor(X[tst,]%*%b,y[tst]))
cat("OLS   accuracy:", sprintf("train=%.2f test=%.2f",acc(bO)[1],acc(bO)[2]),"\n")
cat("Ridge accuracy:", sprintf("train=%.2f test=%.2f",acc(bR)[1],acc(bR)[2]),"\n")

# Panel A: weights OLS vs ridge vs truth
wd <- rbind(data.frame(nm=paste0("nm",1:p),w=as.vector(bO),k="OLS (no penalty)"),
            data.frame(nm=paste0("nm",1:p),w=as.vector(bR),k="Ridge (penalized)"),
            data.frame(nm=paste0("nm",1:p),w=b_true,k="TRUE weights"))
wd$nm<-factor(wd$nm,levels=paste0("nm",1:p))
pA <- ggplot(wd,aes(nm,w,fill=k))+geom_col(position="dodge")+
  scale_fill_manual(values=c("OLS (no penalty)"="#d95f02","Ridge (penalized)"="#1b9e77","TRUE weights"="grey50"),name=NULL)+
  labs(title="Step A: index weights beta (one per wavelength)",
       subtitle="OLS weights explode; ridge shrinks toward the truth",x=NULL,y=expression(beta))+
  theme_bw(base_size=11)+theme(legend.position="top",axis.text.x=element_text(angle=45,hjust=1))

# Panel B: index vs trait on TEST set
dfB <- rbind(data.frame(index=as.vector(X[tst,]%*%bO),y=y[tst],k=sprintf("OLS test r=%.2f",acc(bO)[2])),
             data.frame(index=as.vector(X[tst,]%*%bR),y=y[tst],k=sprintf("Ridge test r=%.2f",acc(bR)[2])))
pB <- ggplot(dfB,aes(index,y))+geom_point(color="#2166AC")+geom_smooth(method="lm",se=FALSE,color="grey50")+
  facet_wrap(~k,scales="free_x")+
  labs(title="Step B: selection index vs trait, on held-out individuals",
       x="selection index  (X %*% beta)",y="trait")+theme_bw(base_size=11)
ggsave("../figures/17_toy_rsi.png",(pA/pB)+plot_layout(heights=c(1,1)),width=8.5,height=7,dpi=120)
cat("Wrote figures/17_toy_rsi.png\n")
