#!/usr/bin/env Rscript
# Toy for Lesson 13: hide a test set, predict it, correlate -> accuracy; repeat -> distribution.
suppressMessages({library(ggplot2); library(patchwork)})
set.seed(4)
n <- 20
g <- round(rnorm(n, 0, 5), 1)               # true genetic values of 20 lines
names(g) <- paste0("L",1:n)
pred <- g*0.8 + rnorm(n, 0, 3)              # a model's predictions (correlated w/ truth, not perfect)

# one split: hide 6 (30%) as test
tst <- c(2,5,8,11,14,17); trn <- setdiff(1:n,tst)
r1 <- cor(pred[tst], g[tst])
cat("One 70/30 split -> test accuracy r =", round(r1,2),"\n")

# Panel A: who is train vs test
dfA <- data.frame(line=factor(names(g),levels=names(g)), g=g,
                  set=ifelse(1:n %in% tst,"TEST (hidden)","train"))
pA <- ggplot(dfA,aes(line,g,fill=set))+geom_col()+
  scale_fill_manual(values=c("train"="grey75","TEST (hidden)"="#d73027"),name=NULL)+
  labs(title="Step A: hide 30% of lines (their phenotype is concealed from the model)",
       x=NULL,y="true value")+theme_bw(base_size=11)+theme(legend.position="top")

# Panel B: predicted vs observed on the hidden test lines
dfB <- data.frame(observed=g[tst], predicted=pred[tst], line=names(g)[tst])
pB <- ggplot(dfB,aes(observed,predicted,label=line))+
  geom_smooth(method="lm",se=FALSE,color="grey60",linetype="dashed")+
  geom_point(size=3,color="#1b9e77")+ggrepel::geom_text_repel(size=4)+
  labs(title=sprintf("Step B: accuracy = cor(predicted, observed) on TEST = %.2f",r1),
       x="observed (truth)",y="predicted")+theme_bw(base_size=11)

# Panel C: repeat the split many times -> a DISTRIBUTION of accuracies (why the paper uses boxplots)
accs <- replicate(100,{t<-sample(1:n,6); cor(pred[t],g[t])})
pC <- ggplot(data.frame(r=accs),aes("",r))+geom_boxplot(fill="#9ecae1",width=.4)+
  geom_jitter(width=.08,alpha=.3,size=.8)+coord_flip()+
  labs(title="Step C: repeat 100 random splits -> spread of accuracy (one box = one model)",
       x=NULL,y="accuracy r across splits")+theme_bw(base_size=11)
ggsave("../figures/18_toy_cv.png",(pA/pB/pC)+plot_layout(heights=c(1,1.2,.7)),width=8,height=8.5,dpi=120)
cat("mean accuracy over 100 splits =", round(mean(accs),2)," (range",
    round(min(accs),2),"to",round(max(accs),2),")\n")
cat("Wrote figures/18_toy_cv.png\n")
