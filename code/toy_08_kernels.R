#!/usr/bin/env Rscript
# Toy for Lesson 8: genetic distance -> Gaussian kernel similarity, bandwidth effect.
suppressMessages({library(ggplot2); library(reshape2); library(patchwork)})
M <- matrix(c(2,2,0,2,0, 2,0,0,2,0, 0,0,2,0,2, 0,2,2,1,2),nrow=4,byrow=TRUE,
            dimnames=list(paste0("L",1:4),paste0("S",1:5)))
D <- as.matrix(dist(M))^2; D <- D/mean(D)            # scaled squared distance
hs <- c(0.02,1,5)
cat("scaled squared-distance matrix D:\n"); print(round(D,2))

# Panel A: K = exp(-h*d^2) vs distance, for 3 bandwidths
dd <- seq(0,max(D)*1.1,length=100)
curv <- do.call(rbind, lapply(hs,function(h) data.frame(d=dd, K=exp(-h*dd), h=paste0("theta=",h))))
pA <- ggplot(curv,aes(d,K,color=h))+geom_line(linewidth=1)+
  scale_color_manual(values=c("theta=0.02"="#1b9e77","theta=1"="#7570b3","theta=5"="#d95f02"),name="bandwidth")+
  labs(title="Step A: similarity decays with distance",
       subtitle="small theta = smooth/global ; large theta = local/peaky",
       x=expression("scaled genetic distance  "*d^2), y=expression(K==e^{-theta*d^2}))+
  theme_bw(base_size=12)+theme(legend.position="top")

# Panel B: the 4x4 kernel matrices at each bandwidth
mats <- do.call(rbind, lapply(hs,function(h){m<-melt(exp(-h*D));m$h<-paste0("theta=",h);m}))
mats$Var1<-factor(mats$Var1,levels=rev(rownames(D)))
pB <- ggplot(mats,aes(Var2,Var1,fill=value))+geom_tile(color="white")+
  geom_text(aes(label=sprintf("%.2f",value)),size=3)+facet_wrap(~h)+
  scale_fill_gradient(low="white",high="#2166AC",limits=c(0,1),name="K")+
  labs(title="Step B: same 4 lines, three bandwidths",x=NULL,y=NULL)+
  theme_minimal(base_size=11)
ggsave("../figures/15_toy_kernels.png",(pA/pB)+plot_layout(heights=c(1,1.1)),width=8,height=8,dpi=120)
cat("Wrote figures/15_toy_kernels.png\n")
