#!/usr/bin/env Rscript
# Toy for Lesson 9: test each SNP for trait association; p-values; Bonferroni; Manhattan.
suppressMessages({library(ggplot2); library(patchwork)})
set.seed(7)
n <- 80; nsnp <- 12
M <- matrix(sample(0:2, n*nsnp, TRUE), n, nsnp, dimnames=list(NULL,paste0("SNP",1:nsnp)))
causal <- 8                                   # only SNP8 truly affects the trait
y <- 10 + 3*M[,causal] + rnorm(n, 0, 3)        # trait = baseline + real effect + noise

pval <- sapply(1:nsnp, function(j) summary(lm(y ~ M[,j]))$coefficients[2,4])
bonf <- 0.05/nsnp
cat("Bonferroni threshold = 0.05/12 =", signif(bonf,3), " -log10 =", round(-log10(bonf),2),"\n")
cat("SNP8 (causal) p =", signif(pval[causal],3), " -> ",
    ifelse(pval[causal]<bonf,"PASSES","fails"),"\n")
cat("smallest p among null SNPs =", signif(min(pval[-causal]),3),"\n")

# Panel A: trait by genotype, causal vs a null SNP
dfA <- rbind(
  data.frame(geno=factor(M[,causal]), y=y, snp="SNP8 (causal): means differ -> association"),
  data.frame(geno=factor(M[,1]),      y=y, snp="SNP1 (null): means overlap -> no association"))
pA <- ggplot(dfA,aes(geno,y,fill=geno))+geom_boxplot(alpha=.7)+facet_wrap(~snp)+
  scale_fill_manual(values=c("0"="#deebf7","1"="#9ecae1","2"="#3182bd"),guide="none")+
  labs(title="Step A: a SNP is 'associated' if genotype groups have different trait means",
       x="genotype (copies of alt allele)",y="trait")+theme_bw(base_size=11)

# Panel B: mini Manhattan
dfB <- data.frame(snp=factor(paste0("SNP",1:nsnp),levels=paste0("SNP",1:nsnp)),
                  logp=-log10(pval), causal=(1:nsnp==causal))
pB <- ggplot(dfB,aes(snp,logp,fill=causal))+geom_col()+
  geom_hline(yintercept=-log10(bonf),linetype="dashed",color="red")+
  annotate("text",x=2,y=-log10(bonf)+.2,label="Bonferroni",color="red",size=3)+
  scale_fill_manual(values=c("FALSE"="grey70","TRUE"="#d73027"),guide="none")+
  labs(title="Step B: mini Manhattan plot (only the causal SNP clears the line)",
       x=NULL,y=expression(-log[10](p)))+theme_bw(base_size=11)+
  theme(axis.text.x=element_text(angle=45,hjust=1))
ggsave("../figures/16_toy_gwas.png",(pA/pB)+plot_layout(heights=c(1,1)),width=8.5,height=7,dpi=120)
cat("Wrote figures/16_toy_gwas.png\n")
