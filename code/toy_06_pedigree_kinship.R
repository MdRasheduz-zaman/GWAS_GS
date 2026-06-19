#!/usr/bin/env Rscript
# Toy for Lesson 6: pedigree kinship (A matrix, by the tabular method) vs marker-based
# realized relationship (G). Markers (a) reveal Mendelian-sampling variation around the
# pedigree EXPECTATION (full sibs are not all exactly 0.5) and (b) catch CRYPTIC relatedness
# that an incomplete pedigree misses (two "unrelated founders" that secretly share an ancestor).
suppressMessages({library(ggplot2); library(reshape2); library(patchwork)})
set.seed(42)

## ---- (1) RECORDED pedigree -> expected relationship matrix A (tabular method) ----
# parents listed before offspring; 0 = unknown parent (treated as an unrelated founder)
ped <- data.frame(
  id   = c("P1","P2","C","F","A","B","D","E"),
  sire = c( 0,   0,   0,  0, "P1","P1","A","B"),
  dam  = c( 0,   0,   0,  0, "P2","P2","C","F"))
idx <- setNames(seq_len(nrow(ped)), ped$id)
n <- nrow(ped); A <- matrix(0, n, n, dimnames=list(ped$id, ped$id))
for (i in 1:n){
  s <- ped$sire[i]; d <- ped$dam[i]
  si <- if (s!="0") idx[[s]] else 0; di <- if (d!="0") idx[[d]] else 0
  A[i,i] <- 1 + if (si>0 && di>0) 0.5*A[si,di] else 0
  for (j in seq_len(i-1)){
    asj <- if (si>0) A[j,si] else 0; adj <- if (di>0) A[j,di] else 0
    A[i,j] <- A[j,i] <- 0.5*(asj+adj)
  }
}
cat("Expected (pedigree) relationship matrix A:\n"); print(round(A,3))
cat(sprintf("\n  full sibs   A,B : A = %.3f  (kinship %.3f)\n", A["A","B"], A["A","B"]/2))
cat(sprintf("  cousins     D,E : A = %.3f  (kinship %.3f)\n", A["D","E"], A["D","E"]/2))
cat(sprintf("  'unrelated' P1,C: A = %.3f   <-- pedigree says ZERO\n", A["P1","C"]))

## ---- (2) TRUE pedigree (with a HIDDEN shared ancestor X) -> gene-drop markers ----
# Truth the recorded pedigree didn't capture: P1 and C are actually half-sibs (share X).
L <- 1500                                  # SNP loci
draw <- function() rbind(rbinom(L,1,.5), rbinom(L,1,.5))   # 2 founder haplotypes
hap <- list(X=draw(), Y1=draw(), Y2=draw(), P2=draw(), Fdr=draw())
meiosis <- function(parent) parent[cbind(sample(1:2,L,TRUE), 1:L)]   # one gamete
cross <- function(s,d) rbind(meiosis(s), meiosis(d))
hap$P1 <- cross(hap$X,  hap$Y1)            # P1 = X x Y1
hap$C  <- cross(hap$X,  hap$Y2)            # C  = X x Y2  -> P1,C half-sibs (share X)
hap$A  <- cross(hap$P1, hap$P2)
hap$B  <- cross(hap$P1, hap$P2)
hap$D  <- cross(hap$A,  hap$C)
hap$E  <- cross(hap$B,  hap$Fdr)
hap_key <- c(P1="P1", P2="P2", C="C", F="Fdr", A="A", B="B", D="D", E="E")
M <- t(sapply(ped$id, function(k) colSums(hap[[ hap_key[[k]] ]])))   # 0/1/2 dosages
rownames(M) <- ped$id
# VanRaden G relative to the BASE population (founder allele freq p = 0.5), NOT the
# sample mean of 8 lines: Z = M - 2p; G = ZZ' / (2*sum p(1-p)).  This recovers the
# pedigree-scale relationships (full sib ~0.5, etc.).
p  <- 0.5
Z  <- M - 2*p
G  <- tcrossprod(Z) / (2 * ncol(M) * p * (1-p))          # realized relationship
cat("\nRealized (marker) relationship matrix G:\n"); print(round(G,3))
cat(sprintf("\n  full sibs   A,B : G = %.3f  (expected 0.50 -> Mendelian sampling)\n", G["A","B"]))
cat(sprintf("  cousins     D,E : G = %.3f  (expected 0.125)\n", G["D","E"]))
cat(sprintf("  'unrelated' P1,C: G = %.3f   <-- markers CATCH the hidden relatedness!\n", G["P1","C"]))

## ---- (3) figure: A heatmap | G heatmap | A-vs-G scatter ----
hm <- function(Mx,ttl) { d<-melt(Mx); d$Var1<-factor(d$Var1,levels=rev(ped$id))
  ggplot(d,aes(Var2,Var1,fill=value))+geom_tile(color="white")+
  geom_text(aes(label=sprintf("%.2f",value)),size=2.7)+
  scale_fill_gradient2(low="#2166AC",mid="white",high="#B2182B",midpoint=0,limits=c(-.3,1.3),name="")+
  labs(title=ttl,x=NULL,y=NULL)+theme_minimal(base_size=10)+
  theme(axis.text.x=element_text(angle=0)) }
pairs <- data.frame(
  A=A[upper.tri(A)], G=G[upper.tri(G)],
  lab=outer(ped$id,ped$id,paste,sep="-")[upper.tri(A)])
pairs$kind <- "other"
pairs$kind[pairs$lab%in%c("A-B")]   <- "full sibs"
pairs$kind[pairs$lab%in%c("D-E")]   <- "cousins"
pairs$kind[pairs$lab%in%c("P1-C")]  <- "'unrelated' (cryptic)"
sc <- ggplot(pairs,aes(A,G,color=kind))+
  geom_abline(slope=1,linetype="dashed",color="grey60")+
  geom_point(size=3)+ggrepel::geom_text_repel(aes(label=ifelse(kind=="other","",lab)),size=3.2,show.legend=FALSE)+
  scale_color_manual(values=c("full sibs"="#1b9e77","cousins"="#7570b3","'unrelated' (cryptic)"="#d73027","other"="grey70"),name=NULL)+
  labs(title="Pedigree expectation (A) vs marker reality (G)",
       subtitle="points off the dashed line = where DNA disagrees with the pedigree",
       x="A  (expected, from pedigree)",y="G  (realized, from markers)")+
  theme_bw(base_size=10)+theme(legend.position="top")
ggsave("../figures/22_toy_pedigree_kinship.png",
       (hm(A,"A — pedigree (expected)") | hm(G,"G — markers (realized)")) / sc + plot_layout(heights=c(1,1.1)),
       width=9, height=8.5, dpi=120)
cat("\nWrote figures/22_toy_pedigree_kinship.png\n")
