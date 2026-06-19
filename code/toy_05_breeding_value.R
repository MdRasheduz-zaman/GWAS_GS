#!/usr/bin/env Rscript
# Toy worked example for Lesson 5: 5 lines x 10 SNPs -> breeding value g = M %*% alpha
# Prints the explicit arithmetic and draws a labelled visual.
suppressMessages(library(ggplot2)); suppressMessages(library(reshape2))
set.seed(3)

lines <- paste0("L", 1:5)
snps  <- paste0("SNP", 1:10)
M <- rbind(
  L1 = c(2,0,1,2,0,2,0,0,1,2),
  L2 = c(0,0,2,2,0,0,1,2,1,0),
  L3 = c(2,1,0,0,2,2,0,0,0,2),
  L4 = c(0,2,2,1,0,0,2,2,0,0),
  L5 = c(1,0,0,2,1,2,0,1,2,2))
colnames(M) <- snps
# additive effects (the "price per copy" of the counted allele); chosen small for mental math
alpha <- c(SNP1=2, SNP2=0, SNP3=-1, SNP4=1, SNP5=0, SNP6=3, SNP7=0, SNP8=-2, SNP9=1, SNP10=0)

g <- as.vector(M %*% alpha)          # breeding value of each line
names(g) <- lines

cat("Genotype matrix M (rows = lines i, cols = SNPs j), entries x_ij in {0,1,2}:\n")
print(M)
cat("\nAdditive effects alpha_j (effect per copy):\n"); print(alpha)
cat("\nBreeding values g_i = sum_j alpha_j * x_ij :\n"); print(g)

cat("\n--- Worked arithmetic for line L1 ---\n")
terms <- M["L1",] * alpha
nz <- which(terms != 0)
cat("g_L1 =", paste(sprintf("%g*%g", alpha[nz], M["L1",nz]), collapse=" + "), "\n")
cat("     =", paste(sprintf("%g", terms[nz]), collapse=" + "), "=", sum(terms), "\n")

# ---- visualization: M heatmap + alpha bar + g bar ----
mm <- melt(M); colnames(mm) <- c("line","snp","x")
mm$line <- factor(mm$line, levels=rev(lines))
p1 <- ggplot(mm, aes(snp, line, fill=factor(x))) +
  geom_tile(color="white") + geom_text(aes(label=x), size=4) +
  scale_fill_manual(values=c("0"="#deebf7","1"="#9ecae1","2"="#3182bd"), name="x (dosage)") +
  labs(title="Step 1: genotypes  x[i,j]", x=NULL, y="line i") +
  theme_minimal(base_size=12) + theme(axis.text.x=element_text(angle=45,hjust=1))
ad <- data.frame(snp=factor(snps,levels=snps), alpha=alpha)
p2 <- ggplot(ad, aes(snp, alpha, fill=alpha>0)) + geom_col() +
  geom_text(aes(label=alpha), vjust=ifelse(alpha>=0,-0.3,1.2), size=4) +
  scale_fill_manual(values=c("TRUE"="#1b9e77","FALSE"="#d95f02"), guide="none") +
  scale_y_continuous(expand=expansion(mult=c(0.18,0.18))) +
  labs(title=expression("Step 2: effect per copy  "*alpha[j]), x=NULL, y=expression(alpha)) +
  theme_minimal(base_size=12) + theme(axis.text.x=element_text(angle=45,hjust=1))
gd <- data.frame(line=factor(lines,levels=lines), g=g)
p3 <- ggplot(gd, aes(line, g, fill=g)) + geom_col() + geom_text(aes(label=g), vjust=-0.3, size=4) +
  scale_fill_gradient2(low="#d95f02", mid="grey90", high="#1b9e77", guide="none") +
  scale_y_continuous(expand=expansion(mult=c(0.12,0.18))) +
  labs(title=expression("Step 3: breeding value  "*g[i]==sum(alpha[j]*x[i*j],j)), x="line i", y="g") +
  theme_minimal(base_size=12)

suppressMessages(library(patchwork))
ggsave("../figures/11_toy_breeding_value.png", (p1/p2/p3)+plot_layout(heights=c(2,1,1)),
       width=8, height=9, dpi=120)
cat("\nWrote figures/11_toy_breeding_value.png\n")
