#!/usr/bin/env Rscript
# ============================================================================
# Module 1 reproduction: Explore the black bean dataset (GB_BLB.RData)
# Produces real numbers + figures used in the course lessons.
# ============================================================================
suppressMessages({library(ggplot2)})

repo <- "../repo"
load(file.path(repo, "GB_BLB.RData"))
dir.create("../figures", showWarnings = FALSE)

pheno <- GB_BLB$pheno
geno  <- as.matrix(GB_BLB$geno)
pos   <- GB_BLB$SNP_Position

cat("=== SHAPES ===\n")
cat("pheno :", nrow(pheno), "lines x", ncol(pheno), "cols\n")
cat("geno  :", nrow(geno), "lines x", ncol(geno), "SNPs\n")
cat("SNP_Position cols:", paste(colnames(pos), collapse=", "), "\n")
cat("pheno cols:", paste(colnames(pheno), collapse=", "), "\n\n")

# --- How many lines actually phenotyped per year? (breeding cycle structure) ---
cat("=== Non-missing phenotype counts (key traits) ===\n")
for (tr in c("yd18","yd19","app18","app19","text18","text19","col18","col19")) {
  cat(sprintf("  %-7s : %d lines\n", tr, sum(!is.na(pheno[[tr]]))))
}

# Breeding cycle 1 = first 272 rows (planted 2018 & 2019); cycle 2 = next 143 (2019 only)
cat("\nRows 1:272 with yd18 measured :", sum(!is.na(pheno$yd18[1:272])), "\n")
cat("Rows 273:415 with yd18 measured:", sum(!is.na(pheno$yd18[273:415])), "\n")
cat("Rows 273:415 with yd19 measured:", sum(!is.na(pheno$yd19[273:415])), "\n\n")

# ============================================================================
# FIGURE 1: Trait distributions across the 2 years (reproduces Supp Fig 1-2 idea)
# ============================================================================
traits_pairs <- list(
  c("yd18","yd19","Seed yield (kg/ha)"),
  c("sw18","sw19","100-seed weight (g)"),
  c("app18","app19","Canning appearance (1-5)"),
  c("text18","text19","Texture (peak force)")
)
long <- do.call(rbind, lapply(traits_pairs, function(p){
  rbind(
    data.frame(trait=p[3], year="2018", value=pheno[[p[1]]]),
    data.frame(trait=p[3], year="2019", value=pheno[[p[2]]]))
}))
long <- long[!is.na(long$value),]
ggplot(long, aes(value, fill=year)) +
  geom_density(alpha=.5) +
  facet_wrap(~trait, scales="free") +
  scale_fill_manual(values=c("2018"="#E69F00","2019"="#0072B2")) +
  labs(title="Trait distributions across 2 Michigan growing seasons",
       x="BLUP value", y="density") +
  theme_bw(base_size=12)
ggsave("../figures/01_trait_distributions.png", width=9, height=6, dpi=120)

# ============================================================================
# FIGURE 2: Trait correlation matrix (reproduces the negative yield-quality story)
# ============================================================================
trait_cols <- c("yd18","sw18","dm18","df18","app18","text18","col18","l18","a18","b18")
cm <- cor(pheno[,trait_cols], use="pairwise.complete.obs")
cat("=== Correlation: yield vs others (2018) ===\n")
print(round(cm["yd18",], 2))
cat("\nyield(yd18) vs seed weight(sw18):", round(cm["yd18","sw18"],2), "\n")
cat("seed weight(sw18) vs texture(text18):", round(cm["sw18","text18"],2), "\n\n")

cmdf <- as.data.frame(as.table(cm))
ggplot(cmdf, aes(Var1, Var2, fill=Freq)) +
  geom_tile(color="white") +
  geom_text(aes(label=sprintf("%.2f",Freq)), size=3) +
  scale_fill_gradient2(low="#2166AC", mid="white", high="#B2182B", midpoint=0, limits=c(-1,1)) +
  labs(title="Trait correlations (2018 BLUPs)", x="", y="", fill="r") +
  theme_minimal(base_size=11) + theme(axis.text.x=element_text(angle=45,hjust=1))
ggsave("../figures/02_trait_correlations.png", width=8, height=7, dpi=120)

# ============================================================================
# Genotype matrix sanity: coding, MAF
# ============================================================================
cat("=== Marker matrix ===\n")
cat("Unique genotype codes (sample):", paste(sort(unique(as.vector(geno[1:50,1:50]))), collapse=", "), "\n")
af <- colMeans(geno, na.rm=TRUE)/2   # markers coded 0/1/2 (allele dosage) -> allele freq = mean/2
maf <- pmin(af, 1-af)
cat("Allele coding range:", paste(range(geno, na.rm=TRUE), collapse=" to "), "\n")
cat("MAF summary:\n"); print(round(summary(maf),3))

png("../figures/03_maf_histogram.png", width=700, height=500, res=110)
hist(maf, breaks=40, col="#4DAF4A", border="white",
     main="Minor allele frequency (2,315 SNPs)", xlab="MAF")
abline(v=0.05, col="red", lwd=2, lty=2); dev.off()

# SNP count per chromosome
cat("\n=== SNPs per chromosome ===\n")
print(table(pos$Chromosome))

saveRDS(list(cm=cm, maf=maf), "../figures/explore_objs.rds")
cat("\nDONE. Figures written to ../figures/\n")
