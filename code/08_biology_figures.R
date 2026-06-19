# ============================================================================
# 08_biology_figures.R  —  Visual biology/breeding figures for Lessons 1, 3-5
# ----------------------------------------------------------------------------
# Produces six figures that make the *biology* click:
#   25 chromosome ideogram        (REAL: GB_BLB$SNP_Position)
#   26 LD decay curve             (REAL: GB_BLB$geno + positions)
#   27 trait antagonism + select  (REAL: GB_BLB$pheno sw18 vs text18)
#   28 selfing -> homozygosity     (computed: the F1..F8 selfing series)
#   29 GBS in-silico digest        (simulated ApeKI cut + size-selection window)
#   30 bean flower: self vs cross  (schematic of why beans self-pollinate)
#
# Run:  Rscript code/08_biology_figures.R
# ============================================================================

load("repo/GB_BLB.RData")
SNPpos <- GB_BLB$SNP_Position
geno   <- as.matrix(GB_BLB$geno)
pheno  <- GB_BLB$pheno
dir.create("figures", showWarnings = FALSE)

# tidy chromosome labels: "01".."11" -> 1..11 ; ensure Position numeric
SNPpos$chr_num  <- as.integer(as.character(SNPpos$Chromosome))
SNPpos$Position <- as.numeric(as.character(SNPpos$Position))
SNPpos$pos_mb   <- SNPpos$Position / 1e6
chrs <- sort(unique(SNPpos$chr_num))

# palette
col_dark  <- "#1b3a2b"; col_seed <- "#7d5a2b"; col_accent <- "#c0392b"
col_blue  <- "#2c6e9c"; col_grey <- "#888888"; col_gold <- "#d4a017"

# ============================================================================
# FIG 25 — CHROMOSOME IDEOGRAM (real SNP positions on 11 chromosomes)
# ============================================================================
png("figures/25_chromosome_map.png", width = 1500, height = 1050, res = 150)
par(mar = c(5, 5, 4, 2))
chrlen <- tapply(SNPpos$pos_mb, SNPpos$chr_num, max)
plot(NA, xlim = c(0, max(chrlen) * 1.14), ylim = c(0.5, length(chrs) + 0.5),
     yaxt = "n", xlab = "Position along chromosome (millions of base pairs, Mb)",
     ylab = "", main = "The bean genome: 11 chromosomes, 2,315 SNP signposts",
     cex.lab = 1.15, cex.main = 1.35, bty = "n")
axis(2, at = chrs, labels = sprintf("Chr %02d", chrs), las = 1, cex.axis = 1.05)
for (k in seq_along(chrs)) {
  c <- chrs[k]; sub <- SNPpos[SNPpos$chr_num == c, ]
  y <- length(chrs) - k + 1                      # chr01 on top
  # chromosome body (rounded capsule)
  L <- chrlen[as.character(c)]
  rect(0, y - 0.28, L, y + 0.28, col = "#eef2ee", border = col_dark, lwd = 1.4)
  # each SNP a thin tick
  segments(sub$pos_mb, y - 0.28, sub$pos_mb, y + 0.28, col = col_dark, lwd = 0.35)
  # count label
  text(L + max(chrlen) * 0.005, y, sprintf("%d SNPs", nrow(sub)),
       adj = 0, cex = 0.8, col = col_grey)
}
mtext("Each vertical tick = one SNP. Dense regions = many signposts; gaps = few (e.g. centromeres, or low-cut-site stretches).",
      side = 1, line = 4.3, cex = 0.9, col = col_grey)
dev.off()
cat("wrote 25_chromosome_map.png\n")

# ============================================================================
# FIG 26 — LD DECAY (real): r^2 between SNP pairs vs distance apart
# ============================================================================
maxd_kb <- 1000                                  # look within 1 Mb
recs <- list()
for (c in chrs) {
  idx <- which(SNPpos$chr_num == c)
  if (length(idx) < 5) next
  G <- geno[, idx, drop = FALSE]
  pos <- SNPpos$Position[idx]
  R <- suppressWarnings(cor(G, use = "pairwise.complete.obs"))
  R2 <- R^2
  np <- length(idx)
  for (i in 1:(np - 1)) {
    j <- (i + 1):np
    d <- abs(pos[j] - pos[i]) / 1000             # kb
    keep <- d <= maxd_kb
    if (any(keep)) recs[[length(recs) + 1]] <-
      data.frame(d = d[keep], r2 = R2[i, j][keep])
  }
}
ld <- do.call(rbind, recs)
ld <- ld[is.finite(ld$r2), ]
# bin and average
brks <- seq(0, maxd_kb, by = 25)
ld$bin <- cut(ld$d, breaks = brks, labels = FALSE)
binmid <- (brks[-1] + brks[-length(brks)]) / 2
binmean <- tapply(ld$r2, ld$bin, mean)

png("figures/26_ld_decay.png", width = 1400, height = 1000, res = 150)
par(mar = c(5, 5, 4, 2))
set.seed(1)
samp <- ld[sample(nrow(ld), min(8000, nrow(ld))), ]
plot(samp$d, samp$r2, pch = 16, col = adjustcolor(col_grey, 0.18), cex = 0.5,
     xlab = "Distance between two SNPs on the same chromosome (kb)",
     ylab = expression("Linkage disequilibrium  r"^2),
     main = "Linkage disequilibrium decays with distance (real bean data)",
     cex.lab = 1.15, cex.main = 1.3, ylim = c(0, 1))
points(binmid[as.integer(names(binmean))], binmean, type = "b",
       pch = 19, col = col_blue, lwd = 2.5, cex = 1.1)
abline(h = 0.9, lty = 2, col = col_accent, lwd = 2)
text(maxd_kb * 0.62, 0.93, expression("LD-prune threshold  r"^2 * " > 0.9 (redundant)"),
     col = col_accent, cex = 0.95, adj = 0)
# half-decay annotation
h0 <- binmean[1]
legend("topright", bty = "n", cex = 0.95,
       legend = c("each dot = one SNP pair", "average r² per 25-kb bin"),
       pch = c(16, 19), col = c(col_grey, col_blue))
mtext("Close SNPs are inherited together (high r²) -> one 'tag' SNP speaks for its neighbours. Far-apart SNPs are nearly independent.",
      side = 1, line = 3.7, cex = 0.88, col = col_grey)
dev.off()
cat("wrote 26_ld_decay.png  (", nrow(ld), "pairs )\n")

# ============================================================================
# FIG 27 — TRAIT ANTAGONISM (real): 100-seed weight vs canning texture
# ============================================================================
d <- pheno[, c("sw18", "text18")]
d <- d[complete.cases(d), ]
r <- cor(d$sw18, d$text18)
png("figures/27_trait_antagonism.png", width = 1400, height = 1000, res = 150)
par(mar = c(5, 5, 4, 2))
plot(d$sw18, d$text18, pch = 21, bg = adjustcolor(col_seed, 0.5), col = "white",
     cex = 1.3, xlab = "100-seed weight (g)  — a YIELD component  → bigger seed",
     ylab = "Canning texture (peak force)  — higher = firmer = better",
     main = sprintf("Yield vs quality pull apart: seed weight ↔ texture  (r = %.2f)", r),
     cex.lab = 1.1, cex.main = 1.25)
ab <- lm(text18 ~ sw18, data = d)
abline(ab, col = col_accent, lwd = 3)
# show the selection tension: select heaviest-seed (top yield component) 20%
thr <- quantile(d$sw18, 0.8)
abline(v = thr, lty = 2, col = col_blue, lwd = 2)
sel <- d[d$sw18 >= thr, ]
points(sel$sw18, sel$text18, pch = 21, bg = adjustcolor(col_blue, 0.6), col = "white", cex = 1.4)
rect(thr, par("usr")[3], par("usr")[2], par("usr")[4],
     col = adjustcolor(col_blue, 0.06), border = NA)
text(thr, par("usr")[4] - 0.04 * diff(par("usr")[3:4]),
     "select heaviest seed\n(chase yield) →", col = col_blue, adj = 1.05, cex = 0.95)
arrows(thr + 0.7 * (par("usr")[2] - thr), mean(d$text18) + 0.5 * sd(d$text18),
       thr + 0.7 * (par("usr")[2] - thr), mean(sel$text18),
       col = col_accent, lwd = 2.5, length = 0.12)
text(thr + 0.7 * (par("usr")[2] - thr), mean(d$text18) + 0.7 * sd(d$text18),
     "...and texture\nfollows DOWN", col = col_accent, adj = 0.5, cex = 0.92)
mtext("The negative slope IS the antagonism: picking big seed for yield silently drags canning texture toward mush. This is why you need multi-trait methods (L12).",
      side = 1, line = 3.7, cex = 0.85, col = col_grey)
dev.off()
cat("wrote 27_trait_antagonism.png\n")

# ============================================================================
# FIG 28 — SELFING -> HOMOZYGOSITY (computed): how a LINE gets fixed
# ============================================================================
gen <- 1:8                                       # F1 (Aa) selfed onward
het <- (1/2)^(gen - 1)                            # heterozygosity after each self
het[1] <- 1                                       # F1 = 100% Aa at a segregating locus
homo <- 1 - het
png("figures/28_selfing_homozygosity.png", width = 1500, height = 750, res = 150)
par(mfrow = c(1, 2), mar = c(5, 5, 4, 1))
# panel A: the selfing series
bx <- barplot(rbind(homo, het) * 100, names.arg = paste0("F", gen),
        col = c(col_dark, col_gold), border = "white",
        ylab = "% of segregating loci", xlab = "Generation of self-pollination",
        main = "Each self-generation halves heterozygosity",
        cex.lab = 1.05, cex.main = 1.15, ylim = c(0, 122))
legend("top", horiz = TRUE, fill = c(col_dark, col_gold), bg = "white",
       box.col = NA, cex = 0.78, inset = 0.005, xpd = NA,
       legend = c("homozygous (AA/aa)", "heterozygous (Aa)"))
text(bx, homo * 100 + 4, sprintf("%.0f%%", het * 100), cex = 0.8, col = col_gold)
# mark the F6 generation used in the study (inside the plot)
arrows(bx[6], 84, bx[6], 97, col = col_accent, lwd = 2.5, length = 0.1)
text(bx[6], 79, "F6\n(study)", col = col_accent, cex = 0.9, font = 2)
# panel B: one-locus Punnett of Aa x Aa (self)
plot(NA, xlim = c(0, 3), ylim = c(0, 3), axes = FALSE, xlab = "", ylab = "",
     main = "Self Aa × Aa  →  1/4 AA : 1/2 Aa : 1/4 aa", cex.main = 1.15)
cells <- matrix(c("AA","Aa","Aa","aa"), 2, byrow = TRUE)
cols  <- matrix(c(col_dark, col_gold, col_gold, col_dark), 2, byrow = TRUE)
for (i in 1:2) for (j in 1:2) {
  rect(j, 3 - i, j + 1, 4 - i, col = adjustcolor(cols[i, j], 0.85), border = "white", lwd = 3)
  text(j + 0.5, 3.5 - i, cells[i, j], col = "white", cex = 1.8, font = 2)
}
text(0.5, 2.5, "A"); text(0.5, 1.5, "a"); text(1.5, 3.2, "A"); text(2.5, 3.2, "a")
text(1.5, 0.4, "pollen (self)", cex = 0.9, col = col_grey)
text(0.2, 2, "ovule", srt = 90, cex = 0.9, col = col_grey)
text(1.5, -0.05, "Half the offspring stay Aa; half become fixed.\nRepeat ≈ 6× → a near-pure, true-breeding LINE.",
     cex = 0.9, col = col_dark, xpd = TRUE)
dev.off()
cat("wrote 28_selfing_homozygosity.png\n")

# ============================================================================
# FIG 29 — GBS IN-SILICO DIGEST: restriction cut -> fragments -> size select
# ============================================================================
set.seed(7)
# ApeKI cuts at GCWGC ~ every few hundred bp; fragment lengths ~ exponential
nfrag <- 60000
frag <- round(rexp(nfrag, rate = 1/350))         # mean ~350 bp
frag <- frag[frag > 20 & frag < 2000]
lo <- 200; hi <- 350                              # size-selection window
png("figures/29_gbs_digest.png", width = 1650, height = 760, res = 150)
par(mfrow = c(1, 2), mar = c(5, 4.5, 5, 1.5), xpd = NA)
# panel A: schematic of cutting one DNA molecule
plot(NA, xlim = c(0, 10), ylim = c(0, 6), axes = FALSE, xlab = "", ylab = "",
     main = "1) Cut DNA at every ApeKI site (GCWGC)", cex.main = 1.1)
segments(0.3, 5, 9.7, 5, lwd = 6, col = col_blue)
# 6 cuts -> 7 fragments with chosen lengths; dark = inside 200-350 window
fx     <- c(0.3, 1.6, 2.5, 4.4, 5.4, 6.1, 8.0, 9.7)
fbp    <- c(280,  90,  320, 150,  70, 210, 130)        # bp of each fragment
cutx   <- fx[2:7]
for (x in cutx) { points(x, 5.55, pch = 25, bg = col_accent, col = col_accent, cex = 1.1)
                  segments(x, 4.7, x, 5.35, col = col_accent, lwd = 2) }
text(cutx[1], 5.95, "↓ cut sites", cex = 0.7, col = col_accent, adj = 0)
text(5, 6.4, "...GCWGC......GCWGC......GCWGC...", cex = 0.78, col = col_grey)
for (k in 1:(length(fx) - 1)) {
  kept <- fbp[k] >= lo && fbp[k] <= hi
  col  <- if (kept) col_dark else col_grey
  rect(fx[k] + 0.06, 3.3, fx[k+1] - 0.06, 3.9, col = adjustcolor(col, 0.8), border = "white")
  text((fx[k] + fx[k+1])/2, 3.05, paste0(fbp[k], "bp"), cex = 0.6, col = col)
}
text(5, 2.4, "2) Keep only fragments in a size window (200-350 bp)",
     cex = 0.92, col = col_dark, font = 2)
text(5, 1.85, "dark = sequenced (in window)   grey = discarded", cex = 0.74, col = col_grey)
text(5, 0.9, "Same enzyme + same window in EVERY plant\n= the SAME ~1% genome slice each time (reproducible).",
     cex = 0.82, col = col_blue)
# panel B: the size distribution + window
h <- hist(frag, breaks = 60, plot = FALSE)
plot(h, col = ifelse(h$mids >= lo & h$mids <= hi, col_dark, "#dddddd"),
     border = "white", xlim = c(0, 1200),
     xlab = "Restriction fragment length (bp)", ylab = "number of fragments",
     main = "3) Fragment sizes:\nwindow keeps ~1% to sequence", cex.lab = 1.0, cex.main = 1.05)
abline(v = c(lo, hi), col = col_accent, lwd = 2, lty = 2)
text((lo+hi)/2, max(h$counts)*0.92, "size-selection\nwindow", col = col_accent, cex = 0.85)
arrows(lo, max(h$counts)*0.72, hi, max(h$counts)*0.72, code = 3, col = col_accent, length = 0.08)
dev.off()
cat("wrote 29_gbs_digest.png\n")

# ============================================================================
# FIG 30 — BEAN FLOWER: why beans self-pollinate (schematic) + crossing
# ============================================================================
draw_flower <- function(cx, cy, s = 1, open = FALSE) {
  # banner (large upper petal)
  banner <- col_gold
  polygon(cx + s*c(-1.1, -0.7, 0, 0.7, 1.1, 0),
          cy + s*c(0.4, 1.5, 1.9, 1.5, 0.4, 0.3),
          col = adjustcolor(col_accent, 0.55), border = col_accent, lwd = 2)
  # wings
  polygon(cx + s*c(-1.1,-1.6,-1.3,-0.6), cy + s*c(0.3,-0.2,-0.7,-0.1),
          col = adjustcolor(col_gold,0.6), border = col_gold, lwd=1.5)
  polygon(cx + s*c(1.1,1.6,1.3,0.6), cy + s*c(0.3,-0.2,-0.7,-0.1),
          col = adjustcolor(col_gold,0.6), border = col_gold, lwd=1.5)
  # keel (boat-shaped, encloses sex organs)
  keelcol <- if (open) adjustcolor(col_seed,0.25) else adjustcolor(col_dark,0.55)
  polygon(cx + s*c(-0.55,0,0.55,0.35,-0.35), cy + s*c(-0.1,-1.0,-0.1,0.15,0.15),
          col = keelcol, border = col_dark, lwd = 2)
  if (open) {
    # stamens + stigma exposed
    for (a in seq(-0.3,0.3,length.out=7))
      segments(cx+s*a, cy-s*0.2, cx+s*a*2.0, cy-s*0.9, col=col_gold, lwd=1.5)
    points(cx, cy-s*0.95, pch=8, col=col_accent, cex=1.4) # stigma exposed
  }
}
png("figures/30_bean_flower_selfing.png", width = 1550, height = 760, res = 150)
par(mfrow = c(1, 2), mar = c(2, 2, 4, 2))
# panel A — selfing (closed keel)
plot(NA, xlim = c(-3, 3), ylim = c(-3, 3), axes = FALSE, xlab="", ylab="",
     main = "DEFAULT: self-pollination")
draw_flower(0, 0.3, 1.1, open = FALSE)
text(0, -1.9, "The KEEL petal encloses the\nanthers + stigma. Pollen is shed\nINSIDE the closed flower, often\nbefore it even opens.", cex = 0.92, col = col_dark)
text(0, 2.7, "→ each plant fertilises itself → offspring ≈ genetically identical to parent",
     cex = 0.82, col = col_blue, xpd = TRUE)
# panel B — controlled cross (breeder intervenes)
plot(NA, xlim = c(-3, 3), ylim = c(-3, 3), axes = FALSE, xlab="", ylab="",
     main = "BREEDER'S CROSS: force an outcross")
draw_flower(-1.2, 0.3, 0.95, open = TRUE)
text(-1.2, -1.7, "Parent 1 (mother):\nopen keel, snip out anthers\n(EMASCULATE) before pollen sheds", cex = 0.78, col = col_dark)
draw_flower(1.6, 1.1, 0.55, open = TRUE)
text(1.7, -0.2, "Parent 2\n(father):\npollen donor", cex = 0.78, col = col_dark)
arrows(1.2, 0.9, -0.4, 0.1, col = col_accent, lwd = 3, length = 0.15)
text(0.6, 0.9, "brush pollen", col = col_accent, cex = 0.8, srt = -20)
text(0, 2.7, "→ F1 hybrid combines TWO parents' alleles → then self for ~6 gens to fix new LINES",
     cex = 0.78, col = col_blue, xpd = TRUE)
dev.off()
cat("wrote 30_bean_flower_selfing.png\n")

cat("\nAll biology figures written to figures/25..30\n")
