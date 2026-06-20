# ============================================================================
# 11_gmatrix_construction_figure.R  —  Visual for Lesson 6.2/6.3 (building G)
# ----------------------------------------------------------------------------
# Connects the toy (4 lines x 3 SNPs, doses 0/1/2) to the two equations:
#   Step 1:  Z_ij = (M_ij - 2p_j) / s_j         (center & scale each SNP column)
#            ...with a full breakdown of how s_j is computed for EACH SNP.
#   Step 2:  G_ik = (1/p) * sum_j Z_ij * Z_kj    (average product of two Z-rows)
#
#   35 G_construction : M -> (per-SNP 2p_j & s_j) -> Z (one cell worked)
#                       -> dot-product of two rows -> G (incl. the colour heatmap).
#
# Run:  Rscript code/11_gmatrix_construction_figure.R
# ============================================================================

dir.create("figures", showWarnings = FALSE)

# palette (shared with 08/09/10 figure scripts)
col_dark  <- "#1b3a2b"; col_seed <- "#7d5a2b"; col_accent <- "#c0392b"
col_blue  <- "#2c6e9c"; col_grey <- "#888888"; col_gold <- "#d4a017"
col_line  <- "#eef2ee"

# ---- the toy: 4 lines x 3 SNPs, now WITH heterozygotes (dose 1) -------------
M <- matrix(c(2,2,2,
              2,1,2,
              1,1,0,
              1,0,0), nrow = 4, byrow = TRUE)
lines <- c("L1","L2","L3","L4"); snps <- c("SNP1","SNP2","SNP3")
n  <- nrow(M); p <- ncol(M)
twop <- colMeans(M)                         # 2p_j  = column mean        (1.5, 1.0, 1.0)
SS   <- apply(M, 2, function(x) sum((x - mean(x))^2))   # sum of squared devs (1, 2, 4)
varj <- SS / (n - 1)                         # /(n-1)                     (0.33,0.67,1.33)
s    <- sqrt(varj)                           # s_j  = column SD           (0.58,0.82,1.15)
Z    <- scale(M)                             # = (M - 2p_j)/s_j
G    <- tcrossprod(Z) / p                    # Z Z' / p

# ---- helpers ---------------------------------------------------------------
divcol <- function(v, maxabs = 1.0) {        # red +, blue -, white ~0
  t <- max(min(v / maxabs, 1), -1)
  if (t >= 0) colorRampPalette(c("white", col_accent))(101)[round(t * 100) + 1]
  else        colorRampPalette(c("white", col_blue ))(101)[round(-t * 100) + 1]
}
divtxt <- function(v) if (abs(v) > 0.6) "white" else col_dark

grid_mat <- function(x0, ytop, V, cw, ch, fmt = "%.2f", fillfun = NULL,
                     fixedfill = col_line, txtcex = 0.9, txtfun = NULL,
                     rowlab = NULL, collab = NULL, title = NULL) {
  nr <- nrow(V); nc <- ncol(V)
  if (!is.null(title))
    text(x0 + nc * cw / 2, ytop + 0.5, title, font = 2, cex = 0.92, col = col_dark)
  if (!is.null(collab)) for (c in 1:nc)
    text(x0 + (c - 0.5) * cw, ytop + 0.16, collab[c], cex = 0.62, col = col_grey)
  if (!is.null(rowlab)) for (r in 1:nr)
    text(x0 - 0.16, ytop - (r - 0.5) * ch, rowlab[r], cex = 0.74, font = 2,
         col = col_dark, adj = 1)
  for (r in 1:nr) for (c in 1:nc) {
    xl <- x0 + (c - 1) * cw; yt <- ytop - (r - 1) * ch
    f  <- if (!is.null(fillfun)) fillfun(V[r, c]) else fixedfill
    rect(xl, yt - ch, xl + cw, yt, col = f, border = "white", lwd = 1.6)
    tc <- if (!is.null(txtfun)) txtfun(V[r, c]) else col_dark
    text(xl + cw / 2, yt - ch / 2, sprintf(fmt, V[r, c]), col = tc, cex = txtcex)
  }
}
vcell <- function(xc, yc, v, cw, ch, fill, tcol = col_dark, fmt = "%+.2f", cex = 0.9) {
  rect(xc - cw/2, yc - ch/2, xc + cw/2, yc + ch/2, col = fill, border = "white", lwd = 1.6)
  text(xc, yc, sprintf(fmt, v), col = tcol, cex = cex)
}

# ============================================================================
# FIG 35 — the construction, step by step
# ============================================================================
png("figures/35_G_construction.png", width = 2150, height = 1780, res = 150)
par(mar = c(0.3, 0.6, 0.3, 0.6))
plot(NA, xlim = c(0, 21.5), ylim = c(0, 18.2), axes = FALSE, xlab = "", ylab = "")

text(10.7, 17.8, "Building G by hand: from genotypes M to relatedness G   (toy: 4 lines, 3 SNPs, doses 0/1/2)",
     font = 2, cex = 1.4, col = col_dark)

# ---- STEP 1 banner ---------------------------------------------------------
text(10.7, 16.95,
     expression(bold("Step 1") ~ "— center & scale every SNP column:" ~~
                  Z[ij] == (M[ij] - 2 * p[j]) / s[j] ~~ "  (" * i * "=line, " * j * "=SNP)"),
     cex = 1.0, col = col_dark)

cw <- 0.92; ch <- 0.78
Mtop <- 15.9; Mx <- 1.9
grid_mat(Mx, Mtop, M, cw, ch, fmt = "%.0f", fixedfill = "#e7ece7",
         rowlab = lines, collab = snps, title = "M  (genotypes 0/1/2)")
# mark the heterozygotes (dose 1) in M
for (r in 1:n) for (c in 1:p) if (M[r, c] == 1)
  rect(Mx + (c-1)*cw + 0.05, Mtop - r*ch + 0.05, Mx + c*cw - 0.05, Mtop - (r-1)*ch - 0.05,
       border = col_gold, lwd = 2.2)
text(Mx + 1.5*cw, Mtop - 4*ch - 0.18, "gold = heterozygote (dose 1)", cex = 0.6, col = col_gold)

# arrow M -> Z
arrows(Mx + 3*cw + 0.3, Mtop - 2*ch, Mx + 3*cw + 1.55, Mtop - 2*ch,
       lwd = 2.6, length = 0.12, col = col_grey)
text(Mx + 3*cw + 0.92, Mtop - 2*ch + 0.45, "use 2p_j, s_j\n(see table)", cex = 0.58, col = col_grey)

# Z matrix
Zx <- Mx + 3*cw + 1.9; Ztop <- 15.9
grid_mat(Zx, Ztop, round(Z, 2), cw, ch, fmt = "%+.2f", fillfun = divcol, txtfun = divtxt,
         rowlab = lines, collab = snps, title = "Z  (centered + scaled)")

# worked-cell callout (top-left cell) + the heterozygote nuance
calx <- Zx + 3*cw + 0.5
rect(Zx + 0.04, Ztop - ch + 0.04, Zx + cw - 0.04, Ztop - 0.04, border = col_dark, lwd = 2.6, xpd = NA)
rect(calx, Ztop - 2.55, calx + 9.6, Ztop + 0.12, col = adjustcolor(col_dark, 0.05),
     border = col_grey, lwd = 1.2)
text(calx + 0.25, Ztop - 0.35, "worked example - one cell:", adj = 0, cex = 0.74, col = col_dark, font = 2)
text(calx + 0.25, Ztop - 0.92,
     bquote(Z[11] == (M[11] - 2*p[1]) / s[1] ~ "= (2 - 1.5) / 0.58 = +0.87"),
     adj = 0, cex = 0.9, col = col_dark)
text(calx + 0.25, Ztop - 1.55, "what about a heterozygote (dose 1)?", adj = 0, cex = 0.72, col = col_gold, font = 2)
text(calx + 0.25, Ztop - 2.05,
     bquote("SNP2 (2p=1.0): dose 1 -> (1 - 1.0)/0.82 = " * bold("0") * "  (exactly average)"),
     adj = 0, cex = 0.8, col = col_dark)
text(calx + 0.25, Ztop - 2.42,
     bquote("SNP1 (2p=1.5): dose 1 -> (1 - 1.5)/0.58 = " * bold("-0.87") * "  (below average)"),
     adj = 0, cex = 0.8, col = col_dark)

# ---- STEP 1b: how s_j is computed, per SNP --------------------------------
text(2.0, 12.25, "How s_j (the SD of each SNP column) is built - shown for all 3 SNPs:",
     adj = 0, cex = 0.92, font = 2, col = col_dark)

# table geometry
labx <- 5.2                                  # right edge of the step-label column
colx <- c(7.0, 11.1, 15.2)                   # centres of SNP1/2/3 columns
ytab <- 11.35; rh <- 0.42
# background
rect(1.7, ytab - 6*rh - 0.05, 18.6, ytab + 0.55, col = adjustcolor(col_line, 0.6), border = col_grey, lwd = 1)
# header row (SNP names)
for (j in 1:3) text(colx[j], ytab + 0.28, snps[j], font = 2, cex = 0.82, col = col_blue)
text(labx, ytab + 0.28, "step", adj = 1, cex = 0.74, font = 3, col = col_grey)
# helper to draw one table row
devstr <- function(j) paste(sprintf("%+.1f", M[, j] - twop[j]), collapse = ", ")
dosstr <- function(j) paste(M[, j], collapse = ", ")
trow <- function(k, label, vals, col = col_dark, font = 1) {
  y <- ytab - (k - 0.5) * rh
  text(labx, y, label, adj = 1, cex = 0.72, col = col_grey)
  for (j in 1:3) text(colx[j], y, vals[j], cex = 0.72, col = col, font = font)
}
trow(1, "doses  M_ij",            sapply(1:3, dosstr))
trow(2, "mean = 2p_j",            sprintf("%.2f", twop), col = col_blue, font = 2)
trow(3, "deviations  M_ij - 2p_j", sapply(1:3, devstr))
trow(4, "sum of squares  S(.)^2",  sprintf("%.2f", SS))
trow(5, "/ (n-1) = / 3",           sprintf("%.3f", varj))
# highlight the final s_j row
yb <- ytab - 6*rh
rect(1.7, yb, 18.6, yb + rh, col = adjustcolor(col_gold, 0.16), border = NA)
trow(6, "s_j = sqrt( . )",         sprintf("%.2f", s), col = col_accent, font = 2)
# tie the s_j row back to Z
text(18.7, yb + rh/2, "->  these s_j\nfeed Z above", adj = 0, cex = 0.6, col = col_grey, xpd = NA)

# divider
segments(0.6, 8.55, 21.0, 8.55, col = adjustcolor(col_grey, 0.5), lwd = 1.2, lty = 3)

# ---- STEP 2 banner ---------------------------------------------------------
text(10.7, 8.15,
     expression(bold("Step 2") ~ "— one entry of G = average product of two lines' Z-rows:" ~~
                  G[ik] == (1/p) * sum(Z[ij] * Z[kj], j == 1, p) ~~ "  (" * p == 3 * " SNPs)"),
     cex = 1.0, col = col_dark)

rx0 <- 3.4; rcw <- 1.3; rch <- 0.78
for (c in 1:3) text(rx0 + (c-1)*rcw, 7.5,
                    bquote(italic(j) == .(c) ~ "(" * .(snps[c]) * ")"), cex = 0.64, col = col_grey)
# row i = L1
text(rx0 - 1.05, 6.85, expression(italic(i) == 1 ~ "(L1)"), cex = 0.78, font = 2, col = col_dark, adj = 1)
for (c in 1:3) vcell(rx0 + (c-1)*rcw, 6.85, Z[1, c], rcw - 0.18, rch, divcol(Z[1, c]), divtxt(Z[1, c]))
for (c in 1:3) text(rx0 + (c-1)*rcw, 6.30, "x", cex = 0.95, col = col_grey)
# row k = L2
text(rx0 - 1.05, 5.75, expression(italic(k) == 2 ~ "(L2)"), cex = 0.78, font = 2, col = col_dark, adj = 1)
for (c in 1:3) vcell(rx0 + (c-1)*rcw, 5.75, Z[2, c], rcw - 0.18, rch, divcol(Z[2, c]), divtxt(Z[2, c]))
for (c in 1:3) arrows(rx0 + (c-1)*rcw, 5.32, rx0 + (c-1)*rcw, 4.92, length = 0.06, lwd = 1.4, col = col_grey)
# products
prod12 <- Z[1, ] * Z[2, ]
for (c in 1:3) vcell(rx0 + (c-1)*rcw, 4.5, round(prod12[c], 2), rcw - 0.18, rch,
                     divcol(prod12[c]), divtxt(prod12[c]))
text(rx0 - 1.05, 4.5, expression(Z[1*j] %.% Z[2*j]), cex = 0.76, col = col_dark, adj = 1)
# call out the heterozygote (SNP2: L2 is dose 1 -> Z=0 -> product 0)
arrows(rx0 + 1*rcw, 4.5 - rch/2 - 0.05, rx0 + 1*rcw + 0.5, 3.75, length = 0.05, lwd = 1.2, col = col_gold)
text(rx0 + 1*rcw + 0.6, 3.65, "L2 is a het here (dose 1) -> Z=0 -> adds nothing",
     adj = 0, cex = 0.62, col = col_gold)

# sum + divide -> G12
text(rx0 + 1.0*rcw, 3.0,
     bquote(sum(Z[1*j]*Z[2*j], j, "") ~ "= +0.75 + 0.00 + 0.75 = +1.50"),
     cex = 0.88, col = col_dark)
text(rx0 + 1.0*rcw, 2.45,
     bquote("divide by p = 3   ->   " * G[12] == 0.50 ~ "(related)"),
     cex = 0.96, col = col_accent, font = 2)
# opposite case G14
text(rx0 - 1.05, 1.55, "the opposite extreme:", cex = 0.72, font = 2, col = col_dark, adj = 0)
text(rx0 - 1.05, 1.05,
     bquote(G[14] ~ "= (1/3)[ (-0.75) + (-1.50) + (-0.75) ] = -3.00/3 = -1.00  (opposite at every SNP)"),
     cex = 0.82, col = col_blue, adj = 0)

# ---- resulting G matrix (right) -------------------------------------------
Gx <- 15.2; Gtop <- 6.7; gcw <- 0.95; gch <- 0.9
grid_mat(Gx, Gtop, round(G, 2), gcw, gch, fmt = "%+.2f", fillfun = divcol, txtfun = divtxt,
         rowlab = lines, collab = lines, title = "G = Z Z' / p   (relatedness)")
rect(Gx + 1*gcw + 0.04, Gtop - gch + 0.04, Gx + 2*gcw - 0.04, Gtop - 0.04, border = col_accent, lwd = 3, xpd = NA)
rect(Gx + 3*gcw + 0.04, Gtop - gch + 0.04, Gx + 4*gcw - 0.04, Gtop - 0.04, border = col_blue,   lwd = 3, xpd = NA)
text(Gx + 2*gcw, Gtop - 4*gch - 0.4, "diagonal G_ii = how far a line sits from average",
     cex = 0.6, col = col_grey)
arrows(rx0 + 2.0*rcw, 2.45, Gx + 1.5*gcw, Gtop - 0.5*gch - 0.05,
       lwd = 1.6, length = 0.09, col = adjustcolor(col_accent, 0.6))

# ---- caption ---------------------------------------------------------------
text(10.7, 0.35,
     "The real G runs these same steps over 415 lines x 2,315 SNPs - only the matrices grow. Every off-diagonal G_ik is one dot product of two standardized rows, averaged over markers.",
     cex = 0.78, col = col_blue)

dev.off()
cat("wrote 35_G_construction.png\n")
