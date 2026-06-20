# ============================================================================
# 10_breeding_value_dna_figure.R  —  Visual for Lesson 5.2 (breeding value)
# ----------------------------------------------------------------------------
# Makes the equation  g_i = sum_j alpha_j * x_ij  click by showing WHERE x_ij
# comes from: the actual DNA.
#
#   33 snp_to_breeding_value : 2 LINES (i), each diploid = 2 chromosome copies
#       drawn as DNA double-strands. 4 SNPs (j) sit at positions along the
#       chromosome. At each SNP we COUNT copies of the effect allele (red) over
#       the two copies -> x_ij in {0,1,2}. Multiply by the per-copy effect
#       alpha_j and sum across SNPs -> the breeding value g_i.
#
# Conceptual diagram (no data needed).  Run:
#   Rscript code/10_breeding_value_dna_figure.R
# ============================================================================

dir.create("figures", showWarnings = FALSE)

# palette (shared with 08/09 figure scripts)
col_dark  <- "#1b3a2b"; col_seed <- "#7d5a2b"; col_accent <- "#c0392b"
col_blue  <- "#2c6e9c"; col_grey <- "#888888"; col_gold <- "#d4a017"
col_line  <- "#eef2ee"

# ---- the worked example (2 lines x 4 SNPs) ---------------------------------
# effect allele counted at each SNP, and per-copy effect alpha_j
eff_allele <- c("A", "G", "C", "T")     # the allele whose copies we COUNT (red)
alpha      <- c( 3,  -2,   1,   2)      # additive effect per copy

# the two chromosome copies (homologs) of each line, as bases at the 4 SNPs
L1_copyA <- c("A", "C", "G", "T")
L1_copyB <- c("A", "C", "C", "T")
L2_copyA <- c("T", "G", "G", "A")
L2_copyB <- c("T", "G", "C", "A")

dosage <- function(a, b) as.integer(a == eff_allele) + as.integer(b == eff_allele)
xL1 <- dosage(L1_copyA, L1_copyB)        # -> 2 0 1 2
xL2 <- dosage(L2_copyA, L2_copyB)        # -> 0 2 1 0
gL1 <- sum(alpha * xL1)                   # -> 11
gL2 <- sum(alpha * xL2)                   # -> -3

# ---- geometry --------------------------------------------------------------
xlo <- 2.4; xhi <- 13.6                   # chromosome span
px  <- c(3.9, 6.5, 9.1, 11.7)             # the 4 SNP positions along it
xdiv <- 14.4                              # divider: DNA (left) | arithmetic (right)

# ---- drawing helpers -------------------------------------------------------
# one chromosome copy drawn as a DNA double-strand (two edges + ladder rungs)
draw_homolog <- function(yc, h = 0.24) {
  yt <- yc + h; yb <- yc - h
  # rounded backbone
  rect(xlo, yb, xhi, yt, col = col_line, border = NA)
  segments(xlo, yt, xhi, yt, col = col_dark, lwd = 2)
  segments(xlo, yb, xhi, yb, col = col_dark, lwd = 2)
  segments(xlo, yb, xlo, yt, col = col_dark, lwd = 2)
  segments(xhi, yb, xhi, yt, col = col_dark, lwd = 2)
  # ladder rungs (the "two strands" texture)
  for (rx in seq(xlo + 0.18, xhi - 0.18, by = 0.34))
    segments(rx, yb + 0.03, rx, yt - 0.03, col = adjustcolor(col_grey, 0.45), lwd = 1)
}

# the allele base sitting at a SNP: red disc if it is the counted effect allele
draw_base <- function(x, y, base, is_eff) {
  fill <- if (is_eff) col_accent else "#c9d3cb"
  tcol <- if (is_eff) "white"     else col_dark
  points(x, y, pch = 21, bg = fill, col = "white", cex = 3.0, lwd = 2)
  text(x, y, base, col = tcol, cex = 1.0, font = 2)
}

# draw one line's diploid pair + its dosage row
draw_line <- function(label, sub, yA, yB, copyA, copyB, xvec, ydose, gval, verdict, vcol) {
  draw_homolog(yA); draw_homolog(yB)
  for (k in 1:4) {
    draw_base(px[k], yA, copyA[k], copyA[k] == eff_allele[k])
    draw_base(px[k], yB, copyB[k], copyB[k] == eff_allele[k])
  }
  # left labels: which LINE (i) + "2 chromosome copies"
  text(xlo - 0.35, (yA + yB) / 2 + 0.18, label, adj = 1, font = 2,
       cex = 1.15, col = col_dark, xpd = NA)
  text(xlo - 0.35, (yA + yB) / 2 - 0.30, sub, adj = 1, cex = 0.66,
       col = col_grey, xpd = NA)
  text(xhi + 0.18, yA, "copy 1", adj = 0, cex = 0.62, col = col_grey, xpd = NA)
  text(xhi + 0.18, yB, "copy 2", adj = 0, cex = 0.62, col = col_grey, xpd = NA)
  # brace tying the two copies together
  segments(xlo - 0.30, yB - 0.05, xlo - 0.30, yA + 0.05, col = col_grey, lwd = 1.2, xpd = NA)
  # dosage boxes: x_ij = count of red copies
  for (k in 1:4) {
    cv <- xvec[k]
    bf <- if (cv == 0) "#eef2ee" else adjustcolor(col_accent, 0.10 + 0.10 * cv)
    rect(px[k] - 0.42, ydose - 0.34, px[k] + 0.42, ydose + 0.34,
         col = bf, border = col_accent, lwd = 1.6)
    text(px[k], ydose, cv, col = col_accent, font = 2, cex = 1.25)
  }
  text(xlo - 0.35, ydose, expression(x[ij]), adj = 1, cex = 0.95, col = col_accent, xpd = NA)
  text(xlo - 0.35, ydose - 0.45, "(# red)", adj = 1, cex = 0.6, col = col_grey, xpd = NA)
  # arrows from each copy-pair down into its dosage box
  for (k in 1:4)
    arrows(px[k], yB - 0.30, px[k], ydose + 0.40, length = 0.06,
           lwd = 1.4, col = adjustcolor(col_accent, 0.55))
  # RIGHT side: the dot product -> g_i
  eq <- paste(sprintf(ifelse(alpha < 0, "(%g)·%d", "%g·%d"),
                      alpha, xvec), collapse = " + ")
  yg <- (yA + yB) / 2
  text(xdiv + 0.2, yg + 0.55, bquote(g[.(sub2(label))] ~ "=" ~ sum(alpha[j] * x[ij], j, "")),
       adj = 0, cex = 0.95, col = col_dark, xpd = NA)
  text(xdiv + 0.2, yg + 0.02, eq, adj = 0, cex = 0.92, col = col_dark, xpd = NA)
  text(xdiv + 0.2, yg - 0.5, sprintf("= %d", gval), adj = 0, cex = 1.25,
       col = vcol, font = 2, xpd = NA)
  text(xdiv + 2.4, yg - 0.5, verdict, adj = 0, cex = 0.95, col = vcol, font = 2, xpd = NA)
}
sub2 <- function(label) gsub("Line ", "", label)   # "Line L1" -> "L1"

# ============================================================================
# FIG 33 — from SNPs on the DNA to the breeding value g_i
# ============================================================================
png("figures/33_snp_to_breeding_value.png", width = 2000, height = 1300, res = 150)
par(mar = c(0.5, 1, 0.5, 1))
plot(NA, xlim = c(0, 20), ylim = c(0, 15.6), axes = FALSE, xlab = "", ylab = "")

# ---- title + legend --------------------------------------------------------
text(10, 15.1, "Where does  x  come from?  Reading the breeding value off the DNA",
     font = 2, cex = 1.5, col = col_dark)
mtext_y <- 14.5
points(1.2, mtext_y, pch = 21, bg = col_accent, col = "white", cex = 1.8)
text(1.55, mtext_y, "= the EFFECT allele we count (red)", adj = 0, cex = 0.8, col = col_dark)
points(7.4, mtext_y, pch = 21, bg = "#c9d3cb", col = "white", cex = 1.8)
text(7.7, mtext_y, "= the other allele (not counted)", adj = 0, cex = 0.8, col = col_dark)
text(13.2, mtext_y, "x = how many red copies on the 2 chromosome copies  -> 0, 1, or 2",
     adj = 0, cex = 0.8, col = col_accent)

# divider: DNA / dosage (left)  vs  arithmetic (right)
segments(xdiv - 0.3, 1.2, xdiv - 0.3, 13.6, col = adjustcolor(col_grey, 0.5),
         lwd = 1.2, lty = 3)
text(xdiv - 0.3, 13.75, "read x off the strands", cex = 0.62, col = col_grey, adj = 1.02)
text(xdiv - 0.1, 13.75, "then do the arithmetic", cex = 0.62, col = col_grey, adj = -0.02)

# ---- Line 1 (top) ----------------------------------------------------------
draw_line("Line L1", "L1", yA = 13.0, yB = 12.2, L1_copyA, L1_copyB,
          xL1, ydose = 11.1, gL1, "-> keep (best)", col_dark)

# ---- shared chromosome axis (j) + effect row (alpha_j) ---------------------
yax <- 9.85
rect(xlo, yax - 0.12, xhi, yax + 0.12, col = "#dfe7df", border = col_dark, lwd = 1.2)
for (k in 1:4) {
  segments(px[k], yax - 0.20, px[k], yax + 0.20, col = col_dark, lwd = 2)
  text(px[k], yax + 0.55, bquote("SNP " ~ italic(j) == .(k)), cex = 0.85, font = 2, col = col_dark)
  # alpha box under the axis
  ay <- 9.0
  text(px[k], ay + 0.30, sprintf("effect allele = %s", eff_allele[k]),
       cex = 0.62, col = col_grey)
  acol <- if (alpha[k] > 0) col_blue else if (alpha[k] < 0) col_accent else col_grey
  rect(px[k] - 0.55, ay - 0.40, px[k] + 0.55, ay + 0.05, col = adjustcolor(acol, 0.12),
       border = acol, lwd = 1.6)
  text(px[k], ay - 0.18, bquote(alpha[.(k)] == .(sprintf("%+d", alpha[k]))),
       cex = 0.95, font = 2, col = acol)
}
text(xlo - 0.35, yax, "chromosome", adj = 1, cex = 0.72, col = col_dark, font = 3, xpd = NA)
text(xlo - 0.35, 9.0 - 0.18, expression(alpha[j]), adj = 1, cex = 0.95, col = col_dark, xpd = NA)
text(xlo - 0.35, 9.0 - 0.6, "(per copy)", adj = 1, cex = 0.6, col = col_grey, xpd = NA)
text(xhi + 0.18, yax, "(11 of these\nin the real genome)", adj = 0, cex = 0.6,
     col = col_grey, xpd = NA)

# ---- Line 2 (bottom) -------------------------------------------------------
draw_line("Line L2", "L2", yA = 7.7, yB = 6.9, L2_copyA, L2_copyB,
          xL2, ydose = 5.8, gL2, "-> cull (worst)", col_accent)

# ---- bottom strip: the equation + what i, j, x mean ------------------------
yb0 <- 4.0
rect(0.4, 0.4, 19.6, yb0, col = adjustcolor(col_dark, 0.05), border = col_grey, lwd = 1)
text(10, yb0 - 0.45,
     expression(bold("The breeding (genetic) value of line ") * italic(i) * bold(":") ~~
                  g[i] == sum(alpha[j] * x[ij], j == 1, p)),
     cex = 1.25, col = col_dark)
text(10, yb0 - 1.15,
     expression(italic(i) ~ "= which LINE (plant)" ~~ "   " ~~
                  italic(j) ~ "= which SNP (position on the chromosome)"),
     cex = 0.86, col = col_grey)
text(10, yb0 - 1.62,
     expression(x[ij] ~ "= effect-allele copies of line " * italic(i) * " at SNP " * italic(j) ~
                  "(0/1/2)" ~~ "   " ~~ alpha[j] ~ "= effect per copy"),
     cex = 0.86, col = col_grey)
text(10, yb0 - 2.35,
     sprintf("g(L1) = 3·2 + (-2)·0 + 1·1 + 2·2 = %d        g(L2) = 3·0 + (-2)·2 + 1·1 + 2·0 = %d",
             gL1, gL2),
     cex = 0.92, col = col_dark, font = 2)
text(10, yb0 - 3.0,
     "Same machinery as the 5x10 toy table in §5.2 - here every x is read straight off the chromosome copies. Stack all p SNPs -> g = Mα.",
     cex = 0.78, col = col_blue)

dev.off()
cat("wrote 33_snp_to_breeding_value.png\n")
