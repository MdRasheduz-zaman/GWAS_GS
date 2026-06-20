# ============================================================================
# 09_experimental_design_figures.R  —  Visuals for Lesson 3.2 (field design)
# ----------------------------------------------------------------------------
# Makes the experimental design click:
#   31 block designs   : INCOMPLETE (2018) vs COMPLETE (2019) block layouts,
#                        with the 3 check cultivars (Eclipse/Zorro/Zenith)
#                        threaded through every block, and a replicated line
#                        traced across two blocks.
#   32 plot + pipeline : anatomy of one 4-row plot (harvest the centre 2) and
#                        how plot -> block -> checks -> model -> BLUP connect.
#
# Schematic (real 2018 trial had 200 lines + 72 once; far more plots). No data
# needed — these are conceptual diagrams.
#
# Run:  Rscript code/09_experimental_design_figures.R
# ============================================================================

dir.create("figures", showWarnings = FALSE)

# palette (shared with 08_biology_figures.R)
col_dark  <- "#1b3a2b"; col_seed <- "#7d5a2b"; col_accent <- "#c0392b"
col_blue  <- "#2c6e9c"; col_grey <- "#888888"; col_gold <- "#d4a017"
col_line  <- "#eef2ee"

# the 3 check cultivars -> one colour each (E=Eclipse, Z=Zorro, N=Zenith)
check_col <- c(E = col_accent, Z = col_blue, N = col_gold)
is_check  <- function(lab) lab %in% names(check_col)

# ----- field-drawing helpers ------------------------------------------------
draw_field <- function(M, x0, y0, cw = 1, ch = 1, cex = 0.82) {
  nr <- nrow(M); nc <- ncol(M)
  for (r in 1:nr) for (cc in 1:nc) {
    lab <- M[r, cc]
    xl  <- x0 + (cc - 1) * cw
    yb  <- y0 + (nr - r) * ch
    chk <- is_check(lab)
    fill <- if (chk) check_col[[lab]] else col_line
    rect(xl, yb, xl + cw, yb + ch, col = fill, border = "white", lwd = 1.6)
    text(xl + cw / 2, yb + ch / 2, lab,
         col = if (chk) "white" else col_dark,
         cex = cex, font = if (chk) 2 else 1)
  }
}
# outline a block spanning rows rtop..(rtop+rh-1), cols cleft..(cleft+cwn-1)
block_rect <- function(x0, y0, nr, rtop, cleft, rh, cwn, cw = 1, ch = 1,
                       border = col_dark, lwd = 3) {
  xl <- x0 + (cleft - 1) * cw
  xr <- x0 + (cleft - 1 + cwn) * cw
  yt <- y0 + (nr - (rtop - 1)) * ch
  yb <- y0 + (nr - (rtop - 1) - rh) * ch
  rect(xl, yb, xr, yt, border = border, lwd = lwd, col = NA, xpd = NA)
}
# highlight one cell (row r from top, col cc)
hl_cell <- function(x0, y0, nr, r, cc, cw = 1, ch = 1, col = col_accent) {
  xl <- x0 + (cc - 1) * cw; yb <- y0 + (nr - r) * ch
  rect(xl + 0.05, yb + 0.05, xl + cw - 0.05, yb + ch - 0.05,
       border = col, lwd = 3.2, xpd = NA)
}

# ============================================================================
# FIG 31 — INCOMPLETE (2018) vs COMPLETE (2019) block designs
# ============================================================================
# Same 9 lines (L1..L9) + 3 checks. INCOMPLETE: six small 2x2 blocks, each
# holds only a SUBSET of lines. COMPLETE: two big blocks (reps), each holds
# EVERY line once. Checks sit in every block in both.
M_inc <- rbind(
  c("L1","L2","L4","L5","L7","L8"),
  c("L3","E" ,"L6","Z" ,"L9","N" ),
  c("L1","L4","L2","L5","L3","L6"),
  c("L7","Z" ,"L8","N" ,"L9","E" ))
M_comp <- rbind(
  c("L3","E" ,"L7","L5","N" ,"L2"),
  c("L1","L5","Z" ,"Z" ,"L8","L4"),
  c("L8","L2","N" ,"L6","L1","E" ),
  c("L6","L9","L4","L9","L3","L7"))

png("figures/31_block_designs.png", width = 1850, height = 980, res = 150)
par(mar = c(3.5, 1, 5.5, 1))
plot(NA, xlim = c(0, 16.6), ylim = c(-1.4, 7.4), axes = FALSE, xlab = "", ylab = "")

mtext("The experimental design: how plots are grouped into blocks",
      side = 3, line = 3.3, cex = 1.45, font = 2)
mtext("3 check cultivars  —  Eclipse (E) · Zorro (Z) · Zenith (N)  —  are planted in EVERY block, so the model can read the field's gradient.",
      side = 3, line = 1.6, cex = 0.92, col = col_grey)

# ---- LEFT: incomplete block design (2018) ----
xL <- 0.6; yL <- 1.1; nr <- 4
draw_field(M_inc, xL, yL)
for (b in list(c(1,1),c(1,3),c(1,5),c(3,1),c(3,3),c(3,5)))
  block_rect(xL, yL, nr, rtop = b[1], cleft = b[2], rh = 2, cwn = 2,
             border = col_dark, lwd = 3)
text(xL + 3, yL + 5.5, "2018 — INCOMPLETE block design", font = 2, cex = 1.05, col = col_dark)
text(xL + 3, yL + 4.95, "six SMALL blocks — each holds only a SUBSET of the lines",
     cex = 0.82, col = col_grey)
# trace one replicated line across two different blocks
hl_cell(xL, yL, nr, 1, 1); hl_cell(xL, yL, nr, 3, 1)
arrows(xL + 0.5, yL + 3.5, xL + 0.5, yL + 1.5, col = col_accent, lwd = 2,
       length = 0.08, code = 3)
text(xL + 3, yL - 0.55, "→ line L1 is grown in two different blocks (replication)",
     col = col_accent, cex = 0.72)

# ---- RIGHT: complete block design (2019) ----
xR <- 9.4; yR <- 1.1
draw_field(M_comp, xR, yR)
block_rect(xR, yR, nr, rtop = 1, cleft = 1, rh = 4, cwn = 3, border = col_blue, lwd = 3.5)
block_rect(xR, yR, nr, rtop = 1, cleft = 4, rh = 4, cwn = 3, border = col_blue, lwd = 3.5)
text(xR + 3, yL + 5.5, "2019 — COMPLETE block design", font = 2, cex = 1.05, col = col_blue)
text(xR + 3, yL + 4.95, "two BIG blocks (Reps) — each holds EVERY line once",
     cex = 0.82, col = col_grey)
text(xR + 1.5, yL + 4.4, "Rep I", font = 2, cex = 0.92, col = col_blue)
text(xR + 4.5, yL + 4.4, "Rep II", font = 2, cex = 0.92, col = col_blue)
# same line L1 appears once per rep
hl_cell(xR, yR, nr, 2, 1); hl_cell(xR, yR, nr, 3, 5)
text(xR + 3, yL - 0.55, "→ here L1 sits once in Rep I and once in Rep II",
     col = col_accent, cex = 0.72)

# ---- shared legend (bottom) ----
ly <- -1.15
sw <- function(x, fill, lab, tcol = col_dark) {
  rect(x, ly, x + 0.5, ly + 0.5, col = fill, border = "white", lwd = 1.4)
  text(x + 0.7, ly + 0.25, lab, adj = 0, cex = 0.82, col = tcol)
}
sw(0.6,  col_line,   "a breeding LINE (one plot)")
sw(5.0,  col_accent, "Eclipse")
sw(7.4,  col_blue,   "Zorro")
sw(9.4,  col_gold,   "Zenith")
text(11.2, ly + 0.25, "( red box = the SAME line, grown twice = replication )",
     adj = 0, cex = 0.78, col = col_accent)
dev.off()
cat("wrote 31_block_designs.png\n")

# ============================================================================
# FIG 32 — one PLOT (anatomy) + how plot -> block -> checks -> model connect
# ============================================================================
png("figures/32_plot_and_pipeline.png", width = 1850, height = 1000, res = 150)
layout(matrix(c(1, 2, 2), nrow = 1))

# ---- PANEL A: anatomy of one 4-row plot ----
par(mar = c(2, 2, 5, 1))
plot(NA, xlim = c(0, 10), ylim = c(0, 12), axes = FALSE, xlab = "", ylab = "")
title(main = "One PLOT = 4 rows\n(harvest the centre 2)", cex.main = 1.15, font.main = 2)
rowx <- c(2.2, 4.0, 5.8, 7.6); ytop <- 9.5; ybot <- 2.5
roles <- c("guard", "HARVEST", "HARVEST", "guard")
fills <- c("#d9d9d9", col_dark, col_dark, "#d9d9d9")
for (i in 1:4) {
  rect(rowx[i] - 0.55, ybot, rowx[i] + 0.55, ytop,
       col = adjustcolor(fills[i], 0.92), border = "white", lwd = 2)
  text(rowx[i], (ytop + ybot) / 2, roles[i], srt = 90,
       col = if (i %in% c(2, 3)) "white" else col_grey, cex = 0.82, font = 2)
}
# row spacing brace
arrows(rowx[2], ytop + 0.7, rowx[3], ytop + 0.7, code = 3, length = 0.07,
       col = col_accent, lwd = 1.8)
text((rowx[2] + rowx[3]) / 2, ytop + 1.25, "50 cm", col = col_accent, cex = 0.85)
# length brace
arrows(1.1, ybot, 1.1, ytop, code = 3, length = 0.07, col = col_accent, lwd = 1.8)
text(0.7, (ytop + ybot) / 2, "4.5 m (trimmed)", srt = 90, col = col_accent, cex = 0.85)
text(5, 1.3, "Only the centre 2 rows are cut\n→ one yield number for this plot",
     col = col_dark, cex = 0.92, font = 2)
text(5, 11.3, "( grey outer rows = guard / border effects )", col = col_grey, cex = 0.78)

# ---- PANEL B: how it all connects (top-down pipeline) ----
par(mar = c(1, 1, 5, 1))
plot(NA, xlim = c(0, 10), ylim = c(0, 14.4), axes = FALSE, xlab = "", ylab = "")
title(main = "How it all connects: from a plot in the field to one fair number per line",
      cex.main = 1.12, font.main = 2)
flowbox <- function(yc, txt, fill, border, tcol = col_dark, cex = 0.86, font = 1,
                    w = 9.2, h = 1.7) {
  rect(5 - w / 2, yc - h / 2, 5 + w / 2, yc + h / 2, col = fill, border = border, lwd = 2.2)
  text(5, yc, txt, col = tcol, cex = cex, font = font)
}
ys <- c(13.3, 11.0, 8.7, 6.4, 4.1, 1.8)
flowbox(ys[1], "FIELD  —  a grid of PLOTS  (rows R × passes P)",
        adjustcolor(col_grey, 0.16), col_grey)
flowbox(ys[2], "PLOTS are grouped into BLOCKS\nINCOMPLETE = small subsets (2018)      COMPLETE = full reps (2019)",
        adjustcolor(col_dark, 0.12), col_dark)
flowbox(ys[3], "CHECKS  E · Z · N  repeat across every block\n= rulers: the same genotype in many spots reveals the field's gradient",
        adjustcolor(col_gold, 0.22), col_gold, tcol = col_dark)
flowbox(ys[4], "Each PLOT: harvest centre 2 rows → one RAW plot value\n(part genetics, part lucky/unlucky spot)",
        adjustcolor(col_blue, 0.12), col_blue)
flowbox(ys[5], "MIXED MODEL (SpATS): subtract block + row/col + 2-D spatial trend\nusing replication + checks to separate field from genetics",
        adjustcolor(col_blue, 0.20), col_blue, tcol = col_dark)
flowbox(ys[6], "BLUP — one fair GENETIC value per line\n→ the clean phenotype every later lesson predicts",
        adjustcolor(col_dark, 0.85), col_dark, tcol = "white", font = 2)
for (i in 1:5)
  arrows(5, ys[i] - 0.9, 5, ys[i + 1] + 0.9, length = 0.11, lwd = 2.4, col = col_grey)
dev.off()
cat("wrote 32_plot_and_pipeline.png\n")

cat("\nExperimental-design figures written to figures/31..32\n")
