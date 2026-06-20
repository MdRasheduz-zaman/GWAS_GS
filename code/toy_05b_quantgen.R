#!/usr/bin/env Rscript
# Toy for Lesson 5: the variance / heritability / breeder's-equation machinery, by hand.
suppressMessages({library(ggplot2); library(patchwork)})
set.seed(7)

## ============================================================================
## (A) Variance components & broad-sense heritability from REPLICATED genotypes
## ----------------------------------------------------------------------------
## We use DETERMINISTIC replicate offsets (mean 0) instead of random noise, so
## every number the reader sees is EXACT and self-consistent: the line means
## equal gv, the grand mean is exactly 8.0, and the variance components are the
## exact 8.4 / 4.3 / 12.7 quoted in the text.
## ============================================================================
G <- 6; reps <- 8
gv  <- c(12, 7, 10, 4, 9, 6)               # genotypic value of each of 6 lines (mean 8)
off <- c(-2.9, -2.1, -1.4, -0.6, 0.6, 1.4, 2.1, 2.9)   # rep offsets: mean 0, var ~ 4.3
dat <- do.call(rbind, lapply(1:G, function(i)
  data.frame(line=paste0("L",i), gv=gv[i], y=gv[i]+off)))
gmean   <- mean(gv)                          # grand mean = 8.0 (exact)
s2g     <- var(gv)                           # BETWEEN-line variance = genetic    = 8.4
s2e     <- var(off)                          # WITHIN-line variance  = environment ~ 4.3
s2P     <- s2g + s2e                         # phenotypic variance               = 12.7
H2      <- s2g / s2P                          # broad-sense heritability          = 0.66
sdg     <- sqrt(s2g); sde_val <- sqrt(s2e)   # the SDs we draw as brackets
m1      <- gv[1]                              # L1's mean, for the within-line bracket
cat(sprintf("(A) s2g=%.1f  s2e=%.1f  s2P=%.1f  ->  H2=%.2f  (grand mean=%.1f)\n",
            s2g, s2e, s2P, H2, gmean))
pA <- ggplot(dat, aes(line, y, color=line)) +
  geom_jitter(width=.12, alpha=.6, size=1.8) +
  stat_summary(fun=mean, geom="point", shape=95, size=14, color="black") +
  geom_hline(yintercept=gmean, linetype="dotted") +
  annotate("text", x=6.5, y=gmean-0.7, label=sprintf("grand mean = %.1f", gmean),
           size=3, hjust=1, color="grey25") +
  # WITHIN-line spread = sigma_e  (capped scale bar, drawn on L1)
  annotate("segment", x=1.34, xend=1.34, y=m1-sde_val, yend=m1+sde_val,
           color="#d95f02", linewidth=0.9, arrow=arrow(ends="both", length=unit(2,"mm"))) +
  annotate("segment", x=1.22, xend=1.46, y=m1-sde_val, yend=m1-sde_val, color="#d95f02", linewidth=0.7) +
  annotate("segment", x=1.22, xend=1.46, y=m1+sde_val, yend=m1+sde_val, color="#d95f02", linewidth=0.7) +
  annotate("text", x=1.52, y=m1, hjust=0, color="#d95f02", size=2.7, lineheight=.9,
           label=sprintf("σ_e = %.1f\n(within\na line)", sde_val)) +
  # BETWEEN-line spread = sigma_g  (capped scale bar in the far-right margin,
  # centred on the grand-mean line so it clearly measures ±σ_g from it)
  annotate("segment", x=7.45, xend=7.45, y=gmean-sdg, yend=gmean+sdg,
           color="#2166AC", linewidth=0.9, arrow=arrow(ends="both", length=unit(2,"mm"))) +
  annotate("segment", x=7.31, xend=7.59, y=gmean-sdg, yend=gmean-sdg, color="#2166AC", linewidth=0.7) +
  annotate("segment", x=7.31, xend=7.59, y=gmean+sdg, yend=gmean+sdg, color="#2166AC", linewidth=0.7) +
  annotate("text", x=7.66, y=gmean, hjust=0, color="#2166AC", size=2.7, lineheight=.9,
           label=sprintf("σ_g = %.1f\n(SD of the\nline means)", sdg)) +
  # the actual numbers, on the plot
  annotate("label", x=4.15, y=15.2, hjust=0.5, size=3, fill="white", label.size=0.4, lineheight=.95,
           label=sprintf("σ²_g = %.1f    σ²_e = %.1f\nσ²_P = σ²_g + σ²_e = %.1f\nH² = σ²_g / σ²_P = %.2f",
                         s2g, s2e, s2P, H2)) +
  scale_x_discrete(expand=expansion(add=c(0.6, 2.9))) +
  coord_cartesian(ylim=c(min(dat$y)-0.4, 16.3), clip="off") +
  labs(title="(A) Splitting the variance: where σ²_g, σ²_e and H² come from",
       subtitle="each dot = one replicate plot · black bar = line mean · BETWEEN lines = genetic, WITHIN a line = environment",
       x="genotype (line)", y="phenotype") +
  theme_bw(base_size=11) + theme(legend.position="none")

## ============================================================================
## (B) One locus: split genetic variance into ADDITIVE + DOMINANCE (Falconer)
## ============================================================================
a <- 10; d <- 4; p <- 0.6; q <- 1-p          # +a (AA), d (Aa), -a (aa); p=freq(A)
alpha <- a + d*(q-p)                          # avg effect of allele substitution
M  <- a*(p-q) + 2*p*q*d                        # population mean (deviation origin = midpoint)
Gval <- c(`0`=-a, `1`=d, `2`=a) - M            # genotypic values (centered)
BV   <- (c(0,1,2) - 2*p) * alpha               # breeding values (the additive line)
Dev  <- Gval - BV                              # dominance deviations
s2A  <- 2*p*q*alpha^2
s2D  <- (2*p*q*d)^2
cat(sprintf("(B) one locus: s2A=%.1f  s2D=%.1f  ->  share additive=%.0f%%  (h2 < H2 because of dominance)\n",
            s2A, s2D, 100*s2A/(s2A+s2D)))
freq <- c(q^2, 2*p*q, p^2)
dfB <- data.frame(x=0:2, Gval=Gval, BV=BV, freq=freq,
                  lab=c("aa","Aa","AA"))
pB <- ggplot(dfB, aes(x, Gval)) +
  geom_abline(intercept=-2*p*alpha, slope=alpha, color="#1b9e77", linewidth=1) +
  geom_segment(aes(xend=x, yend=BV), color="#d95f02", linewidth=1, linetype="dashed") +
  geom_point(aes(size=freq), color="#2166AC") +
  geom_text(aes(label=lab), vjust=-1.3, size=4) +
  annotate("text", x=1.55, y=alpha*1.55-2*p*alpha+0.6, label="additive line\n(slope = α)", color="#1b9e77", size=3, hjust=0) +
  annotate("text", x=1.02, y=mean(c(dfB$Gval[2],dfB$BV[2])), label="  dominance\n  deviation", color="#d95f02", size=3, hjust=0) +
  scale_size(range=c(3,9), guide="none") +
  scale_x_continuous(breaks=0:2, labels=c("0 (aa)","1 (Aa)","2 (AA)")) +
  labs(title=sprintf("Additive vs dominance: σ²_A=%.0f, σ²_D=%.1f", s2A, s2D),
       subtitle="point = genotypic value; green line = additive (breeding value); gap = dominance",
       x="copies of allele A", y="value (deviation from mean)") +
  theme_bw(base_size=11)

## ============================================================================
## (B-of-fig23) Breeder's equation R = h2*S, reusing the SAME toy as panel A
## ----------------------------------------------------------------------------
## No new numbers: the mean, the spread, and the heritability all come straight
## from panel A's variance components (s2g, s2e, H2). We just picture them across
## a big population and select the top 20%.
## ============================================================================
muG  <- gmean                       # 8.0  (panel A's grand mean)
h2br <- H2                          # 0.66 (panel A's heritability; toy ~ additive so h2 ≈ H2)
sdP  <- sqrt(s2P)                   # 3.6  (panel A's phenotypic SD)
set.seed(11)
pop <- rnorm(60000, muG, sdP)       # imagine panel A's mean & spread across a big population
thr <- quantile(pop, 0.8)           # selection threshold: keep the best 20%
selmean <- mean(pop[pop >= thr])    # mean of the selected group
S   <- selmean - muG                # selection differential  (how special the kept ones are)
R   <- h2br * S                     # response: only the heritable fraction is passed on
muN <- muG + R                      # next generation's mean
cat(sprintf("(B) mean=%.1f sel.mean=%.1f  S=%.1f  h2=%.2f  R=%.1f  next-gen %.1f->%.1f\n",
            muG, selmean, S, h2br, R, muG, muN))
ymx <- max(hist(pop, breaks=60, plot=FALSE)$counts)
pBr <- ggplot(data.frame(pop=pop), aes(pop)) +
  geom_histogram(aes(fill = pop >= thr), bins=60, alpha=.85) +
  # the three means, as vertical reference lines
  geom_vline(xintercept=muG,     linetype="dotted", color="grey25") +
  geom_vline(xintercept=selmean, linetype="dashed", color="#1b7837") +
  geom_vline(xintercept=muN,     color="#d73027", linewidth=1) +
  # S = selected mean - population mean  (full arrow, top)
  annotate("segment", x=muG, xend=selmean, y=0.97*ymx, yend=0.97*ymx,
           arrow=arrow(length=unit(2,"mm"), ends="both"), color="black") +
  annotate("text", x=(muG+selmean)/2, y=1.04*ymx, size=3,
           label=sprintf("S = %.1f - %.1f = %.1f", selmean, muG, S)) +
  # R = h2 * S  (shorter arrow below: lands only h2 of the way)
  annotate("segment", x=muG, xend=muN, y=0.73*ymx, yend=0.73*ymx,
           arrow=arrow(length=unit(2,"mm")), color="#d73027") +
  annotate("text", x=(muG+muN)/2, y=0.80*ymx, size=3, color="#d73027",
           label=sprintf("R = h²×S = %.2f×%.1f = %.1f", h2br, S, R)) +
  # name each line (rotated, along the line)
  annotate("text", x=muG,     y=0.50*ymx, angle=90, vjust=-0.4, size=2.7, color="grey25",
           label=sprintf("mean = %.1f", muG)) +
  annotate("text", x=selmean, y=0.34*ymx, angle=90, vjust=-0.4, size=2.7, color="#1b7837",
           label=sprintf("selected mean = %.1f", selmean)) +
  annotate("text", x=muN,     y=0.16*ymx, angle=90, vjust=1.3,  size=2.7, color="#d73027",
           label=sprintf("next-gen mean = %.1f", muN)) +
  scale_fill_manual(values=c("FALSE"="grey75","TRUE"="#1b9e77"),
                    labels=c("not selected","selected (top 20%)"), name=NULL) +
  coord_cartesian(ylim=c(0, 1.08*ymx), clip="off") +
  labs(title="(B) Breeder's equation R = h²×S: the next gen moves only h² of S",
       subtitle=sprintf("same population as panel A (mean %.1f, σ_P=%.1f); h² ≈ H² = %.2f (toy is ~additive, σ²_D≈0)",
                        muG, sdP, h2br),
       x="phenotype", y="count of individuals") +
  theme_bw(base_size=11) + theme(legend.position="top")

ggsave("../figures/23_toy_variance_components.png", pA / pBr + plot_layout(heights=c(1,1)), width=9.8, height=8.4, dpi=120)
ggsave("../figures/24_toy_additive_dominance.png", pB, width=7.5, height=5.2, dpi=120)
cat("Wrote figures/23_toy_variance_components.png and 24_toy_additive_dominance.png\n")
