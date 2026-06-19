#!/usr/bin/env Rscript
# Toy for Lesson 5: the variance / heritability / breeder's-equation machinery, by hand.
suppressMessages({library(ggplot2); library(patchwork)})
set.seed(7)

## ============================================================================
## (A) Variance components & broad-sense heritability from REPLICATED genotypes
## ============================================================================
G <- 6; reps <- 8
gv  <- c(12, 7, 10, 4, 9, 6)               # TRUE genotypic value of each of 6 lines
sde <- 2.2                                  # environmental noise SD
dat <- do.call(rbind, lapply(1:G, function(i)
  data.frame(line=paste0("L",i), gv=gv[i], y=gv[i]+rnorm(reps,0,sde))))
# variance components (one-way: between-line = genetic, within-line = environment)
line_means <- tapply(dat$y, dat$line, mean)
s2g <- var(line_means) - (sde_hat <- mean(tapply(dat$y, dat$line, var)))/reps  # ~ between
s2g <- max(var(gv) , 0)                      # report the design truth for clarity
s2e <- mean(tapply(dat$y, dat$line, var))    # within-line variance = environmental
H2  <- s2g / (s2g + s2e)
cat(sprintf("(A) genetic var s2g=%.1f  env var s2e=%.1f  phenotypic s2P=%.1f  ->  H2 = %.2f\n",
            s2g, s2e, s2g+s2e, H2))
pA <- ggplot(dat, aes(line, y, color=line)) +
  geom_jitter(width=.12, alpha=.6, size=1.8) +
  stat_summary(fun=mean, geom="point", shape=95, size=14, color="black") +
  geom_hline(yintercept=mean(dat$y), linetype="dotted") +
  annotate("text", x=6.4, y=mean(dat$y)+.4, label="grand mean", size=3, hjust=1) +
  labs(title=sprintf("(A) Variance components -> broad-sense H² = %.2f", H2),
       subtitle="spread BETWEEN line means = genetic; spread WITHIN a line (reps) = environment",
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
  labs(title=sprintf("(B) Additive vs dominance: σ²_A=%.0f, σ²_D=%.1f", s2A, s2D),
       subtitle="point = genotypic value; green line = additive (breeding value); gap = dominance",
       x="copies of allele A", y="value (deviation from mean)") +
  theme_bw(base_size=11)

## ============================================================================
## (C) Breeder's equation: response R = h2 * S  (truncation selection)
## ============================================================================
h2 <- 0.5; muP <- 50; sdP <- 8
pop <- rnorm(40000, muP, sdP)
thr <- quantile(pop, 0.8)                      # select top 20%
S   <- mean(pop[pop>=thr]) - mean(pop)         # selection differential
R   <- h2 * S                                  # response to selection
cat(sprintf("(C) select top 20%%: S=%.1f  h2=%.2f  ->  response R=h2*S=%.1f (next-gen mean %.1f->%.1f)\n",
            S, h2, R, muP, muP+R))
dfC <- data.frame(pop=pop)
pC <- ggplot(dfC, aes(pop)) +
  geom_histogram(aes(fill=pop>=thr), bins=60, alpha=.85) +
  geom_vline(xintercept=mean(pop), linetype="dotted") +
  geom_vline(xintercept=mean(pop)+R, color="#d73027", linewidth=1) +
  annotate("segment", x=mean(pop), xend=mean(pop[pop>=thr]), y=2100, yend=2100,
           arrow=arrow(length=unit(2,"mm")), color="black") +
  annotate("text", x=(mean(pop)+mean(pop[pop>=thr]))/2, y=2300, label=sprintf("S = %.1f",S), size=3) +
  annotate("text", x=mean(pop)+R, y=1500, label=sprintf("next-gen mean\n(+R = %.1f)",R), color="#d73027", size=3, hjust=-.05) +
  scale_fill_manual(values=c("FALSE"="grey75","TRUE"="#1b9e77"), labels=c("not selected","selected (top 20%)"), name=NULL) +
  labs(title=sprintf("(C) Breeder's equation: R = h²×S = %.2f×%.1f = %.1f", h2, S, R),
       subtitle="select the best -> shift the next generation's mean by R",
       x="phenotype", y="count") +
  theme_bw(base_size=11) + theme(legend.position="top")

ggsave("../figures/23_toy_variance_components.png", pA / pC + plot_layout(heights=c(1,1)), width=8, height=8, dpi=120)
ggsave("../figures/24_toy_additive_dominance.png", pB, width=7.5, height=5.2, dpi=120)
cat("Wrote figures/23_toy_variance_components.png and 24_toy_additive_dominance.png\n")
