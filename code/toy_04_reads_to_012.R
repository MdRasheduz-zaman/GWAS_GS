#!/usr/bin/env Rscript
# Toy for Lesson 4: how reads at ONE SNP become a 0/1/2 dosage, for 3 samples.
# REF=C, ALT=T. Visualize stacked reads -> allele fraction -> dosage.
suppressMessages({library(ggplot2)})
mk <- function(sample, nC, nT){
  data.frame(sample=sample,
             base=c(rep("C (ref)",nC), rep("T (alt)",nT)),
             read=seq_len(nC+nT))
}
d <- rbind(mk("Sample A",8,0), mk("Sample B",5,4), mk("Sample C",0,7))
# labels: allele fraction + dosage call
lab <- data.frame(
  sample=c("Sample A","Sample B","Sample C"),
  txt=c("alt frac = 0/8 = 0.00  ->  dosage 0  (C/C, hom-ref)",
        "alt frac = 4/9 = 0.44  ->  dosage 1  (C/T, het)",
        "alt frac = 7/7 = 1.00  ->  dosage 2  (T/T, hom-alt)"))
d$sample <- factor(d$sample, levels=c("Sample A","Sample B","Sample C"))
lab$sample <- factor(lab$sample, levels=levels(d$sample))
ggplot(d, aes(read, sample, fill=base)) +
  geom_tile(color="white", height=0.6) +
  geom_text(data=lab, aes(x=0.3, y=sample, label=txt), inherit.aes=FALSE,
            hjust=0, vjust=-1.7, size=3.6) +
  scale_fill_manual(values=c("C (ref)"="#4575b4","T (alt)"="#d73027"), name="read base") +
  scale_x_continuous(limits=c(0,9), breaks=1:9) +
  labs(title="How aligned reads at one SNP become a 0/1/2 dosage",
       subtitle="Each tile = one read covering this position (REF=C, ALT=T)",
       x="reads stacked at this position", y=NULL) +
  theme_minimal(base_size=12) + theme(panel.grid=element_blank())
ggsave("../figures/12_toy_reads_to_dosage.png", width=8.5, height=4.2, dpi=120)
cat("Wrote figures/12_toy_reads_to_dosage.png\n")
