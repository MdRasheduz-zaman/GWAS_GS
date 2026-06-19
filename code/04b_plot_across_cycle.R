#!/usr/bin/env Rscript
# Plots the across-cycle results captured from 04_across_cycle.R (20 reps/point).
# These are the REAL reproduced values (yield 2019, ST GBLUP vs MT yield+seed-weight).
suppressMessages(library(ggplot2))
out <- data.frame(
  prop  = rep(c(0,10,20,30,40), each = 2),
  model = rep(c("ST (yield only)", "MT (yield + seed wt)"), 5),
  acc   = c(0.253,0.308, 0.305,0.360, 0.334,0.380, 0.372,0.419, 0.376,0.443))
ggplot(out, aes(prop, acc, color = model)) +
  geom_line(linewidth = 1) + geom_point(size = 2.5) +
  scale_color_manual(values = c("ST (yield only)" = "#7570B3",
                                "MT (yield + seed wt)" = "#1B9E77")) +
  scale_x_continuous(breaks = c(0,10,20,30,40)) +
  labs(title = "Predicting the NEW breeding cycle (yield 2019)",
       subtitle = "Accuracy rises as cycle-2 lines join training; MT (correlated trait) leads throughout",
       x = "% of cycle-2 lines added to the training set", y = "Prediction accuracy (r)", color = "") +
  theme_bw(base_size = 12) + theme(legend.position = "top")
ggsave("../figures/07_across_cycle.png", width = 8, height = 6, dpi = 120)
cat("DONE -> figures/07_across_cycle.png\n")
