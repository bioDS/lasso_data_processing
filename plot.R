#!/usr/bin/env Rscript
library(ggplot2)
library(reshape2)
library(dplyr)


files <- c(
  # "1k_whinter_working_data-bench_summary.rds",
  "simulated_small_data_sample-bench_summary.rds",
  # "simulated_8k-bench_summary.rds"
  "8k_only-bench_summary.rds",
  #          "removed-bench_summary.rds",
  "wide_only-bench_summary.rds"
)

output_dir <- "bench_plots"

if (!dir.exists(output_dir)) {
  dir.create(output_dir)
}

for (file in files) {
  print(sprintf("reading file %s", file))
  prefix <- strsplit(file, split = ".", fixed = TRUE)[[1]][1]
  prefix <- paste(output_dir, prefix, sep = "/")
  dat <- readRDS(file)
  if (length(dat$glint_times) > 0) {
    print(paste("using glinternet for set", file, sep = " "))
    times <- melt(data.frame(dat$pint_times, dat$whinter_times, dat$glint_times))
    # precs = melt(data.frame(dat$pint_prec, dat$whinter_prec, dat$glint_prec, dat$combo_prec))
    # recs = melt(data.frame(dat$pint_rec, dat$whinter_rec, dat$glint_rec, dat$combo_rec))
    precs <- melt(data.frame(dat$pint_prec, dat$whinter_prec, dat$glint_prec))
    recs <- melt(data.frame(dat$pint_rec, dat$whinter_rec, dat$glint_rec))
  } else {
    print(paste("skipping glinternet for set", file, sep = " "))
    times <- melt(data.frame(dat$pint_times, dat$whinter_times))
    # precs = melt(data.frame(dat$pint_prec, dat$whinter_prec, dat$combo_prec))
    # recs = melt(data.frame(dat$pint_rec, dat$whinter_rec, dat$combo_rec))
    precs <- melt(data.frame(dat$pint_prec, dat$whinter_prec))
    recs <- melt(data.frame(dat$pint_rec, dat$whinter_rec))
  }
  fp <- melt(data.frame(dat$pint_fp, dat$whinter_fp))

  print(sprintf("median time for %s, %f", "pint", median(dat$pint_times)))
  print(sprintf("median time for %s, %f", "whinter", median(dat$whinter_times)))
  print(sprintf("median time for %s, %f", "glint", median(dat$glint_times)))

  names(times) <- c("Method", "Time")
  names(precs) <- c("Method", "Precision")
  names(recs) <- c("Method", "Recall")
  names(fp) <- c("Method", "FP")

  times <- times %>% mutate(
    Method =
      ifelse(Method == "dat.pint_times", "Pint",
        ifelse(Method == "dat.whinter_times", "Whinter",
          ifelse(Method == "dat.glint_times", "Glinternet", NA)
        )
      )
  )
  precs <- precs %>% mutate(
    Method =
      ifelse(Method == "dat.pint_prec", "Pint",
        ifelse(Method == "dat.whinter_prec", "Whinter",
          ifelse(Method == "dat.combo_prec", "Combo",
            ifelse(Method == "dat.glint_prec", "Glinternet", NA)
          )
        )
      )
  )
  recs <- recs %>% mutate(
    Method =
      ifelse(Method == "dat.pint_rec", "Pint",
        ifelse(Method == "dat.whinter_rec", "Whinter",
          ifelse(Method == "dat.combo_rec", "Combo",
            ifelse(Method == "dat.glint_rec", "Glinternet", NA)
          )
        )
      )
  )
  fp <- fp %>% mutate(
    Method =
      ifelse(Method == "dat.pint_fp", "Pint",
        ifelse(Method == "dat.whinter_fp", "Whinter",
          ifelse(Method == "dat.combo_fp", "Combo",
            ifelse(Method == "dat.glint_fp", "Glinternet", NA)
          )
        )
      )
  )

  f1 <- data.frame(Method = precs$Method, F1 = 2 * precs$Precision * recs$Recall / (precs$Precision + recs$Recall))

  time_plot <- ggplot(times, aes(x = Method, y = Time)) +
    geom_boxplot() +
    theme_bw() +
    ylab("Time (s)") +
    scale_y_continuous(trans = "log2") +
    expand_limits(y = 0)
  prec_plot <- ggplot(precs, aes(x = Method, y = Precision)) +
    geom_boxplot() +
    theme_bw() +
    expand_limits(y = 0)
  rec_plot <- ggplot(recs, aes(x = Method, y = Recall)) +
    geom_boxplot() +
    theme_bw() +
    expand_limits(y = 0)
  f1_plot <- ggplot(f1, aes(x = Method, y = F1)) +
    geom_boxplot() +
    theme_bw() +
    expand_limits(y = 0)
  fp_plot <- ggplot(fp, aes(x = Method, y = FP)) +
    geom_boxplot() +
    theme_bw() +
    ylab("False Positives") +
    expand_limits(y = 0)

  print(sprintf("saving plots in %s-*", prefix))

  ggsave(time_plot, file = paste(prefix, "bench_time_plot.pdf", sep = "-"), width = 2, height = 2)
  ggsave(prec_plot, file = paste(prefix, "bench_prec_plot.pdf", sep = "-"), width = 2, height = 2)
  ggsave(rec_plot, file = paste(prefix, "bench_rec_plot.pdf", sep = "-"), width = 2, height = 2)
  ggsave(f1_plot, file = paste(prefix, "bench_f1_plot.pdf", sep = "-"), width = 2, height = 2)
  ggsave(fp_plot, file = paste(prefix, "bench_fp_plot.pdf", sep = "-"), width = 2, height = 2)

  # additional effects found with lasso vs. whinter
  combo_advantage <- ggplot(data.frame(dat$combo_advantage), aes(x = dat.combo_advantage)) +
    geom_histogram(binwidth = 1) +
    xlab("Additionally Identified Interactions")
  ggsave(combo_advantage, file = paste(prefix, "combo_advantage_plot.pdf", sep = "-"), width = 4, height = 3)
}
