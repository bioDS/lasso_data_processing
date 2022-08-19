#!/usr/bin/env Rscript
library(Pint)
library(pROC)
library(dplyr)
library(foreach)
library(doMC)
library(ggplot2)
library(reshape2)

registerDoMC(cores=detectCores())

source("regenerate_data.R")
p100_dir <- "~/work/data/simulated_rerun/simulated_small_data_sample/"
ensure_p100_set_exists(p100_dir)
dir_8k <- "~/work/data/simulated_rerun/8k_only/"
ensure_8k_set_exists(dir_8k)

source("summary_functions.R")

get_all_effects_for_dir <- function(use_data, depth, num_nz, approximate_hierarchy = TRUE) {
    all_files <- list.files(use_data, pattern = "*.rds")
    all_effects <- foreach(file = all_files, .combine = 'c') %do% {
        list(get_file_effects(paste0(use_data, "/", file), depth, num_nz, approximate_hierarchy))
    }
    return(all_effects)
}

effects_dir <- dir_8k

effects_files <- list.files(effects_dir, pattern="*.rds")
num_true_effects <- foreach(file=effects_files, .combine='c') %do% {
    dat <- readRDS(paste0(effects_dir, file))
    return(nrow(dat$bi_ind) + nrow(dat$bij_ind))
}

hierarchy_time <- system.time(
    effects <- get_all_effects_for_dir(effects_dir, 2, 100, approximate_hierarchy = TRUE)
)
base_time <- system.time(
    effects_nohierarchy <- get_all_effects_for_dir(effects_dir, 2, 100, approximate_hierarchy = FALSE)
)

summarise_effects <- function(effects, num_true_effects) {
    summary <- foreach(effect = effects, true_total=num_true_effects, .combine='rbind') %do% {
        num_tp = sum(effect$TP)
        found_total = length(effect$TP)
        prec = num_tp/found_total
        rec = num_tp/true_total
        f1 <- 2*(prec * rec)/(prec + rec)
        df = data.frame(precision = prec, recall=rec, f1=f1)
        return(df)
    }
    return(summary)
}

remove_main <- function(eff) {
  return(eff %>% filter(type != "main"))
}

non_main_effects <- foreach(eff=effects) %do% {
  return(remove_main(eff))
}
non_main_effects_nohierarchy <- foreach(eff=effects_nohierarchy) %do% {
  return(remove_main(eff))
}

summary_with <- summarise_effects(non_main_effects, num_true_effects)
summary_nohierarchy <- summarise_effects(non_main_effects_nohierarchy, num_true_effects)


f1_comparison = melt(data.frame(Hierarchy=summary_with$f1, Plain=summary_nohierarchy$f1))
f1_plot <- ggplot(f1_comparison, aes(x = variable, y = value)) +
    geom_boxplot() +
    theme_bw() +
    xlab("Method") +
    ylab("F1 score") +
    expand_limits(y = 0)
ggsave("plots/hierarchy_f1.pdf", width = 3, height = 3)

precision_comparison = melt(data.frame(Hierarchy=summary_with$precision, Plain=summary_nohierarchy$precision))
precision_plot <- ggplot(precision_comparison, aes(x = variable, y = value)) +
    geom_boxplot() +
    theme_bw() +
    xlab("Method") +
    ylab("Precision") +
    expand_limits(y = 0)
ggsave("plots/hierarchy_precision.pdf", width = 3, height = 3)

recall_comparison = melt(data.frame(Hierarchy=summary_with$recall, Plain=summary_nohierarchy$recall))
recall_plot <- ggplot(recall_comparison, aes(x = variable, y = value)) +
    geom_boxplot() +
    theme_bw() +
    xlab("Method") +
    ylab("Recall") +
    expand_limits(y = 0)
ggsave("plots/hierarchy_recall.pdf", width = 3, height = 3)
