#!/usr/bin/env Rscript
library(foreach)
library(doMC)
library(ggplot2)
library(reshape2)
library(doParallel)
library(Pint)

source("3way_check_functions.R")

# fn = "simulated_data/n1000_p100_SNR5_nbi20_nbij30_nbijk50_nlethals0_12086.rds"
# results = run_23way(fn)
# fn2 = "simulated_data/n1000_p100_SNR10_nbi10_nbij100_nbijk1000_nlethals0_85437.rds"
# results2 = run_23way(fn)
#
# big2_pr = rbind(results$big2_pr, results2$big2_pr)
# big3_pr = rbind(results$big3_pr, results2$big3_pr)

## data_dir <- "3way_data_to_run"
data_dir <- "./data/simulated_rerun/3way/"

big2_pr <- data.frame()
big3_pr <- data.frame()

rds_files <- list.files(data_dir, pattern = "*.rds")

# testfile <- "simulated_data/n2000_p200_SNR2_nbi20_nbij1005_nbijk13333_nlethals0_11829.rds"
# run_23way(testfile)

# registerDoMC(cores=detectCores())

cl <- makeCluster(8)
registerDoParallel(cl)
all_results <- foreach(file = rds_files, .combine = rbind) %do% {
  library(Pint)
  print(sprintf("running %s", file))
  run_23way(sprintf("%s/%s", data_dir, file))
  # big2_pr = rbind(big2_pr, results$big2_pr)
  # big3_pr = rbind(big3_pr, results$big3_pr)
}

big2_pr <- bind_rows(all_results[, "big2_pr"])
big3_pr <- bind_rows(all_results[, "big3_pr"])
all2_pr <- bind_rows(all_results[, "all2_pr"])
all3_pr <- bind_rows(all_results[, "all3_pr"])

big_23 <- list(big2_pr, big3_pr)
all_23 <- list(all2_pr, all3_pr)

plot_pr <- function(use_23, plot_dir) {
  dir.create(plot_dir)
  use2_pr <- use_23[[1]]
  use3_pr <- use_23[[2]]

  main_prec <- data.frame(use2_pr$prec_main, use3_pr$prec_main)
  names(main_prec) <- c("2-way", "3-way")
  main_prec <- main_prec %>% melt()
  names(main_prec) <- c("depth", "Precision")
  main_rec <- data.frame(use2_pr$rec_main, use3_pr$rec_main)
  names(main_rec) <- c("2-way", "3-way")
  main_rec <- main_rec %>% melt()
  names(main_rec) <- c("depth", "Recall")
  main_comp <- cbind(main_prec, main_rec$Recall)
  names(main_comp) <- c("depth", "Precision", "Recall")

  pair_prec <- data.frame(use2_pr$prec_pair, use3_pr$prec_pair)
  names(pair_prec) <- c("2-way", "3-way")
  pair_prec <- pair_prec %>% melt()
  names(pair_prec) <- c("depth", "Precision")
  pair_rec <- data.frame(use2_pr$rec_pair, use3_pr$rec_pair)
  names(pair_rec) <- c("2-way", "3-way")
  pair_rec <- pair_rec %>% melt()
  names(pair_rec) <- c("depth", "Recall")
  pair_comp <- cbind(pair_prec, pair_rec$Recall)
  names(pair_comp) <- c("depth", "Precision", "Recall")

  trip_prec <- data.frame(use3_pr$prec_trip)
  names(trip_prec) <- c("3-way")
  trip_prec <- trip_prec %>% melt()
  names(trip_prec) <- c("depth", "Precision")
  trip_rec <- data.frame(use3_pr$rec_trip)
  names(trip_rec) <- c("3-way")
  trip_rec <- trip_rec %>% melt()
  names(trip_rec) <- c("depth", "Recall")
  trip_comp <- cbind(trip_prec, trip_rec$Recall)
  names(trip_comp) <- c("depth", "Precision", "Recall")

  overall_prec <- data.frame(use2_pr$prec_overall, use3_pr$prec_overall)
  names(overall_prec) <- c("2-way", "3-way")
  overall_prec <- overall_prec %>% melt()
  names(overall_prec) <- c("depth", "Precision")
  overall_rec <- data.frame(use2_pr$rec_overall, use3_pr$rec_overall)
  names(overall_rec) <- c("2-way", "3-way")
  overall_rec <- overall_rec %>% melt()
  names(overall_rec) <- c("depth", "Recall")
  overall_comp <- cbind(overall_prec, overall_rec$Recall)
  names(overall_comp) <- c("depth", "Precision", "Recall")
  # old (no recall) version
  #pair_comp <- data.frame(all2_pr$prec_pair, all3_pr$prec_pair)
  #names(pair_comp) <- c("2-way", "3-way")
  #pair_comp <- pair_comp %>% melt()
  #names(pair_comp) <- c("depth", "Precision")
  #
  #triple_comp <- data.frame(all3_pr$prec_trip)
  #names(triple_comp) <- c("3-way")
  #triple_comp <- triple_comp %>% melt()
  #names(triple_comp) <- c("depth", "Precision")

  ## pairwise_comp <- data.frame(all2_pr$prec_pair, all3_pr$prec_pair) %>% melt()
  ## triple_comp <- data.frame(all2_pr$prec_trip, all3_pr$prec_trip) %>% melt()

  # precision
  ggplot(main_comp, aes(x = depth, y = Precision)) +
    geom_boxplot() +
    xlab("Interaction Depth") +
    theme_bw()
  ggsave(sprintf("%s/main_comp_precision.pdf", plot_dir), height = 2, width = 2)
  ggplot(pair_comp, aes(x = depth, y = Precision)) +
    geom_boxplot() +
    xlab("Interaction Depth") +
    theme_bw()
  ggsave(sprintf("%s/pairwise_comp_precision.pdf", plot_dir), height = 2, width = 2)
  ggplot(trip_comp, aes(x = depth, y = Precision)) +
    geom_boxplot() +
    xlab("Interaction Depth") +
    theme_bw()
  ggsave(sprintf("%s/trip_comp_precision.pdf", plot_dir), height = 2, width = 1.8)
  ggplot(overall_comp, aes(x = depth, y = Precision)) +
    geom_boxplot() +
    xlab("Interaction Depth") +
    theme_bw()
  ggsave(sprintf("%s/overall_comp_precision.pdf", plot_dir), height = 2, width = 2)

  # recall
  ggplot(main_comp, aes(x = depth, y = Recall)) +
    geom_boxplot() +
    xlab("Interaction Depth") +
    theme_bw()
  ggsave(sprintf("%s/main_comp_recall.pdf", plot_dir), height = 2, width = 2)
  ggplot(pair_comp, aes(x = depth, y = Recall)) +
    geom_boxplot() +
    xlab("Interaction Depth") +
    theme_bw()
  ggsave(sprintf("%s/pairwise_comp_recall.pdf", plot_dir), height = 2, width = 2)
  ggplot(trip_comp, aes(x = depth, y = Recall)) +
    geom_boxplot() +
    xlab("Interaction Depth") +
    theme_bw()
  ggsave(sprintf("%s/trip_comp_recall.pdf", plot_dir), height = 2, width = 1.8)
  ggplot(overall_comp, aes(x = depth, y = Recall)) +
    geom_boxplot() +
    xlab("Interaction Depth") +
    theme_bw()
  ggsave(sprintf("%s/overall_comp_recall.pdf", plot_dir), height = 2, width = 2)

  # overall recall


  print("three-way fraction of equivalent-TP:")
  threeway_predictions <- bind_rows(all_results[, "depth3_results"]) %>% filter(!is.na(gene_k))
  print("TP:")
  threeway_predictions$TP %>% summary
  print("equiv_tp:")
  threeway_predictions$equiv_tp %>% summary

  print("overall TP, pairwise vs triple:")
  all2 <- bind_rows(all_results[, "depth2_results"])
  all3 <- bind_rows(all_results[, "depth3_results"])
  summary(all2$TP)
  summary(all3$TP)

  print("> 1 s.d. effects summary")

  all_big3_results <- bind_rows(all_results[, "big3_results"])
  all_big2_results <- bind_rows(all_results[, "big2_results"])

  ## all_big3_results %>%
  ##   filter(type == "three-way") %>%
  ##   select(TP) %>%
  ##   summary()
  ## all_big3_results %>%
  ##   select(TP) %>%
  ##   summary()
  ## all_big2_results %>%
  ##   select(TP) %>%
  ##   summary()


  ## all_big3_results %>%
  ##   filter(type != "three-way") %>%
  ##   select(TP) %>%
  ##   summary()
  ## all_big2_results %>%
  ##   filter(type != "three-way") %>%
  ##   select(TP) %>%
  ##   summary()

  print("2-way tp vs. equivalents:")
  print("TP")
  print(summary(all_big2_results$TP))
  print("equiv_tp")
  print(summary(all_big2_results$equiv_tp))

  print("3-way tp vs. equivalents:")
  print("TP")
  print(summary(all_big3_results$TP))
  print("equiv_tp")
  print(summary(all_big3_results$equiv_tp))

  big2_main <- (all_big2_results %>% filter(type=="main")) %>% select(gene_i, gene_j, TP, equiv_tp)
  big3_main <- (all_big3_results %>% filter(type=="main")) %>% select(gene_i, gene_j, TP, equiv_tp)
  common_main <- intersect(big2_main, big3_main)
  print("difference between big2_main and common_main:")
  print(setdiff(big2_main, common_main))
  print("difference between big3_main and common_main:")
  print(setdiff(big3_main, common_main))
  print("common main:")
  print(nrow(common_main))
  print(summary(common_main$TP))

  big2_pairwise <- (all_big2_results %>% filter(type=="pairwise")) %>% select(gene_i, gene_j, TP, equiv_tp)
  big3_pairwise <- (all_big3_results %>% filter(type=="pairwise")) %>% select(gene_i, gene_j, TP, equiv_tp)
  common_pairwise <- intersect(big2_pairwise, big3_pairwise)
  print("difference between big2_pairwise and common_pairwise:")
  print(setdiff(big2_pairwise, common_pairwise))
  print("difference between big3_pairwise and common_pairwise:")
  print(setdiff(big3_pairwise, common_pairwise))
  print("common pairwise:")
  print(nrow(common_pairwise))
  print(summary(common_pairwise$TP))


  big2_triple <- (all_big2_results %>% filter(type=="three-way")) %>% select(gene_i, gene_j, gene_k, TP, equiv_tp)
  big3_triple <- (all_big3_results %>% filter(type=="three-way")) %>% select(gene_i, gene_j, gene_k, TP, equiv_tp)
  common_triple <- intersect(big2_triple, big3_triple)
  print("difference between big2_triple and common_triple:")
  print(setdiff(big2_triple, common_triple))
  print("difference between big3_triple and common_triple:")
  print(setdiff(big3_triple, common_triple))
  print("common triple:")
  print(nrow(common_triple))
  print(summary(common_triple$TP))
}


plot_pr(big_23, "./plots/big/")
plot_pr(all_23, "./plots/all/")
