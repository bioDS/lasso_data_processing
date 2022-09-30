#!/usr/bin/env Rscript
library(Pint)
library(pROC)
library(dplyr)
library(foreach)
library(doMC)
library(ggplot2)

source("summary_functions.R")

registerDoMC(cores=detectCores())

bench_sets <- c("simulated_small_data_sample", "8k_only", "wide_only")

for (bench_set in bench_sets) {
    dir <- paste("whinter_pint_comparison", bench_set, sep = "/")

    rds_files <- list.files(dir, pattern = "all.rds", recursive = TRUE, full.names = TRUE)

    output <- c()
    for (file in rds_files) {
        output <- cbind(output, readRDS(file))
    }

    output <- data.frame(output)

    all_pint_smrys <- foreach(set = output, .combine = rbind) %do% { set$pint_smry }
    all_whinter_smrys <- foreach(set = output, .combine = rbind) %do% { set$whinter_smry }
    all_glint_smrys <- foreach(set = output, .combine = rbind) %do% { set$glint_smry }

    # Including False Negatives
    pint_roc <- roc(response = all_pint_smrys$TP, predict = abs(all_pint_smrys$strength))
    whinter_roc <- roc(response = all_whinter_smrys$TP, predict = abs(all_whinter_smrys$strength))
    glint_roc <- roc(response = all_glint_smrys$TP, predict = abs(all_glint_smrys$strength))


    roc_plot <-
        ggroc(list(Pint = pint_roc,
                   Whinter = whinter_roc,
                   Glinternet = glint_roc)) +
        theme_bw() +
        geom_segment(x = -1, y = 0, xend = 0, yend = 1, color = "#4c72b0") #+
        ## annotate("text", x = 0.2, y = 0.1, label = sprintf("auc: %.2f", round(roc$auc, digits = 2)))
    ggsave(roc_plot, file = "plots/overall_comparison_roc.pdf", width = 5, height = 3)

    ## main only
    pint_main <- all_pint_smrys |> filter(type == "main")
    glint_main <- all_glint_smrys |> filter(type == "main")
    whinter_main <- all_whinter_smrys |> filter(type == "main")
    pint_main_roc <- roc(response = pint_main$TP, predict = abs(pint_main$strength))
    whinter_main_roc <- roc(response = whinter_main$TP, predict = abs(whinter_main$strength))
    glint_main_roc <- roc(response = glint_main$TP, predict = abs(glint_main$strength))

    roc_plot <-
        ggroc(list(Pint = pint_roc,
                   Whinter = whinter_roc,
                   Glinternet = glint_roc)) +
        theme_bw() +
        geom_segment(x = -1, y = 0, xend = 0, yend = 1, color = "#4c72b0") #+
        ## annotate("text", x = 0.2, y = 0.1, label = sprintf("auc: %.2f", round(roc$auc, digits = 2)))
    ggsave(roc_plot, file = "plots/main_comparison_roc.pdf", width = 5, height = 3)

    ## interaction only
    pint_interaction <- all_pint_smrys |> filter(type == "interaction")
    glint_interaction <- all_glint_smrys |> filter(type == "interaction")
    whinter_interaction <- all_whinter_smrys |> filter(type == "interaction")
    pint_interaction_roc <- roc(response = pint_interaction$TP, predict = abs(pint_interaction$strength))
    whinter_interaction_roc <- roc(response = whinter_interaction$TP, predict = abs(whinter_interaction$strength))
    glint_interaction_roc <- roc(response = glint_interaction$TP, predict = abs(glint_interaction$strength))

    roc_plot <-
        ggroc(list(Pint = pint_interaction_roc,
                   Whinter = whinter_interaction_roc,
                   Glinternet = glint_interaction_roc)) +
        theme_bw() +
        geom_segment(x = -1, y = 0, xend = 0, yend = 1, color = "#4c72b0") #+
        ## annotate("text", x = 0.2, y = 0.1, label = sprintf("auc: %.2f", round(roc$auc, digits = 2)))
    ggsave(roc_plot, file = "plots/interaction_comparison_roc.pdf", width = 5, height = 3)
}
