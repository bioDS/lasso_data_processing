#!/usr/bin/env Rscript
library(Pint)
library(pROC)
library(dplyr)
library(foreach)
library(doMC)
library(ggplot2)

source("summary_functions.R")

registerDoMC(cores=detectCores())

bench_sets <- c("simulated_small_data_sample", "8k_only", "wide_only", "3way")

plot_rocs <- function(all_pint_smrys, all_glint_smrys, all_whinter_smrys, use_glint, use_type, include_FN, output_filename) {
    ## main only
    if (use_type != "all") {
        pint <- all_pint_smrys |> filter(type == use_type)
        if (use_glint) {
            glint <- all_glint_smrys |> filter(type == use_type)
        }
        whinter <- all_whinter_smrys |> filter(type == use_type)
    } else {
        pint <- all_pint_smrys
        if (use_glint) {
            glint <- all_glint_smrys
        }
        whinter <- all_whinter_smrys
    }

    if (include_FN == FALSE) {
      pint <- pint |> filter(found == TRUE)
      glint <- glint |> filter(found == TRUE)
      whinter <- whinter |> filter(found == TRUE)
    }

    pint_roc <- roc(response = pint$TP, predict = abs(pint$strength))
    whinter_roc <- roc(response = whinter$TP, predict = abs(whinter$strength))
    if (use_glint) {
        glint_roc <- roc(response = glint$TP, predict = abs(glint$strength))
    }

    if (use_glint) {
      roc_list = list(Pint = pint_roc,
                   Whinter = whinter_roc,
                   Glinternet = glint_roc)
      glint_annotation =  sprintf("Glinternet auc: %.2f", round(glint_roc$auc, digits = 2))
    } else {
      roc_list = list(Pint = pint_roc,
                   Whinter = whinter_roc)
      glint_annotation = ""
    }
    roc_plot <-
        ggroc(roc_list) +
        theme_bw() +
        geom_segment(x = -1, y = 0, xend = 0, yend = 1, color = "#4c72b0") +
        annotate("text", x = 0.2, y = 0.3, label = sprintf("Pint auc: %.2f", round(pint_roc$auc, digits = 2))) +
        annotate("text", x = 0.2, y = 0.2, label = sprintf("Whinter auc: %.2f", round(whinter_roc$auc, digits = 2))) +
        annotate("text", x = 0.2, y = 0.1, label = glint_annotation)
    ggsave(roc_plot, file = output_filename, width = 5, height = 3)
}

for (bench_set in bench_sets) {
    if (bench_set == "wide_only") {
        use_glint = FALSE
    } else {
        use_glint = TRUE
    }

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

    plot_dir <- "plots/rocs/"
    if (! dir.exists(plot_dir)) {
      dir.create(plot_dir, recursive=TRUE)
    }
    plot_rocs(all_pint_smrys, all_glint_smrys, all_whinter_smrys, use_glint, "all", TRUE, sprintf("plots/rocs/all_comparison_roc_%s.pdf", bench_set))
    plot_rocs(all_pint_smrys, all_glint_smrys, all_whinter_smrys, use_glint, "main", TRUE, sprintf("plots/rocs/main_comparison_roc_%s.pdf", bench_set))
    plot_rocs(all_pint_smrys, all_glint_smrys, all_whinter_smrys, use_glint, "interaction", TRUE, sprintf("plots/rocs/int_comparison_roc_%s.pdf", bench_set))

    plot_rocs(all_pint_smrys, all_glint_smrys, all_whinter_smrys, use_glint, "all", FALSE, sprintf("plots/rocs/all_comparison_roc_%s_FoundOnly.pdf", bench_set))
    plot_rocs(all_pint_smrys, all_glint_smrys, all_whinter_smrys, use_glint, "main", FALSE, sprintf("plots/rocs/main_comparison_roc_%s_FoundOnly.pdf", bench_set))
    plot_rocs(all_pint_smrys, all_glint_smrys, all_whinter_smrys, use_glint, "interaction", FALSE, sprintf("plots/rocs/int_comparison_roc_%s_FoundOnly.pdf", bench_set))
}
