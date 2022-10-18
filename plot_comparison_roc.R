#!/usr/bin/env Rscript
library(Pint)
library(pROC)
library(dplyr)
library(foreach)
library(doMC)
library(ggplot2)
library(reshape2)

source("summary_functions.R")

registerDoMC(cores=detectCores())

#bench_sets <- c("simulated_small_data_sample", "8k_only", "wide_only", "3way")
#bench_sets <- c("simulated_small_data_sample", "8k_only", "wide_only_10k")
#bench_sets <- c("wide_only_10k")
bench_sets <- c("3way")

plot_rocs <- function(all_pint_smrys, all_pint_hierarchy_smrys, all_pint_dedup_smrys,
                      all_pint_unbiased_smrys, all_glint_smrys, all_whinter_smrys,
                      use_glint, use_type, include_FN, output_filename, use_equiv_tp=FALSE) {
    ## main only
    gc()
    if (use_type != "all") {
        pint <- all_pint_smrys |> filter(type == use_type)
        pint_hierarchy <- all_pint_hierarchy_smrys |> filter(type == use_type)
        pint_dedup <- all_pint_dedup_smrys |> filter(type == use_type)
        pint_unbiased <- all_pint_unbiased_smrys |> filter(type == use_type)
        if (use_glint) {
            glint <- all_glint_smrys |> filter(type == use_type)
        }
        whinter <- all_whinter_smrys |> filter(type == use_type)
    } else {
        pint <- all_pint_smrys
        pint_hierarchy <- all_pint_hierarchy_smrys
        pint_dedup <- all_pint_dedup_smrys
        pint_unbiased <- all_pint_unbiased_smrys
        if (use_glint) {
            glint <- all_glint_smrys
        }
        whinter <- all_whinter_smrys
    }

    if (include_FN == FALSE) {
      pint <- pint |> filter(found == TRUE)
      pint_hierarchy <- pint_hierarchy |> filter(found == TRUE)
      pint_dedup <- pint_dedup |> filter(found == TRUE)
      pint_unbiased <- pint_unbiased |> filter(found == TRUE)
      if (use_glint) {
      glint <- glint |> filter(found == TRUE)
      }
      whinter <- whinter |> filter(found == TRUE)
    }

    if (use_equiv_tp) {
        pint_roc <- roc(response = pint$equiv_tp & pint$TP, predict = abs(pint$strength))
        gc()
        pint_hierarchy_roc <- roc(response = pint_hierarchy$equiv_tp & pint_hierarchy$TP, predict = abs(pint_hierarchy$strength))
        gc()
        pint_dedup_roc <- roc(response = pint_dedup$equiv_tp & pint_dedup$TP, predict = abs(pint_dedup$strength))
        gc()
        pint_unbiased_roc <- roc(response = pint_unbiased$equiv_tp & pint_unbiased$TP, predict = abs(pint_unbiased$strength))
        gc()
        whinter_roc <- roc(response = whinter$equiv_tp & whinter$TP, predict = abs(whinter$strength))
        if (use_glint) {
            glint_roc <- roc(response = glint$equiv_tp & glint$TP, predict = abs(glint$strength))
        }
    } else {
        pint_roc <- roc(response = pint$TP, predict = abs(pint$strength))
        gc()
        pint_hierarchy_roc <- roc(response = pint_hierarchy$TP, predict = abs(pint_hierarchy$strength))
        gc()
        pint_dedup_roc <- roc(response = pint_dedup$TP, predict = abs(pint_dedup$strength))
        gc()
        pint_unbiased_roc <- roc(response = pint_unbiased$TP, predict = abs(pint_unbiased$strength))
        gc()
        whinter_roc <- roc(response = whinter$TP, predict = abs(whinter$strength))
        if (use_glint) {
            glint_roc <- roc(response = glint$TP, predict = abs(glint$strength))
        }
    }
    gc()

    if (use_glint) {
      roc_list = list(Pint = pint_roc,
                    "Pint hierarchy" = pint_hierarchy_roc,
                    #Pint_Dedup = pint_dedup_roc,
                    #Pint_Unbiased = pint_unbiased_roc,
                    Whinter = whinter_roc,
                    Glinternet = glint_roc)
      glint_annotation =  sprintf("Glinternet auc: %.2f", round(glint_roc$auc, digits = 2))
    } else {
      roc_list = list(Pint = pint_roc,
                    "Pint hierarchy" = pint_hierarchy_roc,
                   Whinter = whinter_roc)
      glint_annotation = ""
    }
    roc_plot <-
        ggroc(roc_list) +
        theme_bw() +
        geom_segment(x = -1, y = 0, xend = 0, yend = 1, color = "#4c72b0") +
        annotate("text", x = 0.2, y = 0.4, label = sprintf("Pint auc: %.2f", round(pint_roc$auc, digits = 2))) +
        annotate("text", x = 0.2, y = 0.3, label = sprintf("Pint hierarchy auc: %.2f", round(pint_hierarchy_roc$auc, digits = 2))) +
        ## annotate("text", x = 0.2, y = 0.4, label = sprintf("Pint_Dedup auc: %.2f", round(pint_dedup_roc$auc, digits = 2))) +
        ## annotate("text", x = 0.2, y = 0.3, label = sprintf("Pint_Unbiased auc: %.2f", round(pint_unbiased_roc$auc, digits = 2))) +
        annotate("text", x = 0.2, y = 0.2, label = sprintf("Whinter auc: %.2f", round(whinter_roc$auc, digits = 2))) +
        annotate("text", x = 0.2, y = 0.1, label = glint_annotation)
    ggsave(roc_plot, file = output_filename, width = 7, height = 5)
}

plot_pint_3way_rocs <- function(pint, pint_dedup, pint_hierarchy, pint_unbiased, include_FN, output_filename, use_equiv_tp) {
    pint <- all_pint_smrys |> filter(type == "3way")
    pint_hierarchy <- all_pint_hierarchy_smrys |> filter(type == "3way")
    pint_dedup <- all_pint_dedup_smrys |> filter(type == "3way")
    pint_unbiased <- all_pint_unbiased_smrys |> filter(type == "3way")
    if (include_FN == FALSE) {
      pint <- pint |> filter(found == TRUE)
      pint_hierarchy <- pint_hierarchy |> filter(found == TRUE)
      pint_dedup <- pint_dedup |> filter(found == TRUE)
      pint_unbiased <- pint_unbiased |> filter(found == TRUE)
    }
    if (use_equiv_tp) {
        pint_roc <- roc(response = pint$equiv_tp & pint$TP, predict = abs(pint$strength))
        pint_hierarchy_roc <- roc(response = pint_hierarchy$equiv_tp & pint_hierarchy$TP, predict = abs(pint_hierarchy$strength))
        pint_dedup_roc <- roc(response = pint_dedup$equiv_tp & pint_dedup$TP, predict = abs(pint_dedup$strength))
        pint_unbiased_roc <- roc(response = pint_unbiased$equiv_tp & pint_unbiased$TP, predict = abs(pint_unbiased$strength))
    } else {
        pint_roc <- roc(response = pint$TP, predict = abs(pint$strength))
        pint_hierarchy_roc <- roc(response = pint_hierarchy$TP, predict = abs(pint_hierarchy$strength))
        pint_dedup_roc <- roc(response = pint_dedup$TP, predict = abs(pint_dedup$strength))
        pint_unbiased_roc <- roc(response = pint_unbiased$TP, predict = abs(pint_unbiased$strength))
    }
    gc()

    roc_list = list(Pint = pint_roc,
                    "Pint hierarchy" = pint_hierarchy_roc,
                    "Pint deduplicated" = pint_dedup_roc
                )
    roc_plot <-
        ggroc(roc_list) +
        theme_bw() +
        geom_segment(x = -1, y = 0, xend = 0, yend = 1, color = "#4c72b0") +
        annotate("text", x = 0.2, y = 0.4, label = sprintf("Pint auc: %.2f", round(pint_roc$auc, digits = 2))) +
        annotate("text", x = 0.2, y = 0.3, label = sprintf("Pint hierarchy auc: %.2f", round(pint_hierarchy_roc$auc, digits = 2))) +
        annotate("text", x = 0.2, y = 0.2, label = sprintf("Pint deduplicated auc: %.2f", round(pint_dedup_roc$auc, digits = 2)))
        ## annotate("text", x = 0.2, y = 0.3, label = sprintf("Pint_Unbiased auc: %.2f", round(pint_unbiased_roc$auc, digits = 2))) +
    ggsave(roc_plot, file = output_filename, width = 7, height = 5)
}

plot_pint_equiv_vs_tp <- function(pint, include_FN, output_filename) {
    pint <- all_pint_smrys |> filter(type == "3way")
    if (include_FN == FALSE) {
      pint <- pint |> filter(found == TRUE)
      pint_hierarchy <- pint_hierarchy |> filter(found == TRUE)
      pint_dedup <- pint_dedup |> filter(found == TRUE)
      pint_unbiased <- pint_unbiased |> filter(found == TRUE)
    }
    pint_tp_roc <- roc(response = pint$TP, predict = abs(pint$strength))
    pint_equiv_roc <- roc(response = (pint$equiv_tp | pint$TP), predict = abs(pint$strength))
    gc()

    roc_list = list(Pint = pint_tp_roc,
                    "Pint (equiv_tp)" = pint_equiv_roc)
    roc_plot <-
        ggroc(roc_list) +
        theme_bw() +
        geom_segment(x = -1, y = 0, xend = 0, yend = 1, color = "#4c72b0") +
        annotate("text", x = 0.2, y = 0.4, label = sprintf("Pint auc: %.2f", round(pint_tp_roc$auc, digits = 2))) +
        annotate("text", x = 0.2, y = 0.3, label = sprintf("Pint equiv_tp auc: %.2f", round(pint_equiv_roc$auc, digits = 2)))
    ggsave(roc_plot, file = output_filename, width = 7, height = 5)
}

#all_sets_pint_summaries <- c()
#all_sets_pint_hierarchy_summaries <- c()
#all_sets_pint_dedup_summaries <- c()
#all_sets_pint_unbiased_summaries <- c()
#all_sets_whinter_summaries <- c()
#all_sets_glint_summaries <- c()
#all_sets_times <- c()
for (bench_set in bench_sets) {
    if (bench_set == "wide_only_10k") {
        use_glint = FALSE
    } else {
        use_glint = TRUE
    }

    dir <- paste("whinter_pint_comparison", bench_set, sep = "/")

    rds_files <- list.files(dir, pattern = "all.rds", recursive = TRUE, full.names = TRUE)

    output <- c()
    for (file in rds_files) {
        print(sprintf("reading file %s", file))
        output <- cbind(output, readRDS(file))
        print("done")
    }

    output <- data.frame(output)

    all_pint_smrys <- foreach(set = output, .combine = rbind) %do% { set$pint_smry }
    all_pint_hierarchy_smrys <- foreach(set = output, .combine = rbind) %do% { set$pint_hierarchy_smry }
    all_pint_dedup_smrys <- foreach(set = output, .combine = rbind) %do% { set$pint_dedup_smry }
    all_pint_unbiased_smrys <- foreach(set = output, .combine = rbind) %do% { set$pint_unbiased_smry }
    all_whinter_smrys <- foreach(set = output, .combine = rbind) %do% { set$whinter_smry }
    if (use_glint) {
        all_glint_smrys <- foreach(set = output, .combine = rbind) %do% { set$glint_smry }
    } else {
      all_glint_smrys = NA
    }

    #all_sets_pint_summaries <- rbind(all_sets_pint_summaries, all_pint_smrys)
    #gc()
    #all_sets_pint_hierarchy_summaries <- rbind(all_sets_pint_hierarchy_summaries, all_pint_hierarchy_smrys)
    #all_sets_pint_dedup_summaries <- rbind(all_sets_pint_dedup_summaries, all_pint_dedup_smrys)
    #all_sets_pint_unbiased_summaries <- rbind(all_sets_pint_unbiased_summaries, all_pint_unbiased_smrys)
    #all_sets_whinter_summaries <- rbind(all_sets_whinter_summaries, all_whinter_smrys)
    #if (use_glint) {
    #    all_sets_glint_summaries <- rbind(all_sets_glint_summaries, all_glint_smrys)
    #}

    if (use_glint) {
        all_times <- foreach(out = output, .combine = rbind) %do% {
            data.frame(pint = out$pint_time[3],
                        pint_hierarchy = out$pint_hierarchy_time[3],
                        pint_dedup = out$pint_dedup_time[3],
                        pint_unbiased = out$pint_unbiased_time[3],
                        whinter = out$whinter_time[3],
                        glint = out$glint_time[3])
            }
    } else {
        all_times <- foreach(out = output, .combine = rbind) %do% {
            data.frame(pint = out$pint_time[3],
                        pint_hierarchy = out$pint_hierarchy_time[3],
                        pint_dedup = out$pint_dedup_time[3],
                        pint_unbiased = out$pint_unbiased_time[3],
                        whinter = out$whinter_time[3])
        }
    }
    #all_sets_times <- rbind(all_sets_times, all_times)
    rm(output)
    gc()

    write.csv(summary(all_times), file=sprintf("times_summary_%s.csv", bench_set))
    if (use_glint) {
        mean_effects_found = data.frame(
            whinter = sum(all_whinter_smrys$found)/length(rds_files),
            pint = sum(all_pint_smrys$found)/length(rds_files),
            pint_hierarchy = sum(all_pint_hierarchy_smrys$found)/length(rds_files),
            glinternet = sum(all_glint_smrys$found)/length(rds_files)
        )
    } else {
        mean_effects_found = data.frame(
            whinter = sum(all_whinter_smrys$found)/length(rds_files),
            pint = sum(all_pint_smrys$found)/length(rds_files),
            pint_hierarchy = sum(all_pint_hierarchy_smrys$found)/length(rds_files)
        )
    }
    write.csv(mean_effects_found, file=sprintf("mean_effects_found_%s.csv", bench_set))

    plot_dir <- "plots/rocs/"
    if (! dir.exists(plot_dir)) {
      dir.create(plot_dir, recursive=TRUE)
    }
    plot_rocs(all_pint_smrys, all_pint_hierarchy_smrys, all_pint_dedup_smrys, all_pint_unbiased_smrys, all_glint_smrys, all_whinter_smrys, use_glint, "all", TRUE, sprintf("plots/rocs/all_comparison_roc_%s.pdf", bench_set))
    plot_rocs(all_pint_smrys, all_pint_hierarchy_smrys, all_pint_dedup_smrys, all_pint_unbiased_smrys, all_glint_smrys, all_whinter_smrys, use_glint, "main", TRUE, sprintf("plots/rocs/main_comparison_roc_%s.pdf", bench_set))
    plot_rocs(all_pint_smrys, all_pint_hierarchy_smrys, all_pint_dedup_smrys, all_pint_unbiased_smrys, all_glint_smrys, all_whinter_smrys, use_glint, "interaction", TRUE, sprintf("plots/rocs/int_comparison_roc_%s.pdf", bench_set))

    plot_rocs(all_pint_smrys, all_pint_hierarchy_smrys, all_pint_dedup_smrys, all_pint_unbiased_smrys, all_glint_smrys, all_whinter_smrys, use_glint, "all", FALSE, sprintf("plots/rocs/all_comparison_roc_%s_FoundOnly.pdf", bench_set))
    plot_rocs(all_pint_smrys, all_pint_hierarchy_smrys, all_pint_dedup_smrys, all_pint_unbiased_smrys, all_glint_smrys, all_whinter_smrys, use_glint, "main", FALSE, sprintf("plots/rocs/main_comparison_roc_%s_FoundOnly.pdf", bench_set))
    plot_rocs(all_pint_smrys, all_pint_hierarchy_smrys, all_pint_dedup_smrys, all_pint_unbiased_smrys, all_glint_smrys, all_whinter_smrys, use_glint, "interaction", FALSE, sprintf("plots/rocs/int_comparison_roc_%s_FoundOnly.pdf", bench_set))

    # equiv_tp
    ## all
    plot_rocs(all_pint_smrys, all_pint_hierarchy_smrys, all_pint_dedup_smrys, all_pint_unbiased_smrys, all_glint_smrys, all_whinter_smrys, use_glint,
              "all", TRUE, sprintf("plots/rocs/all_comparison_roc_equiv_%s.pdf", bench_set), use_equiv_tp=TRUE)
    ## found_only
    plot_rocs(all_pint_smrys, all_pint_hierarchy_smrys, all_pint_dedup_smrys, all_pint_unbiased_smrys, all_glint_smrys, all_whinter_smrys, use_glint,
              "all", FALSE, sprintf("plots/rocs/all_comparison_roc_equiv_%s_FoundOnly.pdf", bench_set), use_equiv_tp=TRUE)
    # threeway
    plot_pint_3way_rocs(all_pint_smrys, all_pint_hierarchy_smrys, all_pint_dedup_smrys, all_pint_unbiased_smrys, TRUE,
                        sprintf("plots/rocs/3way_comparison_roc_%s.pdf", bench_set), use_equiv_tp=FALSE)
    plot_pint_3way_rocs(all_pint_smrys, all_pint_hierarchy_smrys, all_pint_dedup_smrys, all_pint_unbiased_smrys, TRUE,
                        sprintf("plots/rocs/3way_comparison_roc_equiv_%s.pdf", bench_set), use_equiv_tp=TRUE)
    plot_pint_3way_rocs(all_pint_smrys, all_pint_hierarchy_smrys, all_pint_dedup_smrys, all_pint_unbiased_smrys, FALSE,
                        sprintf("plots/rocs/3way_comparison_roc_equiv_%s_FoundOnly.pdf", bench_set), use_equiv_tp=TRUE)
    plot_pint_equiv_vs_tp(all_pint_smrys, FALSE, sprintf("plots/rocs/3way_comparison_roc_equiv_vs_tp_%s_FoundOnly.pdf", bench_set))

    # set times
    if (use_glint) {
        melted_times <- melt(all_times |> select(pint, pint_hierarchy, whinter, glint), value.name = "Time", variable.name = "Method")
    } else {
        melted_times <- melt(all_times |> select(pint, pint_hierarchy, whinter), value.name = "Time", variable.name = "Method")
    }
    time_plot <- ggplot(melted_times, aes(x = Method, y = Time)) +
        geom_boxplot() +
        theme_bw() +
        scale_y_continuous(trans = "log2") +
        ylab("Time (s)") +
        expand_limits(y = 0)
    ggsave(time_plot, file=sprintf("plots/bench_times_%s.pdf", bench_set), width = 4, height = 4)

}

#use_glint=TRUE
#plot_rocs(all_sets_pint_summaries, all_sets_pint_hierarchy_smrys, all_sets_pint_dedup_smrys, all_sets_pint_unbiased_smrys, all_sets_glint_summaries, all_sets_whinter_summaries, use_glint, "all", TRUE, sprintf("plots/rocs/all_sets_comparison_roc.pdf"))
#plot_rocs(all_sets_pint_summaries, all_sets_pint_hierarchy_smrys, all_sets_pint_dedup_smrys, all_sets_pint_unbiased_smrys, all_sets_glint_summaries, all_sets_whinter_summaries, use_glint, "main", TRUE, sprintf("plots/rocs/main_comparison_roc.pdf"))
#plot_rocs(all_sets_pint_summaries, all_sets_pint_hierarchy_smrys, all_sets_pint_dedup_smrys, all_sets_pint_unbiased_smrys, all_sets_glint_summaries, all_sets_whinter_summaries, use_glint, "interaction", TRUE, sprintf("plots/rocs/int_comparison_roc.pdf"))
#
#plot_rocs(all_sets_pint_summaries, all_sets_pint_hierarchy_smrys, all_sets_pint_dedup_smrys, all_sets_pint_unbiased_smrys, all_sets_glint_summaries, all_sets_whinter_summaries, use_glint, "all", FALSE, sprintf("plots/rocs/all_sets_comparison_roc_FoundOnly.pdf"))
#plot_rocs(all_sets_pint_summaries, all_sets_pint_hierarchy_smrys, all_sets_pint_dedup_smrys, all_sets_pint_unbiased_smrys, all_sets_glint_summaries, all_sets_whinter_summaries, use_glint, "main", FALSE, sprintf("plots/rocs/main_comparison_roc_FoundOnly.pdf"))
#plot_rocs(all_sets_pint_summaries, all_sets_pint_hierarchy_smrys, all_sets_pint_dedup_smrys, all_sets_pint_unbiased_smrys, all_sets_glint_summaries, all_sets_whinter_summaries, use_glint, "interaction", FALSE, sprintf("plots/rocs/int_comparison_roc_FoundOnly.pdf"))
#
## plot times
#melted_times <- melt(all_sets_times, value.name = "Time", variable.name = "Method")
#time_plot <- ggplot(melted_times, aes(x = Method, y = Time)) +
#    geom_boxplot() +
#    theme_bw() +
#    scale_y_continuous(trans = "log2") +
#    ylab("Time (s)") +
#    expand_limits(y = 0)
#ggsave(time_plot, file="plots/all_bench_times.pdf", width = 6, height = 4)
#
#write.csv(summary(all_sets_times), file=sprintf("times_summary_all_sets.csv", bench_set))
