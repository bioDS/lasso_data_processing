#!/usr/bin/env Rscript
library(Pint)
library(pROC)
library(dplyr)
library(foreach)
library(doMC)
library(ggplot2)
library(reshape2)
library(future)

source("summary_functions.R")

registerDoMC(cores = 5)

#bench_sets <- c("simulated_small_data_sample", "3way", "wide_only_10k", "8k_only")
# bench_sets <- c("simulated_small_data_sample", "8k_only", "wide_only_10k")
# bench_sets <- c("wide_only_10k")
bench_sets <- c("wide_only_10k")

options(future.globals.maxSize = 16 * 1024^3)
get_roc_set_future <- function(response, predict) {
    f <- future({
        roc(response = response, predict = predict)
    }) %plan% sequential
    return(f)
}

plot_rocs <- function(all_pint_smrys, all_pint_hierarchy_smrys, all_pint_dedup_smrys,
                      all_pint_unbiased_smrys, all_glint_smrys, all_whinter_smrys,
                      use_glint, use_type, include_FN, output_filename, use_equiv_tp = FALSE) {
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
        pint_roc_future <- get_roc_set_future(response = pint$equiv_tp & pint$TP, predict = abs(pint$strength))
        gc()
        pint_hierarchy_roc_future <- get_roc_set_future(response = pint_hierarchy$equiv_tp & pint_hierarchy$TP, predict = abs(pint_hierarchy$strength))
        gc()
        #pint_dedup_roc_future <- get_roc_set_future(response = pint_dedup$equiv_tp & pint_dedup$TP, predict = abs(pint_dedup$strength))
        #gc()
        #pint_unbiased_roc_future <- get_roc_set_future(response = pint_unbiased$equiv_tp & pint_unbiased$TP, predict = abs(pint_unbiased$strength))
        #gc()
        whinter_roc_future <- get_roc_set_future(response = whinter$equiv_tp & whinter$TP, predict = abs(whinter$strength))
        if (use_glint) {
            glint_roc_future <- get_roc_set_future(response = glint$equiv_tp & glint$TP, predict = abs(glint$strength))
        }
    } else {
        #pint_roc <- roc(response = pint$TP, predict = abs(pint$strength))
        pint_roc_future <- get_roc_set_future(response = pint$TP, predict = abs(pint$strength))
        gc()
        pint_hierarchy_roc_future <- get_roc_set_future(response = pint_hierarchy$TP, predict = abs(pint_hierarchy$strength))
        gc()
        #pint_dedup_roc <- roc(response = pint_dedup$TP, predict = abs(pint_dedup$strength))
        #gc()
        #pint_unbiased_roc <- roc(response = pint_unbiased$TP, predict = abs(pint_unbiased$strength))
        #gc()
        whinter_roc_future <- get_roc_set_future(response = whinter$TP, predict = abs(whinter$strength))
        if (use_glint) {
            glint_roc_future <- get_roc_set_future(response = glint$TP, predict = abs(glint$strength))
        }
    }
    gc()

    if (use_glint) {
        roc_list <- list(
            Pint = value(pint_roc_future),
            "Pint hierarchy" = value(pint_hierarchy_roc_future),
            # Pint_Dedup = value(pint_dedup_roc_future),
            # Pint_Unbiased = value(pint_unbiased_roc_future),
            Whinter = value(whinter_roc_future),
            Glinternet = value(glint_roc_future)
        )
        glint_annotation <- sprintf("Glinternet auc: %.2f", round(value(glint_roc_future)$auc, digits = 2))
    } else {
        roc_list <- list(
            Pint = value(pint_roc_future),
            "Pint hierarchy" = value(pint_hierarchy_roc_future),
            Whinter = value(whinter_roc_future)
        )
        glint_annotation <- ""
    }
    roc_plot <-
        ggroc(roc_list) +
        theme_bw() +
        geom_segment(x = -1, y = 0, xend = 0, yend = 1, color = "#4c72b0") +
        annotate("text", x = 0.2, y = 0.4, label = sprintf("Pint auc: %.2f", round(value(pint_roc_future)$auc, digits = 2))) +
        annotate("text", x = 0.2, y = 0.3, label = sprintf("Pint hierarchy auc: %.2f", round(value(pint_hierarchy_roc_future)$auc, digits = 2))) +
        ## annotate("text", x = 0.2, y = 0.4, label = sprintf("Pint_Dedup auc: %.2f", round(pint_dedup_roc$auc, digits = 2))) +
        ## annotate("text", x = 0.2, y = 0.3, label = sprintf("Pint_Unbiased auc: %.2f", round(pint_unbiased_roc$auc, digits = 2))) +
        annotate("text", x = 0.2, y = 0.2, label = sprintf("Whinter auc: %.2f", round(value(whinter_roc_future)$auc, digits = 2))) +
        annotate("text", x = 0.2, y = 0.1, label = glint_annotation)
    ggsave(roc_plot, file = output_filename, width = 6, height = 4)
    roc_df = data.frame(
        Pint = value(pint_roc_future)$auc,
        Pint_Hierarchy = value(pint_hierarchy_roc_future)$auc,
        Whinter = value(whinter_roc_future)$auc,
        Glinternet = value(glint_roc_future)$auc
    )
    write.csv(roc_df, sprintf("%s-values.csv", output_filename))
}

plot_rocs_anyfound <- function(all_pint_smrys, all_pint_hierarchy_smrys, all_pint_dedup_smrys,
                      all_pint_unbiased_smrys, all_pint_pair_smrys, all_glint_smrys, all_whinter_smrys,
                      use_glint, use_pint_pair, use_type, output_filename, use_equiv_tp = FALSE) {
    ## main only
    gc()

    if (use_type != "all") {
        pint <- all_pint_smrys |> filter(type == use_type)
        pint_hierarchy <- all_pint_hierarchy_smrys |> filter(type == use_type)
        # pint_dedup <- all_pint_dedup_smrys |> filter(type == use_type)
        #pint_unbiased <- all_pint_unbiased_smrys |> filter(type == use_type)
        if (use_glint) {
            glint <- all_glint_smrys |> filter(type == use_type)
        }
        if (use_pint_pair) {
            pint_pair <- all_pint_pair_smrys |> filter(type == use_type)
        }
        whinter <- all_whinter_smrys |> filter(type == use_type)
    } else {
        pint <- all_pint_smrys
        pint_hierarchy <- all_pint_hierarchy_smrys
        # pint_dedup <- all_pint_dedup_smrys
        #pint_unbiased <- all_pint_unbiased_smrys
        if (use_glint) {
            glint <- all_glint_smrys
        }
        if (use_pint_pair) {
            pint_pair <- all_pint_pair_smrys
        }
        whinter <- all_whinter_smrys
    }

    if ("gene_k" %in% names(pint)) {
        # use 3way
        use_genes <- c("gene_i", "gene_j", "gene_k")
    } else {
        use_genes <- c("gene_i", "gene_j")
    }

    if (use_pint_pair) {
        anyfound_ids <- rbind(
            pint            |> filter(found==TRUE) |> select(all_of(use_genes)),
            pint_hierarchy  |> filter(found==TRUE) |> select(all_of(use_genes)),
            # pint_dedup      |> filter(found==TRUE) |> select(all_of(use_genes)),
            pint_pair       |> filter(found==TRUE) |> select(all_of(use_genes)),
            whinter         |> filter(found==TRUE) |> select(all_of(use_genes)),
            glint           |> filter(found==TRUE) |> select(all_of(use_genes))
        ) |> unique()
    } else {
        anyfound_ids <- rbind(
            pint            |> filter(found==TRUE) |> select(all_of(use_genes)),
            pint_hierarchy  |> filter(found==TRUE) |> select(all_of(use_genes)),
            # pint_dedup      |> filter(found==TRUE) |> select(all_of(use_genes)),
            whinter         |> filter(found==TRUE) |> select(all_of(use_genes)),
            glint           |> filter(found==TRUE) |> select(all_of(use_genes))
        ) |> unique()
    }

    pint_anyfound <- left_join(anyfound_ids, pint, by=all_of(use_genes))
    # pint_dedup_anyfound <- left_join(anyfound_ids, pint_dedup, by=all_of(use_genes))
    pint_hierarchy_anyfound <- left_join(anyfound_ids, pint_hierarchy, by=all_of(use_genes))
    if (use_pint_pair) {
        pint_pair_anyfound <- left_join(anyfound_ids, pint_pair, by=all_of(use_genes))
    }
    whinter_anyfound <- left_join(anyfound_ids, whinter, by=all_of(use_genes))
    glint_anyfound <- left_join(anyfound_ids, glint, by=all_of(use_genes))

    #all_ids <- full_join(
    #    pint |> select("gene_i", "gene_j", "gene_k"),
    #    pint_hierarchy |> select("gene_i", "gene_j", "gene_k"),
    #    by=c("gene_i", "gene_j", "gene_k")
    #)

    if (use_equiv_tp) {
        pint_anyfound$TP <- pint_anyfound$equiv_tp & pint_anyfound$TP
        pint_hierarchy_anyfound$TP <- pint_hierarchy_anyfound$equiv_tp & pint_hierarchy_anyfound$TP
        # pint_dedup_anyfound$TP <- pint_dedup_anyfound$equiv_tp & pint_dedup_anyfound$TP
        whinter_anyfound$TP <- whinter_anyfound$equiv_tp & whinter_anyfound$TP
        glint_anyfound$TP <- glint_anyfound$equiv_tp & glint_anyfound$TP
        if (use_pint_pair) {
            pint_pair_anyfound$TP <- pint_pair_anyfound$equiv_tp & pint_pair_anyfound$TP
        }
    }

    pint_roc_future <- get_roc_set_future(response = pint_anyfound$TP, predict = abs(pint_anyfound$strength))
    gc()
    pint_hierarchy_roc_future <- get_roc_set_future(response = pint_hierarchy_anyfound$TP, predict = abs(pint_hierarchy_anyfound$strength))
    gc()
    # pint_dedup_roc_future <- get_roc_set_future(response = pint_dedup_anyfound$TP, predict = abs(pint_dedup_anyfound$strength))
    gc()
    whinter_roc_future <- get_roc_set_future(response = whinter_anyfound$TP, predict = abs(whinter_anyfound$strength))
    if (use_glint) {
        glint_roc_future <- get_roc_set_future(response = glint_anyfound$TP, predict = abs(glint_anyfound$strength))
    }
    if (use_pint_pair) {
        pint_pair_roc_future <- get_roc_set_future(response = pint_pair_anyfound$TP, predict = abs(pint_pair_anyfound$strength))
    }
    gc()

    roc_list <- list(
        Pint = value(pint_roc_future),
        # "Pint deduplicated" = value(pint_dedup_roc_future),
        "Pint hierarchy" = value(pint_hierarchy_roc_future),
        Whinter = value(whinter_roc_future)
    )
    glint_annotation <- ""
    pint_pair_annotation <- ""
    if (use_glint) {
        roc_list <- append(
            roc_list,
            list(Glinternet = value(glint_roc_future))
        )
        glint_annotation <- sprintf("Glinternet auc: %.2f", round(value(glint_roc_future)$auc, digits = 2))
    }
    if (use_pint_pair) {
        roc_list <- append(
            roc_list,
            list("Pint (pairwise only)" = value(pint_pair_roc_future))
        )
        pint_pair_annotation <- sprintf("Pint (pairwise) auc: %.2f", round(value(pint_pair_roc_future)$auc, digits = 2))
    }
    roc_plot <-
        ggroc(roc_list) +
        theme_bw() +
        geom_segment(x = -1, y = 0, xend = 0, yend = 1, color = "#4c72b0") +
        annotate("text", x = 0.2, y = 0.5, label = sprintf("Pint auc: %.2f", round(value(pint_roc_future)$auc, digits = 2))) +
        # annotate("text", x = 0.2, y = 0.4, label = sprintf("Pint (dedup) auc: %.2f", round(pint_dedup_roc$auc, digits = 2))) +
        annotate("text", x = 0.2, y = 0.4, label = sprintf("Pint hierarchy auc: %.2f", round(value(pint_hierarchy_roc_future)$auc, digits = 2))) +
        ## annotate("text", x = 0.2, y = 0.3, label = sprintf("Pint_Unbiased auc: %.2f", round(pint_unbiased_roc$auc, digits = 2))) +
        annotate("text", x = 0.2, y = 0.3, label = sprintf("Whinter auc: %.2f", round(value(whinter_roc_future)$auc, digits = 2))) +
        annotate("text", x = 0.2, y = 0.2, label = glint_annotation) +
        annotate("text", x = 0.2, y = 0.1, label = pint_pair_annotation)
    ggsave(roc_plot, file = output_filename, width = 6, height = 4)
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

    roc_list <- list(
        Pint = pint_roc,
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
    ggsave(roc_plot, file = output_filename, width = 6, height = 4)
    roc_df <- data.frame(
        Pint = pint_roc$auc,
        pint_hierarchy = pint_hierarchy_roc$auc,
        pint_dedup = pint_dedup_roc$auc
    )
    write.csv(roc_df, sprintf("%s.csv", output_filename))
}

plot_pint_equiv_vs_tp <- function(pint, include_FN, output_filename) {
    pint <- all_pint_smrys |> filter(type == "3way")
    if (include_FN == FALSE) {
        pint <- pint |> filter(found == TRUE)
    }
    pint_tp_roc <- roc(response = pint$TP, predict = abs(pint$strength))
    pint_equiv_roc <- roc(response = (pint$equiv_tp | pint$TP), predict = abs(pint$strength))
    gc()

    roc_list <- list(
        Pint = pint_tp_roc,
        "Pint (equiv_tp)" = pint_equiv_roc
    )
    roc_plot <-
        ggroc(roc_list) +
        theme_bw() +
        geom_segment(x = -1, y = 0, xend = 0, yend = 1, color = "#4c72b0") +
        annotate("text", x = 0.2, y = 0.4, label = sprintf("Pint auc: %.2f", round(pint_tp_roc$auc, digits = 2))) +
        annotate("text", x = 0.2, y = 0.3, label = sprintf("Pint equiv_tp auc: %.2f", round(pint_equiv_roc$auc, digits = 2)))
    ggsave(roc_plot, file = output_filename, width = 6, height = 4)
}

fix_na_strengths <- function(smry) {
    if (sum(is.na(smry$strength)) > 0) {
        smry[is.na(smry$strength), ]$strength <- 0
    }
    return(smry)
}

convert_time <- function(time) {
    if (length(time) > 1) {
        return(time[[3]])
    } else {
        return(time[[1]])
    }
}

for (bench_set in bench_sets) {
    use_glint <- TRUE
    if (bench_set == "3way") {
        use_pint_pair <- TRUE
    } else {
        use_pint_pair <- FALSE
        all_pint_pair_smrys <- NA
    }

    dir <- paste("whinter_pint_comparison", bench_set, sep = "/")

    rds_files <- list.files(dir, pattern = "all.rds", recursive = TRUE, full.names = TRUE)

    output <- c()
    output <- foreach(file = rds_files, .combine=cbind) %dopar% {
        print(sprintf("reading file %s", file))
        new_out <- readRDS(file)
        new_out$file <- file
        print("done")
        return(new_out)
    }
    output <- data.frame(output)

    if (bench_set == "3way") {
        use_values <- c("gene_i", "gene_j", "gene_k", "TP", "strength", "found", "type")
    } else {
        use_values <- c("gene_i", "gene_j", "TP", "strength", "found", "type")
    }
    all_pint_smrys <- foreach(set = output, .combine = rbind) %do% {
        set$pint_smry |> select(all_of(use_values)) |> mutate(file=set$file)
    }
    gc()
    all_pint_hierarchy_smrys <- foreach(set = output, .combine = rbind) %do% {
        set$pint_hierarchy_smry |> select(all_of(use_values)) |> mutate(file=set$file)
    }
    gc()
    all_pint_dedup_smrys <- foreach(set = output, .combine = rbind) %do% {
        set$pint_dedup_smry |> select(all_of(use_values)) |> mutate(file=set$file)
    }
    gc()
    all_pint_unbiased_smrys <- foreach(set = output, .combine = rbind) %do% {
        set$pint_unbiased_smry |> select(all_of(use_values)) |> mutate(file=set$file)
    }
    gc()
    if (use_pint_pair) {
        all_pint_pair_smrys <- foreach(set = output, .combine = rbind) %do% {
            set$pint_pair_smry |> select(all_of(use_values)) |> mutate(file=set$file)
        }
        gc()
    }
    all_whinter_smrys <- foreach(set = output, .combine = rbind) %do% {
        set$whinter_smry |> select(all_of(use_values)) |> mutate(file=set$file)
    }
    gc()
    if (use_glint) {
        all_glint_smrys <- foreach(set = output, .combine = rbind) %do% {
            set$glint_smry |> select(all_of(use_values)) |> mutate(file=set$file)
        }
        gc()
    } else {
        all_glint_smrys <- NA
    }
    #f <- future({
    #    roc(response = response, predict = predict)
    #}) %plan% multicore
    #all_pint_smrys <- future({fix_na_strengths(all_pint_smrys)}) %plan% multicore
    #all_pint_dedup_smrys <- future({fix_na_strengths(all_pint_dedup_smrys)}) %plan% multicore
    #all_pint_hierarchy_smrys <- future({fix_na_strengths(all_pint_hierarchy_smrys)}) %plan% multicore
    #all_pint_unbiased_smrys <- future({fix_na_strengths(all_pint_unbiased_smrys)}) %plan% multicore
    #all_whinter_smrys <- future({fix_na_strengths(all_whinter_smrys)}) %plan% multicore
    #all_glint_smrys <- future({fix_na_strengths(all_glint_smrys)}) %plan% multicore
    #if (use_pint_pair) {
    #    all_pint_pair_smrys <- future({fix_na_strengths(all_pint_pair_smrys)}) %plan% multicore
    #}
    #all_pint_smrys <- value(all_pint_smrys)
    #all_pint_dedup_smrys <- value(all_pint_dedup_smrys)
    #all_pint_hierarchy_smrys <- value(all_pint_hierarchy_smrys)
    #all_pint_unbiased_smrys <- value(all_pint_unbiased_smrys)
    #all_whinter_smrys <- value(all_whinter_smrys)
    #all_glint_smrys <- value(all_glint_smrys)
    #if (use_pint_pair) {
    #    all_pint_pair_smrys <- value(all_pint_pair_smrys)
    #}
    all_pint_smrys <- fix_na_strengths(all_pint_smrys)
    all_pint_dedup_smrys <- fix_na_strengths(all_pint_dedup_smrys)
    all_pint_hierarchy_smrys <- fix_na_strengths(all_pint_hierarchy_smrys)
    all_pint_unbiased_smrys <- fix_na_strengths(all_pint_unbiased_smrys)
    all_whinter_smrys <- fix_na_strengths(all_whinter_smrys)
    all_glint_smrys <- fix_na_strengths(all_glint_smrys)
    if (use_pint_pair) {
        all_pint_pair_smrys <- fix_na_strengths(all_pint_pair_smrys)
    }
    gc()

    if (use_glint) {
        all_times <- foreach(out = output, .combine = rbind) %do% {
            data.frame(
                Pint = convert_time(out$pint_time),
                "Pint hierarchy" = convert_time(out$pint_hierarchy_time),
                "Pint deduplicated" = convert_time(out$pint_dedup_time),
                "Pint unbiased" = convert_time(out$pint_unbiased_time),
                WHInter = convert_time(out$whinter_time),
                Glinternet = convert_time(out$glint_time)
            )
        }
    } else {
        all_times <- foreach(out = output, .combine = rbind) %do% {
            data.frame(
                Pint = convert_time(out$pint_time),
                "Pint hierarchy" = convert_time(out$pint_hierarchy_time),
                "Pint deduplicated" = convert_time(out$pint_dedup_time),
                "Pint unbiased" = convert_time(out$pint_unbiased_time),
                WHInter = convert_time(out$whinter_time),
            )
        }
    }

    rm(output)
    gc()

    write.csv(summary(all_times), file = sprintf("times_summary_%s.csv", bench_set))
    if (use_glint) {
        mean_effects_found <- data.frame(
            whinter = sum(all_whinter_smrys$found) / length(rds_files),
            pint = sum(all_pint_smrys$found) / length(rds_files),
            pint_hierarchy = sum(all_pint_hierarchy_smrys$found) / length(rds_files),
            glinternet = sum(all_glint_smrys$found) / length(rds_files)
        )
    } else {
        mean_effects_found <- data.frame(
            whinter = sum(all_whinter_smrys$found) / length(rds_files),
            pint = sum(all_pint_smrys$found) / length(rds_files),
            pint_hierarchy = sum(all_pint_hierarchy_smrys$found) / length(rds_files)
        )
    }
    write.csv(mean_effects_found, file = sprintf("mean_effects_found_%s.csv", bench_set))

    plot_dir <- "plots/rocs/"
    if (!dir.exists(plot_dir)) {
        dir.create(plot_dir, recursive = TRUE)
    }
    plot_rocs(
        all_pint_smrys, all_pint_hierarchy_smrys, all_pint_dedup_smrys,
        all_pint_unbiased_smrys, all_glint_smrys, all_whinter_smrys, use_glint, "all",
        TRUE, sprintf("plots/rocs/all_comparison_roc_%s.pdf", bench_set)
    )
    plot_rocs(
        all_pint_smrys, all_pint_hierarchy_smrys, all_pint_dedup_smrys,
        all_pint_unbiased_smrys, all_glint_smrys, all_whinter_smrys, use_glint, "main",
        TRUE, sprintf("plots/rocs/main_comparison_roc_%s.pdf", bench_set)
    )
    plot_rocs(
        all_pint_smrys, all_pint_hierarchy_smrys, all_pint_dedup_smrys,
        all_pint_unbiased_smrys, all_glint_smrys, all_whinter_smrys, use_glint, "interaction",
        TRUE, sprintf("plots/rocs/int_comparison_roc_%s.pdf", bench_set)
    )

    plot_rocs(
        all_pint_smrys, all_pint_hierarchy_smrys, all_pint_dedup_smrys,
        all_pint_unbiased_smrys, all_glint_smrys, all_whinter_smrys, use_glint, "all",
        FALSE, sprintf("plots/rocs/all_comparison_roc_%s_FoundOnly.pdf", bench_set)
    )
    plot_rocs(
        all_pint_smrys, all_pint_hierarchy_smrys, all_pint_dedup_smrys,
        all_pint_unbiased_smrys, all_glint_smrys, all_whinter_smrys, use_glint, "main",
        FALSE, sprintf("plots/rocs/main_comparison_roc_%s_FoundOnly.pdf", bench_set)
    )
    plot_rocs(
        all_pint_smrys, all_pint_hierarchy_smrys, all_pint_dedup_smrys,
        all_pint_unbiased_smrys, all_glint_smrys, all_whinter_smrys, use_glint, "interaction",
        FALSE, sprintf("plots/rocs/int_comparison_roc_%s_FoundOnly.pdf", bench_set)
    )

    # plot rocs including only the effects found by at least one method
    #plot_rocs_anyfound(
    #    all_pint_smrys, all_pint_hierarchy_smrys, all_pint_dedup_smrys,
    #    all_pint_unbiased_smrys, all_pint_pair_smrys, all_glint_smrys, all_whinter_smrys, use_glint, use_pint_pair, "all",
    #    sprintf("plots/rocs/all_comparison_roc_%s_anyfound.pdf", bench_set)
    #)
    #plot_rocs_anyfound(
    #    all_pint_smrys, all_pint_hierarchy_smrys, all_pint_dedup_smrys,
    #    all_pint_unbiased_smrys, all_pint_pair_smrys, all_glint_smrys, all_whinter_smrys, use_glint, use_pint_pair, "main",
    #    sprintf("plots/rocs/main_comparison_roc_%s_anyfound.pdf", bench_set)
    #)
    #plot_rocs_anyfound(
    #    all_pint_smrys, all_pint_hierarchy_smrys, all_pint_dedup_smrys,
    #    all_pint_unbiased_smrys, all_pint_pair_smrys, all_glint_smrys, all_whinter_smrys, use_glint, use_pint_pair, "interaction",
    #    sprintf("plots/rocs/int_comparison_roc_%s_anyfound.pdf", bench_set)
    #)
    ### anyfound using equiv_tp
    #plot_rocs_anyfound(
    #    all_pint_smrys, all_pint_hierarchy_smrys, all_pint_dedup_smrys,
    #    all_pint_unbiased_smrys, all_pint_pair_smrys, all_glint_smrys, all_whinter_smrys, use_glint, use_pint_pair, "all",
    #    sprintf("plots/rocs/all_comparison_roc_%s_equiv_tp_anyfound.pdf", bench_set), use_equiv_tp = TRUE
    #)
    #plot_rocs_anyfound(
    #    all_pint_smrys, all_pint_hierarchy_smrys, all_pint_dedup_smrys,
    #    all_pint_unbiased_smrys, all_pint_pair_smrys, all_glint_smrys, all_whinter_smrys, use_glint, use_pint_pair, "main",
    #    sprintf("plots/rocs/main_comparison_roc_%s_equiv_tp_anyfound.pdf", bench_set), use_equiv_tp = TRUE
    #)
    #plot_rocs_anyfound(
    #    all_pint_smrys, all_pint_hierarchy_smrys, all_pint_dedup_smrys,
    #    all_pint_unbiased_smrys, all_pint_pair_smrys, all_glint_smrys, all_whinter_smrys, use_glint, use_pint_pair, "interaction",
    #    sprintf("plots/rocs/int_comparison_roc_%s_equiv_tp_anyfound.pdf", bench_set), use_equiv_tp = TRUE
    #)

    # equiv_tp
    ## all
    #plot_rocs(all_pint_smrys, all_pint_hierarchy_smrys, all_pint_dedup_smrys, all_pint_unbiased_smrys, all_glint_smrys, all_whinter_smrys, use_glint,
    #    "all", TRUE, sprintf("plots/rocs/all_comparison_roc_equiv_%s.pdf", bench_set),
    #    use_equiv_tp = TRUE
    #)
    ### found_only
    #plot_rocs(all_pint_smrys, all_pint_hierarchy_smrys, all_pint_dedup_smrys, all_pint_unbiased_smrys, all_glint_smrys, all_whinter_smrys, use_glint,
    #    "all", FALSE, sprintf("plots/rocs/all_comparison_roc_equiv_%s_FoundOnly.pdf", bench_set),
    #    use_equiv_tp = TRUE
    #)
    # threeway
    if (bench_set == "3way") {
        plot_pint_3way_rocs(all_pint_smrys, all_pint_hierarchy_smrys, all_pint_dedup_smrys, all_pint_unbiased_smrys, TRUE,
            sprintf("plots/rocs/3way_comparison_roc_%s.pdf", bench_set),
            use_equiv_tp = FALSE
        )
        #plot_pint_3way_rocs(all_pint_smrys, all_pint_hierarchy_smrys, all_pint_dedup_smrys, all_pint_unbiased_smrys, TRUE,
        #    sprintf("plots/rocs/3way_comparison_roc_equiv_%s.pdf", bench_set),
        #    use_equiv_tp = TRUE
        #)
        #plot_pint_3way_rocs(all_pint_smrys, all_pint_hierarchy_smrys, all_pint_dedup_smrys, all_pint_unbiased_smrys, FALSE,
        #    sprintf("plots/rocs/3way_comparison_roc_equiv_%s_FoundOnly.pdf", bench_set),
        #    use_equiv_tp = TRUE
        #)
        #plot_pint_equiv_vs_tp(all_pint_smrys, FALSE, sprintf("plots/rocs/3way_comparison_roc_equiv_vs_tp_%s_FoundOnly.pdf", bench_set))
    }

    # set times
    if (use_glint) {
        melted_times <- melt(all_times |> select(Pint, "Pint.hierarchy", WHInter, Glinternet), value.name = "Time", variable.name = "Method")
    } else {
        melted_times <- melt(all_times |> select(Pint, "Pint.hierarchy", WHInter), value.name = "Time", variable.name = "Method")
    }
    time_plot <- ggplot(melted_times, aes(x = Method, y = Time)) +
        geom_boxplot() +
        theme_bw() +
        scale_y_continuous(trans = "log2") +
        ylab("Time (s)") +
        scale_x_discrete(labels=c('Pint', 'Pint (hier)', 'WHInter', 'Glinternet')) +
        expand_limits(y = 0)
    ggsave(time_plot, file = sprintf("plots/bench_times_%s.pdf", bench_set), width = 3, height = 3)
}

# use_glint=TRUE
# plot_rocs(all_sets_pint_summaries, all_sets_pint_hierarchy_smrys, all_sets_pint_dedup_smrys, all_sets_pint_unbiased_smrys, all_sets_glint_summaries, all_sets_whinter_summaries, use_glint, "all", TRUE, sprintf("plots/rocs/all_sets_comparison_roc.pdf"))
# plot_rocs(all_sets_pint_summaries, all_sets_pint_hierarchy_smrys, all_sets_pint_dedup_smrys, all_sets_pint_unbiased_smrys, all_sets_glint_summaries, all_sets_whinter_summaries, use_glint, "main", TRUE, sprintf("plots/rocs/main_comparison_roc.pdf"))
# plot_rocs(all_sets_pint_summaries, all_sets_pint_hierarchy_smrys, all_sets_pint_dedup_smrys, all_sets_pint_unbiased_smrys, all_sets_glint_summaries, all_sets_whinter_summaries, use_glint, "interaction", TRUE, sprintf("plots/rocs/int_comparison_roc.pdf"))
#
# plot_rocs(all_sets_pint_summaries, all_sets_pint_hierarchy_smrys, all_sets_pint_dedup_smrys, all_sets_pint_unbiased_smrys, all_sets_glint_summaries, all_sets_whinter_summaries, use_glint, "all", FALSE, sprintf("plots/rocs/all_sets_comparison_roc_FoundOnly.pdf"))
# plot_rocs(all_sets_pint_summaries, all_sets_pint_hierarchy_smrys, all_sets_pint_dedup_smrys, all_sets_pint_unbiased_smrys, all_sets_glint_summaries, all_sets_whinter_summaries, use_glint, "main", FALSE, sprintf("plots/rocs/main_comparison_roc_FoundOnly.pdf"))
# plot_rocs(all_sets_pint_summaries, all_sets_pint_hierarchy_smrys, all_sets_pint_dedup_smrys, all_sets_pint_unbiased_smrys, all_sets_glint_summaries, all_sets_whinter_summaries, use_glint, "interaction", FALSE, sprintf("plots/rocs/int_comparison_roc_FoundOnly.pdf"))
#
## plot times
# melted_times <- melt(all_sets_times, value.name = "Time", variable.name = "Method")
# time_plot <- ggplot(melted_times, aes(x = Method, y = Time)) +
#    geom_boxplot() +
#    theme_bw() +
#    scale_y_continuous(trans = "log2") +
#    ylab("Time (s)") +
#    expand_limits(y = 0)
# ggsave(time_plot, file="plots/all_bench_times.pdf", width = 6, height = 4)
#
# write.csv(summary(all_sets_times), file=sprintf("times_summary_all_sets.csv", bench_set))

# overall summary
all_sets_pint_summaries <- c()
all_sets_pint_hierarchy_summaries <- c()
#all_sets_pint_dedup_summaries <- c()
#all_sets_pint_unbiased_summaries <- c()
all_sets_whinter_summaries <- c()
all_sets_glint_summaries <- c()
all_sets_times <- c()
summarise_bench_sets <- c("8k_only", "wide_10k", "simulated_small_data_sample")
#foreach (bench_set=summarise_bench_sets) %dor% {
#    }
#library(future)
for (bench_set in summarise_bench_sets) {
    gc()
    use_glint <- TRUE
    if (bench_set == "3way") {
        use_pint_pair <- TRUE
    } else {
        use_pint_pair <- FALSE
        all_pint_pair_smrys <- NA
    }

    dir <- paste("whinter_pint_comparison.safe", bench_set, sep = "/")

    rds_files <- list.files(dir, pattern = "all.rds", recursive = TRUE, full.names = TRUE)

    use_values <- c("gene_i", "gene_j", "TP", "strength", "found", "type")

    output <- c()
    output <- foreach(file = rds_files, .combine=cbind) %dopar% {
        print(sprintf("reading file %s", file))
        new_out <- readRDS(file)
        print("done")
        return(new_out)
    }

    output <- data.frame(output)

    all_pint_smrys <- foreach(set = output, .combine = rbind) %do% {
        set$pint_smry |> select(all_of(use_values))
    }
    all_pint_hierarchy_smrys <- foreach(set = output, .combine = rbind) %do% {
        set$pint_hierarchy_smry |> select(all_of(use_values))
    }
    #all_pint_dedup_smrys <- foreach(set = output, .combine = rbind) %do% {
    #    set$pint_dedup_smry
    #}
    #all_pint_unbiased_smrys <- foreach(set = output, .combine = rbind) %do% {
    #    set$pint_unbiased_smry
    #}
    if (use_pint_pair) {
        all_pint_pair_smrys <- foreach(set = output, .combine = rbind) %do% {
            set$pint_pair_smry |> select(all_of(use_values))
        }
    }
    all_whinter_smrys <- foreach(set = output, .combine = rbind) %do% {
        set$whinter_smry |> select(all_of(use_values))
    }
    if (use_glint) {
        all_glint_smrys <- foreach(set = output, .combine = rbind) %do% {
            set$glint_smry |> select(all_of(use_values))
        }
    } else {
        all_glint_smrys <- NA
    }
    fix_na_strengths(all_pint_smrys)
    #fix_na_strengths(all_pint_dedup_smrys)
    #fix_na_strengths(all_pint_unbiased_smrys)
    fix_na_strengths(all_pint_hierarchy_smrys)
    fix_na_strengths(all_whinter_smrys)
    fix_na_strengths(all_glint_smrys)

    all_sets_pint_summaries <- rbind(all_sets_pint_summaries, all_pint_smrys)
    all_sets_pint_hierarchy_summaries <- rbind(all_sets_pint_hierarchy_summaries, all_pint_hierarchy_smrys)
    #all_sets_pint_dedup_summaries <- rbind(all_sets_pint_dedup_summaries, all_pint_dedup_smrys)
    #all_sets_pint_unbiased_summaries <- rbind(all_sets_pint_unbiased_summaries, all_pint_unbiased_smrys)
    all_sets_whinter_summaries <- rbind(all_sets_whinter_summaries, all_whinter_smrys)
    if (use_glint) {
       all_sets_glint_summaries <- rbind(all_sets_glint_summaries, all_glint_smrys)
    }

    #if (use_glint) {
    #    all_times <- foreach(out = output, .combine = rbind) %do% {
    #        data.frame(
    #            pint = out$pint_time,
    #            pint_hierarchy = out$pint_hierarchy_time,
    #            pint_dedup = out$pint_dedup_time,
    #            pint_unbiased = out$pint_unbiased_time,
    #            whinter = out$whinter_time,
    #            glint = out$glint_time
    #        )
    #    }
    #} else {
    #    all_times <- foreach(out = output, .combine = rbind) %do% {
    #        data.frame(
    #            pint = out$pint_time,
    #            pint_hierarchy = out$pint_hierarchy_time,
    #            pint_dedup = out$pint_dedup_time,
    #            pint_unbiased = out$pint_unbiased_time,
    #            whinter = out$whinter_time
    #        )
    #    }
    #}
    all_sets_times <- rbind(all_sets_times, all_times)
    rm(output)
    gc(verbose=FALSE)
}

all_sets_pint_dedup_summaries <- c()
all_sets_pint_unbiased_summaries <- c()
plot_rocs(
    all_sets_pint_summaries, all_sets_pint_hierarchy_summaries, all_sets_pint_dedup_summaries,
    all_sets_pint_unbiased_summaries, all_sets_glint_summaries, all_sets_whinter_summaries, use_glint, "all",
    TRUE, sprintf("plots/rocs/overall_comparison_rocs.pdf")
)
plot_rocs_anyfound(
    all_pint_smrys, all_pint_hierarchy_smrys, all_pint_dedup_smrys,
    all_pint_unbiased_smrys, all_pint_pair_smrys, all_glint_smrys, all_whinter_smrys, use_glint, use_pint_pair, "all",
    sprintf("plots/rocs/overall_comparison_roc_anyfound.pdf")
)
