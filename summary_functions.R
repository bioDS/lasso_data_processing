
get_file_effects <- function(file, depth, num_nz, approximate_hierarchy) {
    print(sprintf("running %s", file))

    f <- readRDS(file)
    X <- f$X
    Y <- f$Y
    bi_ind <- f$bi_ind
    bij_ind <- f$bij_ind
    bijk_ind <- f$bijk_ind
    #obs <- f$obs # not present for threeway simulations
    lethal_ind <- f$lethal_ind

    full_lasso_results <- interaction_lasso(X, Y, max_nz_beta = num_nz,
                                            check_duplicates = TRUE, depth = depth,
                                            estimate_unbiased = TRUE, num_threads = 8,
                                            approximate_hierarchy = approximate_hierarchy,
                                            verbose = FALSE)
    ## lasso_results <- full_lasso_results$estimate_unbiased
    lasso_results <- full_lasso_results
    pint_fx_main <- data.frame(gene_i = as.numeric(lasso_results$main$effects$i),
                            unbiased_strength = as.numeric(lasso_results$estimate_unbiased$main$effects$strength),
                            strength = as.numeric(lasso_results$main$effects$strength)) %>%
    arrange(gene_i) %>%
    mutate(type = "main", gene_j = NA, gene_k = NA, TP = (gene_i %in% bi_ind[["gene_i"]])) %>%
    select(gene_i, gene_j, gene_k, type, TP, strength, unbiased_strength) %>%
    arrange(desc(TP)) %>%
    mutate(found=TRUE) %>%
    tbl_df()

    true_effects = bi_ind %>%
      select(gene_i) %>%
      mutate(found=(gene_i %in% pint_fx_main[["gene_i"]]))
    not_found = true_effects %>% filter(found == FALSE) %>%
      mutate(gene_j = NA, gene_k = NA) %>%
      mutate(TP=TRUE, strength=0.0, unbiased_strength=0.0, type="main")

    pint_fx_main <- rbind(pint_fx_main, not_found)
    if (length(lasso_results$pairwise$effects$strength) > 0) {
        pint_fx_int <- data.frame(
            gene_i = as.numeric(lasso_results$pairwise$effects$i),
            gene_j = as.numeric(lasso_results$pairwise$effects$j),
            unbiased_strength = as.numeric(lasso_results$estimate_unbiased$pairwise$effects$strength),
            strength = as.numeric(lasso_results$pairwise$effects$strength)) %>%
        arrange(gene_i, gene_j) %>%
        #left_join(., obs, by = c("gene_i", "gene_j")) %>%
        rowwise() %>%
        full_join(., rbind(bij_ind, lethal_ind), by = c("gene_i", "gene_j")) %>%
        ungroup() %>%
        mutate(TP = !is.na(coef)) %>%
        mutate(type = "interaction", gene_k = NA) %>%
        mutate(found = !is.na(strength)) %>%
        arrange(desc(TP)) %>%
        select(gene_i, gene_j, gene_k, type, TP, strength, unbiased_strength, found) %>%
        distinct(gene_i, gene_j, .keep_all = TRUE) %>%
        tbl_df()
        pint_fx_int[is.na(pint_fx_int$strength),]$strength <- 0


    pint_fx_main <- rbind(pint_fx_main, not_found)

    } else {
        pint_fx_int <- NA
    }
    if (length(lasso_results$triple$effects$strength) > 0) {
        pint_fx_trip <- data.frame(
            gene_i = as.numeric(lasso_results$trip$effects$i),
            gene_j = as.numeric(lasso_results$trip$effects$j),
            gene_k = as.numeric(lasso_results$trip$effects$k),
            unbiased_strength = as.numeric(lasso_results$estimate_unbiased$trip$effects$strength),
            strength = as.numeric(lasso_results$trip$effects$strength)) %>%
        arrange(gene_i, gene_j, gene_k) %>%
        ## left_join(., obs, by = c("gene_i", "gene_j")) %>%
        rowwise() %>%
        full_join(., rbind(bijk_ind, lethal_ind), by = c("gene_i", "gene_j", "gene_k")) %>%
        ungroup() %>%
        mutate(TP = !is.na(coef)) %>%
        mutate(found = !is.na(strength)) %>%
        arrange(desc(TP)) %>%
        mutate(type = "3way") %>%
        select(gene_i, gene_j, gene_k, type, TP, strength, unbiased_strength, found) %>%
        distinct(gene_i, gene_j, .keep_all = TRUE) %>%
        tbl_df()
        pint_fx_trip[is.na(pint_fx_trip$strength),]$strength <- 0
    } else {
        pint_fx_trip <- NA
    }
    if (is.na(pint_fx_trip)) {
        this_run_effects <- rbind(pint_fx_main, pint_fx_int)
    } else {
        this_run_effects <- rbind(pint_fx_main, pint_fx_int, pint_fx_trip)
    }
    # add rows for missed effects (false negatives)
    return(this_run_effects)
}

get_all_effects_for_dir <- function(use_data, depth, num_nz, approximate_hierarchy = TRUE) {
    all_files <- list.files(use_data, pattern = "*.rds")
    all_effects <- foreach(file = all_files, .combine = rbind) %do% {
        get_file_effects(paste0(use_data, "/", file), depth, num_nz, approximate_hierarchy)
    }
    return(data.frame(all_effects))
}

save_roc <- function(roc, filename) {
    roc_plot <- ggroc(roc) +
        theme_bw() +
        geom_segment(x = -1, y = 0, xend = 0, yend = 1, color = "#4c72b0") +
        annotate("text", x = 0.2, y = 0.1, label = sprintf("auc: %.2f", round(roc$auc, digits = 2)))
    ggsave(roc_plot, file = filename, width = 3, height = 3)
    print(roc)
}

generate_sets <- function(path, num, size = "large") {
    viol = 100
    num_lethals = 0
    if (size == "large") {
        n = 40000
        p = 4000
        snr = 5
        num_bi = 10
        num_bij = 100
        num_bijk = 1000
    } else if (size == "small") {
        n = 4000
        p = 400
        snr = 5
        num_bi = 10
        num_bij = 100
        num_bijk = 1000
    } else if (size == "wide") {
        n = 1000
        p = 10000
        snr = 5
        num_bi = 10
        num_bij = 100
        num_bijk = 1000
    }
    else {
        print("size must be either large, small, or wide")
    }
    for (e in 1:num) {
        dataset = generate_set(n, p, snr, num_bi, num_bij, num_bijk, num_lethals)
        saveRDS(dataset, file = paste0(path, sprintf("n%d_p%d_nbi%d_nbij%d_nbijk%d_nlethals%d_viol%d_snr%d_%d.rds",
                        n, p, num_bi, num_bij, num_bijk, num_lethals, viol, snr, (runif(1) * 1e5) %>% floor)))
    }
}

ensure_3way_set_exists <- function(path) {
    if(!dir.exists(path)) {
        dir.create(path, recursive = TRUE)
    }
    existing_files = list.files(path)
    if (length(existing_files) == 0) {
        generate_sets(path, 10, "large")
    }
}

ensure_small_set_exists <- function(path) {
    if(!dir.exists(path)) {
        dir.create(path, recursive = TRUE)
    }
    existing_files = list.files(path)
    if (length(existing_files) == 0) {
        generate_sets(path, 10, "small")
    }
}
