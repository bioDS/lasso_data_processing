#!/usr/bin/env Rscript
library(Pint)
library(pROC)
library(dplyr)
library(foreach)
library(doMC)
library(ggplot2)

registerDoMC(cores=detectCores())

source("generation_functions.R")

generate_sets <- function(path, num, size = "large", threeway = FALSE) {
    viol = 100
    num_lethals = 0
    if (size == "8k") {
        n = 8000
        p = 4000
        snr = 5
        num_bi = 40
        num_bij = 200
        num_bijk = 0
    } else if (size == "p100") {
        n = 1000
        p = 100
        snr = 5
        num_bi = 10
        num_bij = 50
        num_bijk = 0
    } else if (size == "wide") {
        n = 1000
        p = 20000
        snr = 5
        num_bi = 100
        num_bij = 500
        num_bijk = 0
    } else if (size == "large_3way") {
        n = 40000
        p = 4000
        snr = 5
        num_bi = 10
        num_bij = 100
        num_bijk = 1000
    } else if (size == "small_3way") {
        n = 4000
        p = 400
        snr = 5
        num_bi = 10
        num_bij = 100
        num_bijk = 1000
    } else if (size == "wide_3way") {
        n = 1000
        p = 10000
        snr = 5
        num_bi = 10
        num_bij = 100
        num_bijk = 1000
    }
    else {
        print("invalid size")
    }
    for (e in 1:num) {
        dataset = generate_set(n, p, snr, num_bi, num_bij, num_bijk, num_lethals)
        saveRDS(dataset, file = paste0(path, sprintf("n%d_p%d_nbi%d_nbij%d_nbijk%d_nlethals%d_viol%d_snr%d_%d.rds",
                        n, p, num_bi, num_bij, num_bijk, num_lethals, viol, snr, (runif(1) * 1e5) %>% floor)))
    }
}

ensure_8k_set_exists <- function(path) {
    if(!dir.exists(path)) {
        dir.create(path, recursive = TRUE)
    }
    existing_files = list.files(path)
    if (length(existing_files) == 0) {
        generate_sets(path, 10, "8k")
    }
}

ensure_p100_set_exists <- function(path) {
    if(!dir.exists(path)) {
        dir.create(path, recursive = TRUE)
    }
    existing_files = list.files(path)
    if (length(existing_files) == 0) {
        generate_sets(path, 50, "p100")
    }
}

ensure_wide_set_exists <- function(path) {
    if(!dir.exists(path)) {
        dir.create(path, recursive = TRUE)
    }
    existing_files = list.files(path)
    if (length(existing_files) == 0) {
        generate_sets(path, 10, "wide")
    }
}

ensure_3way_set_exists <- function(path) {
    if(!dir.exists(path)) {
        dir.create(path, recursive = TRUE)
    }
    existing_files = list.files(path)
    if (length(existing_files) == 0) {
        generate_sets(path, 10, "small_3way", threeway = TRUE)
    }
}

ensure_p100_set_exists("./data/simulated_rerun/simulated_small_data_sample/")
ensure_8k_set_exists("./data/simulated_rerun/8k_only/")
ensure_wide_set_exists("./data/simulated_rerun/wide_only/")
ensure_3way_set_exists("./data/simulated_rerun/3way")
