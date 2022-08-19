#!/usr/bin/env Rscript
source("generation_functions.R")

n <- 4000
p <- 400
snr <- 5
num_bi <- 10
num_bij <- 100
num_bijk <- 0
viol <- 100
num_lethals <- 0

dataset = generate_set(n, p, snr, num_bi, num_bij, num_bijk, num_lethals)

saveRDS(dataset, file = "latest_generator_version.rds")
