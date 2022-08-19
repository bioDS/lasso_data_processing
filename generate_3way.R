#!/usr/bin/env Rscript
library(foreach)
library(doMC)

source("generation_functions.R")

registerDoMC(cores=detectCores())

foreach (SNR=c(2, 4, 8)) %:%
    foreach (p=c(100, 200)) %dopar% {
        n <- 10*p;

        p_2 <- p*(p+1)/2
        p_3 <- p*(p+1)*(p-1)/6

        num_bi <- 10
        num_bij <- 20
        num_bijk <- 20
        ## num_bi <-    as.integer(p/10)
        #num_bi <-   as.integer(p/10)
        #num_bij <-  as.integer(p_2/20)
        #num_bijk <- as.integer(p_3/100)
        num_lethals=0

        set = generate_set(n, p, SNR, num_bi, num_bij, num_bijk, num_lethals)

        saveRDS(set, file = sprintf("~/work/data/simulated_rerun/3way/n%d_p%d_SNR%d_nbi%d_nbij%d_nbijk%d_nlethals%d_%d.rds",
                           n, p, SNR, num_bi, num_bij, num_bijk, num_lethals, (runif(1) * 1e5) %>% floor))
    }
