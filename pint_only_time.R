#!/usr/bin/Rscript
library(Pint)

num_features <- 500

d <- readRDS("./data/8k_only/n8000_p4000_SNR5_nbi40_nbij200_nlethals0_viol100_31386.rds")
#d <- readRDS("./data/simulated_small_data_sample/n1000_p100_SNR10_nbi0_nbij100_nlethals0_viol0_33859.rds")
X <- d$X
Y <- d$Y
run_times = c()
for (i in 0:5) {
  pint_time <- system.time(
    fit <- interaction_lasso(X, Y, max_nz_beta = num_features, depth = 2)
  )
  run_times = c(run_times, pint_time[3])
}
print(run_times)
pint_time <- median(run_times)
#print(sprintf("time: %s", pint_time))
print(pint_time)
saveRDS(file = "time.rds", pint_time)
