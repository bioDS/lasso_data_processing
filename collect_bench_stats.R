#!/usr/bin/env Rscript
library(foreach)

# whinter_pint_comparison/simulated_large_data_sample/summaries/n10000_p1000_SNR10_nbi0_nbij500_nlethals0_viol0_19337
library(dplyr)

# bench_sets <- c("1k_whinter_working_data", "simulated_small_data_sample", "simulated_8k", "removed")
bench_sets <- c("simulated_small_data_sample", "8k_only", "wide_only")
# bench_sets <- c("8k_only")
include_main <- TRUE

for (bench_set in bench_sets) {
#bench_set <- "8k_only"
dir <- paste("whinter_pint_comparison", bench_set, sep = "/")

rds_files <- list.files(dir, pattern = "all.rds", recursive = TRUE, full.names = TRUE)

output <- c()
for (file in rds_files) {
  output <- cbind(output, readRDS(file))
}

output <- data.frame(output)

pint_times <- c()
pint_prec <- c()
pint_rec <- c()
whinter_times <- c()
whinter_prec <- c()
whinter_rec <- c()
glint_times <- c()
glint_prec <- c()
glint_rec <- c()
combo_prec <- c()
combo_rec <- c()
combo_total_tp <- c()
whinter_total_tp <- c()
pint_fp <- c()
whinter_fp <- c()
print(sprintf("output contains %d items", length(output)))
#set <- output[1, ]
ind  <-  1
for (set in output) {
  if (include_main) {
    true_effects <- nrow(set$bij) + nrow(set$bi)
  } else  {
    true_effects <- nrow(set$bij)
  }
  print(sprintf("true_effects %d, %d %d \t %s", true_effects, nrow(set$bij), nrow(set$bi), rds_files[ind]))
  ind  <-  ind + 1

  combo <- union(set$pint_fx_int, set$whinter_fx_int) %>%
    filter(type == "interaction") %>%
    select(gene_i, gene_j, type, TP) %>%
    unique()

  if (include_main) {
    pint_ints <- set$pint_smry
    whinter_ints <- set$whinter_smry
  } else  {
    pint_ints <- set$pint_smry %>% filter(type == "interaction")
    whinter_ints <- set$whinter_smry %>% filter(type == "interaction")
  }

  combo_tp <- nrow(combo %>% filter(TP == TRUE))
  combo_total <- nrow(combo)
  pint_tp <- nrow(pint_ints %>% filter(TP == TRUE))
  pint_fp <- c(pint_fp, nrow(pint_ints %>% filter(TP == FALSE)))
  pint_total <- nrow(pint_ints)
  whinter_tp <- nrow(whinter_ints %>% filter(TP == TRUE))
  whinter_fp <- c(whinter_fp, nrow(whinter_ints %>% filter(TP == FALSE)))
  whinter_total <- nrow(whinter_ints)

  combo_total_tp <- c(combo_total_tp, combo_tp)
  whinter_total_tp <- c(whinter_total_tp, whinter_tp)

  pint_times <- cbind(pint_times, set$pint_time[3])
  pint_prec <- cbind(pint_prec, pint_tp / pint_total)
  pint_rec <- cbind(pint_rec, pint_tp / true_effects)

  whinter_times <- cbind(whinter_times, set$whinter_time[3])
  whinter_prec <- cbind(whinter_prec, whinter_tp / whinter_total)
  whinter_rec <- cbind(whinter_rec, whinter_tp / true_effects)

  combo_prec <- cbind(combo_prec, combo_tp / combo_total)
  combo_rec <- cbind(combo_rec, combo_tp / true_effects)

  if (length(set$glint_smry) > 0) {
    if (include_main) {
        glint_ints <- set$glint_smry
    } else {
        glint_ints <- set$glint_smry %>% filter(type == "interaction")
    }
    glint_tp <- nrow(glint_ints %>% filter(TP == TRUE))
    glint_total <- nrow(glint_ints)
    print(sprintf("glint %d:%d %d", glint_tp, glint_total, true_effects))
    glint_times <- cbind(glint_times, set$glint_time[3])
    glint_prec <- cbind(glint_prec, glint_tp / glint_total)
    glint_rec <- cbind(glint_rec, glint_tp / true_effects)
  }
  print(sprintf("whinter %d:%d %d", whinter_tp, whinter_total, true_effects))
  print(sprintf("pint %d:%d %d", pint_tp, pint_total, true_effects))
}

pint_big_smry <- foreach(set=output, .combine=rbind) %do% {
  pint_smry <- set$pint_smry
  ## large <- pint_smry %>% filter(abs(coef.est) > mean(abs(coef.est) + 2*sd(abs(coef.est))))
  ## large <- pint_smry %>% arrange(abs(pval)) %>% head(n=10)
  large <- pint_smry %>% arrange(-abs(coef.est)) %>% head(n=10)
  return(large)
}

pint_times <- pint_times[1, ]
pint_prec <- pint_prec[1, ]
pint_rec <- pint_rec[1, ]
whinter_times <- whinter_times[1, ]
whinter_prec <- whinter_prec[1, ]
whinter_rec <- whinter_rec[1, ]
if (length(glint_times) > 0) {
  glint_times <- glint_times[1, ]
  glint_prec <- glint_prec[1, ]
  glint_rec <- glint_rec[1, ]
}
combo_prec <- combo_prec[1, ]
combo_rec <- combo_rec[1, ]

combo_advantage <- combo_total_tp - whinter_total_tp

print(sprintf("saving to %s-%s", bench_set, "bench_summary.rds"))
saveRDS(
  file = paste(bench_set, "bench_summary.rds", sep = "-"),
  list(
    whinter_times = whinter_times,
    whinter_prec = whinter_prec,
    whinter_rec = whinter_rec,
    whinter_fp = whinter_fp,
    pint_times = pint_times,
    pint_prec = pint_prec,
    pint_rec = pint_rec,
    pint_fp = pint_fp,
    glint_times = glint_times,
    glint_prec = glint_prec,
    glint_rec = glint_rec,
    combo_prec = combo_prec,
    combo_rec = combo_rec,
    combo_advantage = combo_advantage
  )
)
# print("poor sets for pint:")
# print(rds_files[pint_rec <= 0.25])
# print("poor sets for whint:")
# print(rds_files[whinter_prec <= 0.05])
# print("whinter precision on pint poor cases")
# print(whinter_prec[pint_rec <= 0.25])
}
