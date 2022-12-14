#!/usr/bin/env Rscript
library(Pint)
library(pROC)
library(dplyr)
library(foreach)
library(doMC)
library(ggplot2)

registerDoMC(cores=detectCores())

source("generation_functions.R")
source("summary_functions.R")

path <- "./data/simulated_rerun/3way/"
ensure_small_set_exists(path)
all_effects = get_all_effects_for_dir(path, 3, -1, approximate_hierarchy = FALSE)

main_effects = all_effects %>% filter(type=="main")
main_roc <- roc(response = main_effects$TP, predict = abs(main_effects$strength))
save_roc(main_roc, "plots/roc/resim_1way_roc.pdf")

# roc including effects that were not found
pair_effects = all_effects %>% filter(type=="interaction")
pair_roc <- roc(response = pair_effects$TP, predict = abs(pair_effects$strength))
save_roc(pair_roc, "plots/roc/resim_2way_roc.pdf")

trip_effects = all_effects %>% filter(type=="3way")
trip_roc <- roc(response = trip_effects$TP, predict = abs(trip_effects$strength))
save_roc(trip_roc, "plots/roc/resim_3way_roc.pdf")

all_roc <- roc(response = all_effects$TP, predict = abs(all_effects$strength))
save_roc(all_roc, "plots/roc/resim_all_roc.pdf")

# roc summarising found effects only (i.e. ignoring false negatives)
main_effects_found = main_effects %>% filter(found==TRUE)
main_found_roc <- roc(response = main_effects_found$TP, predict = abs(main_effects_found$strength))
save_roc(main_found_roc, "plots/roc/resim_1way_foundonly_roc.pdf")

pair_effects_found = pair_effects %>% filter(found==TRUE)
pair_found_roc <- roc(response = pair_effects_found$TP, predict = abs(pair_effects_found$strength))
save_roc(pair_found_roc, "plots/roc/resim_2way_foundonly_roc.pdf")

trip_effects_found = trip_effects %>% filter(found==TRUE)
trip_found_roc <- roc(response = trip_effects_found$TP, predict = abs(trip_effects_found$strength))
save_roc(trip_found_roc, "plots/roc/resim_3way_foundonly_roc.pdf")

all_found <- all_effects %>% filter(found == TRUE)
all_roc <- roc(response = all_found$TP, predict = abs(all_found$strength))
save_roc(all_roc, "plots/roc/resim_all_found_roc.pdf")
