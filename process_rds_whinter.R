#!/usr/bin/env Rscript
require(dplyr)
require(Pint)
require(glinternet)
require(Matrix)
require(torch)
library(digest)
library(doParallel)

no_cores <- detectCores()
cl <- makeCluster(no_cores)
registerDoParallel(cl)

args <- commandArgs(trailingOnly = TRUE)
file <- as.character(args[1])
output_dir <- as.character(args[2])
methods <- as.character(args[3])
num_features <- as.numeric(args[4])
depth <- as.numeric(args[5])
num_threads <- as.numeric(args[6])
print(sprintf("using %d features", num_features))

#file="~/work/data/simulated_rerun/3way//n4000_p400_nbi10_nbij100_nbijk1000_nlethals0_viol100_snr5_10948.rds"
#output_dir="whinter_pint_comparison/3way//summaries/n4000_p400_nbi10_nbij100_nbijk1000_nlethals0_viol100_snr5_10948"
#methods="noglint"
#num_features=5000
#depth=3
#num_threads=48
#
#file="/home/kel63/work/data/simulated_rerun/8k_only//n8000_p4000_nbi40_nbij200_nbijk0_nlethals0_viol100_snr5_12924.rds"
#output_dir="whinter_pint_comparison/8k_only//summaries/n8000_p4000_nbi40_nbij200_nbijk0_nlethals0_viol100_snr5_12924"
#methods="all"
#num_features=5000
#depth=2
#num_threads=96

#file="/home/kel63/work/data/simulated_rerun/simulated_small_data_sample/n1000_p100_nbi10_nbij50_nbijk0_nlethals0_viol100_snr5_10918.rds"
#output_dir="whinter_pint_comparison/simulated_small_data_sample/summaries/n1000_p100_nbi10_nbij50_nbijk0_nlethals0_viol100_snr5_10918.rds"
#methods="all"
#num_features=1000
#depth=2
#num_threads=1

lethal_coef <- -1000
lambda_min_ratio <- 0.01 # for glinternet only
max_lambdas = 200

whinter_bin <- "~/work/lasso_data_processing/WHInter/src/train_WHInter"
convert_format <- "~/work/lasso_data_processing/convert_format.py"

cat("reading from: ", file, "\n")
cat("writing to dir: ", output_dir, "\n")

# if (!file.exists(paste(output_dir,"converted.tsv", sep='/'))) {
data <- readRDS(file)

X <- data$X
Y <- data$Y
obs <- data$obs
bi_ind <- data$bi_ind
bij_ind <- data$bij_ind
bijk_ind <- data$bijk_ind
lethal_ind <- data$lethal_ind
p <- ncol(X)

if (is.na(sum(Y))) {
  print("Y contains NA, skipping this set")
  q()
}

if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive=TRUE)
}
setwd(output_dir)

if (!file.exists("converted.tsv")) {
  # write tsv for whinter
  x_filename <- "X.csv"
  y_filename <- "Y.csv"

  cat("writing to ", output_dir, "{", x_filename, ",", y_filename, "}\n")

  write.table(X, file = x_filename, sep = ",", col.names = FALSE)
  write.table(Y, file = y_filename, sep = ",", col.names = FALSE)

  print("converting to tsv")
  system(sprintf("%s %s %s converted.tsv", convert_format, x_filename, y_filename))
}

print("running whinter")
command <- sprintf("%s -pathResults ./results/ -useBias 0 -lambdaMinRatio 0.01 -nlambda %d -eps 0.01 -maxSelectedFeatures %d converted.tsv", whinter_bin, max_lambdas, num_features)
cat(getwd(), ":", command, "\n")

if (!dir.exists("results")) {
  dir.create("results")
}
whinter_time <- system.time(
  system(command)
)

model <- list.files("results", pattern = "*model.csv")

print(sprintf("reading from output file %s\n", model))
output <- read.delim(paste("results", model, sep = "/"), sep = ",")

output_last <- output %>% filter(X0 == max(output[1]))

names(output_last) <- c("X0", "lambda", "gene_i", "gene_j", "strength")
whinter_effects <- output_last %>% select(gene_i, gene_j, strength)
whinter_effects$gene_i <- whinter_effects$gene_i + 1
whinter_effects$gene_j <- whinter_effects$gene_j + 1


print("checking accuracy of results")

## handy functions

get_3way_hash <- function(i, j, k, X) {
  return(digest::sha1(X[,i] * X[,j] * X[,k]))
}

tp_main_hashes <- foreach(i=1:nrow(bi_ind), .combine=rbind) %do% {
  get_3way_hash(bi_ind[i,]$gene_i, bi_ind[i,]$gene_i, bi_ind[i,]$gene_i, X)
}
tp_pair_hashes <- foreach(i=1:nrow(bij_ind), .combine=rbind) %do% {
  get_3way_hash(bij_ind[i,]$gene_i, bij_ind[i,]$gene_j, bij_ind[i,]$gene_i, X)
}
if (depth == 3) {
  tp_trip_hashes <- foreach(i=1:nrow(bijk_ind), .combine=rbind) %do% {
    get_3way_hash(bijk_ind[i,]$gene_i, bijk_ind[i,]$gene_j, bijk_ind[i,]$gene_k, X)
  }
  tp_hashes <- rbind(tp_main_hashes, tp_pair_hashes, tp_trip_hashes)
  rm(tp_trip_hashes)
} else {
  tp_hashes <- rbind(tp_main_hashes, tp_pair_hashes)
}
rm(tp_main_hashes)
rm(tp_pair_hashes)

get_equiv_tp <- function(fx_set, X, tp_hashes, depth) {
  found_set <- fx_set |> filter(found==TRUE)
  i_inds <- found_set$gene_i
  if (is.na(found_set$gene_j[1])) {
    j_inds <- i_inds
  } else {
    j_inds <- found_set$gene_j
  }
  if (depth < 3) {
    k_inds <- i_inds
  }else {
    if (is.na(found_set$gene_k[1])) {
      k_inds <- i_inds
    } else {
      k_inds <- found_set$gene_k
    }
  }

  hashes <- foreach(i=1:nrow(found_set), .combine=rbind) %do% {
    get_3way_hash(i_inds[i], j_inds[i], k_inds[i], X)
  }
  fx_set$equiv_tp <- fx_set$TP

  equiv_tp <- hashes %in% tp_hashes

  fx_set[fx_set$found == TRUE,]$equiv_tp <- equiv_tp
  return(fx_set$equiv_tp)
}

## WHInter ********************************************************************************************************


all_i <- data.frame(gene_i = 1:p) |> mutate(gene_j = NA)

whinter_fx_main <- whinter_effects %>%
  filter(gene_i == gene_j) %>%
  mutate(gene_j = NA) %>%
  mutate(found = TRUE) %>%
  full_join(bi_ind |> select(gene_i, coef), by = "gene_i") %>%
  mutate(TP = !is.na(coef)) %>%
  mutate(found = !is.na(strength)) %>%
  mutate(lethal = gene_i %in% lethal_ind[["gene_i"]]) %>%
  full_join(all_i, by=c("gene_i", "gene_j")) %>%
  mutate(type = "main") %>%
  select(gene_i, gene_j, type, lethal, found, TP, strength)
  whinter_fx_main[is.na(whinter_fx_main$strength),]$strength <- 0
  whinter_fx_main[is.na(whinter_fx_main$TP),]$TP <- FALSE
  whinter_fx_main[is.na(whinter_fx_main$found),]$found <- FALSE
  whinter_fx_main[is.na(whinter_fx_main$lethal),]$lethal <- FALSE
whinter_fx_main$gene_j <- NA
whinter_fx_main$equiv_tp <- get_equiv_tp(whinter_fx_main, X, tp_hashes, 2)

whinter_fx_int  <- whinter_effects %>%
  filter(gene_i != gene_j)

if (nrow(whinter_fx_int) > 0) {
  all_ij <- as_array(torch_combinations(torch_tensor(1:p), 2))
  all_ij <- data.frame(gene_i = all_ij[,1], gene_j = all_ij[,2])
  whinter_fx_int <- whinter_fx_int %>% mutate(found = TRUE) %>%
    full_join(bij_ind, by = c("gene_i","gene_j")) %>%
    mutate(TP = !is.na(coef)) %>%
    mutate(found = !is.na(strength)) %>%
    mutate(lethal = gene_i %in% lethal_ind[["gene_i"]]) %>%
    full_join(all_ij, by=c("gene_i", "gene_j")) %>%
    mutate(type = "interaction") %>%
    select(gene_i, gene_j, type, lethal, found, TP, strength)
  whinter_fx_int[is.na(whinter_fx_int$strength),]$strength <- 0
  whinter_fx_int[is.na(whinter_fx_int$TP),]$TP <- FALSE
  whinter_fx_int[is.na(whinter_fx_int$found),]$found <- FALSE
  whinter_fx_int[is.na(whinter_fx_int$lethal),]$lethal <- FALSE
  whinter_fx_int$equiv_tp <- get_equiv_tp(whinter_fx_int, X, tp_hashes, 2)
} else {
  whinter_fx_int <- NA
}
#found_fx_main <- whinter_fx_main |> filter(found==TRUE)
#found_fx_int <- whinter_fx_int |> filter(found==TRUE)
#Z <- cbind(X[, (found_fx_main)[["gene_i"]]])
#if (nrow(found_fx_int) > 0) {
#  for (i in 1:nrow(found_fx_int)) {
#    Z <- cbind(Z, X[, found_fx_int[i, ][["gene_i"]], drop = FALSE] * X[, found_fx_int[i, ][["gene_j"]], drop = FALSE])
#  }
#}
#Z <- as.matrix(Z)
#colnames(Z) <- rownames(Z) <- NULL
#Ynum <- as.numeric(Y)
#ols_time <- system.time(fit_red <- lm(Ynum ~ Z))

#pvals <- data.frame(id = 1:ncol(Z), coef = coef(fit_red)[-1]) %>%
#  filter(!is.na(coef)) %>%
#  data.frame(., pval = summary(fit_red)$coef[-1, 4]) %>%
#  tbl_df()
#whinter_smry <- left_join(rbind(whinter_fx_main, whinter_fx_int) %>% data.frame(id = 1:nrow(.), .), pvals, by = "id") %>%
#  mutate(pval = ifelse(is.na(pval), 1, pval)) %>%
#  rename(coef.est = coef) #%>%
  ## left_join(., obs, by = c("gene_i", "gene_j"))
whinter_smry <- rbind(whinter_fx_main, whinter_fx_int) %>% data.frame(id = 1:nrow(.), .)


saveRDS(list(whinter_fx_int = whinter_fx_int, time = time, smry = whinter_smry), "whinter_fx_int.rds")

print("running Pint for comparison")

## Pint ********************************************************************************************************

get_pint_smry <- function(fit, depth) {
  all_i$gene_j <- NA
  pint_fx_main <- data.frame(gene_i = as.numeric(fit$main$effects$i), strength = fit$main$effects$strength) %>%
    arrange(gene_i) %>%
    mutate(gene_j = NA) %>%
    full_join(bi_ind, by = "gene_i") %>%
    mutate(TP = !is.na(coef)) %>%
    mutate(found = !is.na(strength)) %>%
    mutate(lethal = gene_i %in% lethal_ind[["gene_i"]]) %>%
    select(gene_i, gene_j, TP, found, lethal, strength) %>%
    arrange(desc(TP)) %>%
    arrange(desc(lethal)) %>%
    full_join(all_i, by=c("gene_i", "gene_j")) %>%
    mutate(type = "main") %>%
    tbl_df()
  pint_fx_main[is.na(pint_fx_main$strength),]$strength <- 0
  pint_fx_main[is.na(pint_fx_main$TP),]$TP <- FALSE
  pint_fx_main[is.na(pint_fx_main$found),]$found <- FALSE
  pint_fx_main[is.na(pint_fx_main$lethal),]$lethal <- FALSE
  pint_fx_main$equiv_tp <- get_equiv_tp(pint_fx_main, X, tp_hashes, 2)
  if (length(fit$pairwise$effects$strength) > 0) {
    pint_fx_int <- data.frame(
      gene_i = as.numeric(fit$pairwise$effects$i), gene_j = as.numeric(fit$pairwise$effects$j),
      strength = fit$pairwise$effects$strength %>% unlist()
    ) %>%
      arrange(gene_i) %>%
      ## left_join(., obs, by = c("gene_i", "gene_j")) %>%
      rowwise() %>%
      full_join(., rbind(bij_ind, lethal_ind), by = c("gene_i", "gene_j")) %>%
      ungroup() %>%
      mutate(TP = !is.na(coef)) %>%
      mutate(found = !is.na(strength)) %>%
      mutate(lethal = (coef == lethal_coef)) %>%
      arrange(desc(TP)) %>%
      arrange(desc(lethal)) %>%
      select(gene_i, gene_j, TP, found, lethal, strength) %>%
      distinct(gene_i, gene_j, .keep_all = TRUE) %>%
      full_join(all_ij, by=c("gene_i", "gene_j")) %>%
      mutate(type = "interaction") %>%
      tbl_df()
    pint_fx_int[is.na(pint_fx_int$strength),]$strength <- 0
    pint_fx_int[is.na(pint_fx_int$TP),]$TP <- FALSE
    pint_fx_int[is.na(pint_fx_int$found),]$found <- FALSE
    pint_fx_int[is.na(pint_fx_int$lethal),]$lethal <- FALSE
    pint_fx_int$equiv_tp <- get_equiv_tp(pint_fx_int, X, tp_hashes, 2)
  } else {
    pint_fx_int <- NA
  }
  if (depth == 3) {
    # three-way
    bijk_ind <- data$bijk_ind
    all_ijk <- as_array(torch_combinations(torch_tensor(1:p), 3))
    all_ijk <- data.frame(gene_i = all_ijk[,1], gene_j = all_ijk[,2], gene_k = all_ijk[,3])
    pint_fx_trip <- data.frame(gene_i = as.numeric(fit$triple$effects$i),
                          gene_j = as.numeric(fit$triple$effects$j),
                          gene_k = as.numeric(fit$triple$effects$k),
                          strength = fit$triple$effects$strength)
    hashes <- foreach(i=1:nrow(pint_fx_trip), .combine=rbind) %do% get_3way_hash(pint_fx_trip[i,]$gene_i, pint_fx_trip[i,]$gene_j, pint_fx_trip[i,]$gene_k, X)
    pint_fx_trip$hash <- hashes
    pint_fx_trip$equiv_tp <- pint_fx_trip$hash %in% tp_hashes
    pint_fx_trip <- full_join(pint_fx_trip, all_ijk, by=c("gene_i", "gene_j", "gene_k")) |>
      full_join(bijk_ind, by=c("gene_i", "gene_j", "gene_k")) |>
      mutate(found = !is.na(strength), TP = !is.na(coef)) |>
      mutate(lethal = (coef == lethal_coef), type="3way") |>
      select(gene_i, gene_j, gene_k, TP, found, lethal, strength, type, equiv_tp)
    pint_fx_trip[is.na(pint_fx_trip$equiv_tp),]$equiv_tp <- FALSE
    pint_fx_main$gene_k <- NA
    pint_fx_int$gene_k <- NA
    pint_smry <- rbind(pint_fx_main, pint_fx_int, pint_fx_trip) %>% data.frame(id = 1:nrow(.), .)
  } else {
    pint_smry <- rbind(pint_fx_main, pint_fx_int) %>% data.frame(id = 1:nrow(.), .)
  }
  return(pint_smry)
}


pint_time <- system.time(fit <- interaction_lasso(X, Y, max_nz_beta = num_features, max_lambdas = max_lambdas, depth = depth, num_threads = num_threads))
pint_smry <- get_pint_smry(fit, depth)

## Pint w/ hierarchy assumption
pint_hierarchy_time <- system.time(fit <- interaction_lasso(X, Y, max_nz_beta = num_features, max_lambdas = max_lambdas, depth = depth, num_threads = num_threads, approximate_hierarchy = TRUE))
pint_hierarchy_smry <- get_pint_smry(fit, depth)

## Pint w/ duplicate removal
pint_dedup_time <- system.time(fit <- interaction_lasso(X, Y, max_nz_beta = num_features, max_lambdas = max_lambdas, depth = depth, num_threads = num_threads, check_duplicates = TRUE))
pint_dedup_smry <- get_pint_smry(fit, depth)

## Pint w/ estimate unbiased
pint_unbiased_time <- system.time(fit <- interaction_lasso(X, Y, max_nz_beta = num_features, max_lambdas = max_lambdas, depth = depth, num_threads = num_threads, estimate_unbiased = TRUE))
pint_unbiased_smry <- get_pint_smry(fit$estimate_unbiased, depth)

if (depth == 3) {
  pint_pair_time <- system.time(fit <- interaction_lasso(X, Y, max_nz_beta = num_features, max_lambdas = max_lambdas, depth = 2, num_threads = num_threads))
  pint_pair_smry <- get_pint_smry(fit, 2)
}


## glinternet ********************************************************************************************************

if (methods == "noglint") {
  print("skipping glinternet")
  glint_smry <- c()
  glint_time <- c()
  glint_fx_int <- c()
  glint_fx_main <- c()
} else {
  print("glinternet for comparison")

  glint_time <- system.time(fit <- glinternet(
    X = X %>% as.matrix(),
    Y = Y %>% as.numeric(),
    numLevels = rep(1, p),
    family = "gaussian",
    numToFind = num_features,
    nLambda = max_lambdas, numCores = num_threads, lambdaMinRatio = lambda_min_ratio, verbose = TRUE
  ))

  cf <- coef(fit, lambdaType = "lambdaHat") # lambdaIndex = 50)#
  cf <- cf[[length(cf)]]

  ## Collect coefficients
  glint_fx_main <- data.frame(
    gene_i = cf$mainEffects$cont,
    strength = cf$mainEffectsCoef$cont %>% lapply(., function(x) x[[1]]) %>% unlist()
  ) %>%
    arrange(gene_i) %>%
    mutate(gene_j = NA) %>%
    full_join(bi_ind, by = "gene_i") %>%
    mutate(TP = !is.na(coef)) %>%
    mutate(found = !is.na(strength)) %>%
    mutate(lethal = gene_i %in% lethal_ind[["gene_i"]]) %>%
    select(gene_i, gene_j, TP, found, lethal, strength) %>%
    arrange(desc(TP)) %>%
    arrange(desc(lethal)) %>%
    full_join(all_i, by=c("gene_i", "gene_j")) %>%
    mutate(type = "main") %>%
    tbl_df()
  glint_fx_main[is.na(glint_fx_main$strength),]$strength <- 0
  glint_fx_main[is.na(glint_fx_main$TP),]$TP <- FALSE
  glint_fx_main[is.na(glint_fx_main$found),]$found <- FALSE
  glint_fx_main[is.na(glint_fx_main$lethal),]$lethal <- FALSE
  glint_fx_main$equiv_tp <- get_equiv_tp(glint_fx_main, X, tp_hashes, 2)

  glint_fx_int <- data.frame(
    gene_i = cf$interactions$contcont[, 1], gene_j = cf$interactions$contcont[, 2],
    strength = cf$interactionsCoef$contcont %>% unlist()
  ) %>%
    arrange(gene_i) %>%
    ## left_join(., obs, by = c("gene_i", "gene_j")) %>%
    rowwise() %>%
    full_join(., rbind(bij_ind, lethal_ind), by = c("gene_i", "gene_j")) %>%
    ungroup() %>%
    mutate(TP = !is.na(coef)) %>%
    mutate(found = !is.na(strength)) %>%
    mutate(lethal = (coef == lethal_coef)) %>%
    select(gene_i, gene_j, TP, found, lethal, strength) %>%
    arrange(desc(TP)) %>%
    arrange(desc(lethal)) %>%
    full_join(all_ij, by=c("gene_i", "gene_j")) %>%
    mutate(type = "interaction") %>%
    tbl_df()
  glint_fx_int[is.na(glint_fx_int$strength),]$strength <- 0
  glint_fx_int[is.na(glint_fx_int$TP),]$TP <- FALSE
  glint_fx_int[is.na(glint_fx_int$found),]$found <- FALSE
  glint_fx_int[is.na(glint_fx_int$lethal),]$lethal <- FALSE
  glint_fx_int$equiv_tp <- get_equiv_tp(glint_fx_int, X, tp_hashes, 2)


  #found_fx_main <- glint_fx_main |> filter(found==TRUE)
  #found_fx_int <- glint_fx_int |> filter(found==TRUE)
  ### Statistical test if b_i and b_ij are sig. > 0
  #Z <- cbind(X[, found_fx_main[["gene_i"]]])
  #for (i in 1:nrow(found_fx_int)) {
  #  Z <- cbind(Z, X[, found_fx_int[i, ][["gene_i"]], drop = FALSE] * X[, found_fx_int[i, ][["gene_j"]], drop = FALSE])
  #}
  #Z <- as.matrix(Z)
  #colnames(Z) <- rownames(Z) <- NULL
  #Ynum <- as.numeric(Y)
  #fit_red <- lm(Ynum ~ Z)


  #glint_pvals <- data.frame(id = 1:ncol(Z), coef = coef(fit_red)[-1]) %>%
  #  filter(!is.na(coef)) %>%
  #  data.frame(., pval = summary(fit_red)$coef[-1, 4]) %>%
  #  tbl_df()

  #glint_smry <- left_join(rbind(glint_fx_main, glint_fx_int) %>% data.frame(id = 1:nrow(.), .), glint_pvals, by = "id") %>%
  #  mutate(pval = ifelse(is.na(pval), 1, pval)) %>%
  #  rename(coef.est = coef) #%>%
  #  ## left_join(., obs, by = c("gene_i", "gene_j"))
  glint_smry <- rbind(glint_fx_main, glint_fx_int) %>% data.frame(id = 1:nrow(.), .)
}


## summary ***********************************************************************************************************

#print("whinter")
#whinter_fx_int %>% data.frame()
#print("pint")
#pint_fx_int %>% data.frame()

saveRDS(list(
  whinter_time = whinter_time,
  whinter_smry = whinter_smry,
  pint_time = pint_time,
  pint_unbiased_time = pint_unbiased_time,
  pint_dedup_time = pint_dedup_time,
  pint_hierarchy_time = pint_hierarchy_time,
  pint_smry = pint_smry,
  pint_hierarchy_smry = pint_hierarchy_smry,
  pint_dedup_smry = pint_dedup_smry,
  pint_unbiased_smry = pint_unbiased_smry,
  glint_time = glint_time,
  glint_smry = glint_smry,
  bij = bij_ind,
  bi = bi_ind
), "all.rds")

print("summarising")

print(sprintf("there were %d true interactions\n", nrow(bij_ind)))
print(sprintf("Pint suggested %d in %.3fs, %d of which were correct", nrow(pint_smry |> filter(found == TRUE)), pint_time[3], nrow(pint_smry %>% filter(TP == TRUE) |> filter(found == TRUE))))
print(sprintf("Pint (hierarchy) suggested %d in %.3fs, %d of which were correct", nrow(pint_hierarchy_smry |> filter(found == TRUE)), pint_hierarchy_time[3], nrow(pint_hierarchy_smry %>% filter(TP == TRUE) |> filter(found == TRUE))))
print(sprintf("Pint (dedup) suggested %d in %.3fs, %d of which were correct", nrow(pint_dedup_smry |> filter(found == TRUE)), pint_dedup_time[3], nrow(pint_dedup_smry %>% filter(TP == TRUE) |> filter(found == TRUE))))
print(sprintf("Pint (unbiased) suggested %d in %.3fs, %d of which were correct", nrow(pint_unbiased_smry |> filter(found == TRUE)), pint_unbiased_time[3], nrow(pint_unbiased_smry %>% filter(TP == TRUE) |> filter(found == TRUE))))
print(sprintf("Whinter suggested %d in %.3fs, %d of which were correct", nrow(whinter_smry |> filter(found == TRUE)), whinter_time[3], nrow(whinter_smry  |> filter(found == TRUE) %>% filter(TP == TRUE))))
if (methods != "noglint") {
  print(sprintf("glinternet suggested %d in %.3fs, %d of which were correct", nrow(glint_smry |> filter(found == TRUE)), glint_time[3], nrow(glint_smry  |> filter(found == TRUE) %>% filter(TP == TRUE))))
}
if (depth == 3) {
  print(sprintf("Pint (pairwise) suggested %d in %.3fs, %d of which were correct", nrow(pint_pair_smry |> filter(found == TRUE)), pint_pair_time[3], nrow(pint_pair_smry %>% filter(TP == TRUE) |> filter(found == TRUE))))
}
