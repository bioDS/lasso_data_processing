#!/usr/bin/Rscript

require(Matrix)
require(Pint)
require(dplyr)
require(normInfectX)
require(biomaRt)
require(stringr)
require(data.table)
require(future.apply)

plan(multicore)

library(progressr)
# handlers("progress", "beepr")

ko_min <- 0
risearch_energy <- -20
cutoff <- 0 # not relevant for RIsearch2, energy level is specified in the command
## data <- "vaccinia"
data <- "mock"
# nzbeta=2000
nzbeta <- 100
depth <- 3
x_binary_file <- sprintf("infx_X_binary_risearch_energy%d_%s.rds", risearch_energy, data)
x_filtered_file <- sprintf("infx_X_binary_risearch_energy%d_%s_filtered.rds", risearch_energy, data)
yfile <- sprintf("infx_Y_%s.rds", data)
xfile <- sprintf("infx_%s_X.rds", data)

intermediate_result_filename <- sprintf("infectX_ifx_only_cutoff%s_%s_nzbeta%s_depth%s_komin%d_risearch2_energy%d.rds", cutoff, data, nzbeta, depth, ko_min, risearch_energy)
result_filename <- sprintf("infectX_ifx_full_stats_cutoff%s_%s_nzbeta%s_depth%s_risearch2_energy%d.rds", cutoff, data, nzbeta, depth, risearch_energy)
if (file.exists(x_binary_file)) {
  print("reusing X file")
  X <- readRDS(x_binary_file)
} else {
  base_dir <- getwd()
  rs2_processing_dir <- sprintf("%s/risearch2_processing_energy%d", base_dir, risearch_energy)
  if (!dir.exists(rs2_processing_dir)) {
      dir.create(rs2_processing_dir)
  }
  print("binary X file does not exist, regenerating")
  if (data == "vaccinia") {
    nz_data <- vaccinia
  } else if (data == "mock") {
    nz_data <- mock
  }
  print("Collecting relevant InfectX data")
  nz_data <- nz_data[!(duplicated(nz_data$Catalog_number)), ] %>%
    filter(!is.na(ID)) %>%
    filter(WellType == "SIRNA") %>%
    filter(!is.na(Catalog_number)) %>%
    filter(!is.na(Sequence_antisense_5_3))
  # nz_data = head(nz_data, 10000)
  vc_sirna_names <- unlist(lapply(nz_data$Catalog_number, as.character))
  ids <- as.character(nz_data$ID)
  ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

  print("Finding interactions with RIsearch2")

  setwd(rs2_processing_dir)
  mrna_suf_file <- "mrna.suf"
  risearch2_bin <- "/home/kieran/work/infx_lasso_data/risearch2/RIsearch-2.1/bin/risearch2.x"
  if (!file.exists(mrna_suf_file)) {
    system(sprintf("%s -c %s/GRCh38_latest_rna.fna -o %s", risearch2_bin, base_dir, mrna_suf_file))
  }

  with_progress(
    {
      # p <- progressor(along=nz_data$Catalog_number)
      suppressed_mrnas <- future_apply(nz_data, 1, future.seed = TRUE, function(well) {
        setwd(rs2_processing_dir)
        sequence <- as.character(well["Sequence_antisense_5_3"])
        name <- as.character(well["Catalog_number"])
        fasta_file <- sprintf("%s.fasta", name)
        output_file <- sprintf("risearch_%s.out.gz", name)
        # p(sprintf("x=%s", name))

        write(sprintf("> %s", name), file = fasta_file)
        write(sequence, file = fasta_file, append = TRUE)

        if (!file.exists(output_file)) {
          # print(sprintf("/home/kieran/work/infx_lasso_data/risearch2/RIsearch-2.1/bin/risearch2.x --index=mrna.suf --query=%s --seed=2:8 --energy=-20", fasta_file))
          system(sprintf("%s --index=mrna.suf --query=%s --seed=2:8 --energy=%d", risearch2_bin, fasta_file, risearch_energy))
        }

        # print(output_file)
        # empty .gz file is 20 bytes
        if (file.info(output_file)$size > 20) {
          print(sprintf("reading %s\n", output_file))
          dt <- fread(output_file)
          dt$mrna <- str_extract(dt$V4, "[^.]*")
          dt$mrna
        } else {
          NULL
        }
      })
    },
    enable = TRUE
  )
  setwd(base_dir)

  print("Finding gene ids for interacting mrna")
  all_mrnas <- unlist(suppressed_mrnas) %>% unique()
  mrna_genes <- getBM(attributes = c("refseq_mrna", "entrezgene_id"), filters = "refseq_mrna", values = all_mrnas, mart = ensembl, uniqueRows = FALSE)
  mrna_to_gene_map <- as.vector(mrna_genes$entrezgene_id)
  names(mrna_to_gene_map) <- mrna_genes$refseq_mrna


  all_genes_affected <- mrna_to_gene_map[!is.na(mrna_to_gene_map)] %>% unique()
  mat <- lapply(suppressed_mrnas, function(set) {
    all_genes_affected %in% mrna_to_gene_map[set]
  })
  X <- t(sapply(mat, unlist))
  rownames(X) <- nz_data$Catalog_number
  colnames(X) <- all_genes_affected
  storage.mode(X) <- "numeric"


  ## not required for RIsearch2
  #X[cbind(vc_sirna_names, ids)] <- 1 # add on-target effects. These are often outside the 3' UTR, so TargetScan doesn't find them.
  # present_mat = X[present_vc_sirnas,]
  # X = present_mat
  # X[X > cutoff] <- 1
  # X[X <= cutoff] <- 0
  saveRDS(X, file = x_binary_file)
}

if (file.exists(yfile) && file.exists(x_filtered_file)) {
  print("reusing Y file & filtered X file")
  Y <- readRDS(yfile)
  X <- readRDS(x_filtered_file)
} else {
  print("Y file (or filtered X) does not exist, regenerating")
  #
  # generate Y
  if (data == "mock") {
    ds <- mock
  } else if (data == "adeno") {
    ds <- adeno
  } else if (data == "vaccinia") {
    # vaccinia$controls <- ifelse(vaccinia$WellType == "CONTROL","yes","no")
    ds <- clean(pathogen = "vaccinia", killers = FALSE, controls = FALSE, quality = "bad")[["vaccinia"]] %>%
      filter(!is.na(eCount_oCells))
  }
  sirnas_in_matrix <- rownames(X)
  names_in_both <- (ds %>% filter(Catalog_number %in% sirnas_in_matrix))$Catalog_number
  X_filtered <- X[as.character(names_in_both), ]
  saveRDS(X_filtered, file = x_filtered_file)
  X <- X_filtered
  # ds = ds %>% filter(!is.na(ID))
  # control_data = ds %>% filter(WellType == "CONTROL")
  # present_sirnas = unlist(dimnames(X)[[1]])
  # ds = ds %>% filter(Catalog_number %in% present_sirnas)
  # ds = ds %>% filter(WellType == "SIRNA")
  ## mock_data = mock %>% filter(WellType != "CONTROL") %>% filter(!is.na(ID))
  ## adeno_data = adeno %>% filter(WellType != "CONTROL") %>% filter(!is.na(ID))

  if (data == "mock") {
    Y <- (ds %>% filter(Catalog_number %in% sirnas_in_matrix))$eCount_oCells_nZScore
  } else {
    Y <- (ds %>% filter(Catalog_number %in% sirnas_in_matrix))$eCount_oCells_nBScore_nZ
  }


  # no_sirna_effect = mean(control_data$eCount_oCells)
  # Y = vector()

  # for (id in sirnas_in_matrix) {
  #    #entries = mock_data %>% filter(Catalog_number == id)
  #    entries = ds %>% filter(Catalog_number == id)
  #    #mean = log2(mean(entries$eCount_oCells))
  #    #mean = log2(mean(entries$eCount_oCells)) - log2(600)
  #    mean = log2(mean(entries$eCount_oCells)) - log2(no_sirna_effect)
  #    Y = append(Y, mean)
  # }
  Y <- as.numeric(Y)
  Y[Y < -100] <- -100
  saveRDS(Y, file = yfile)
}
if (!file.exists(intermediate_result_filename)) {
  num_kos = apply(X, 2, sum)
  X <- X[,num_kos >= ko_min]
  print("no results found for current settings, re-running lasso")
  print("running lasso dom X w/ dimensions")
  print(dim(X))
  gc()
  time <- system.time(
    result <- interaction_lasso(X, Y, max_nz_beta = nzbeta, depth = depth, use_intercept = FALSE, estimate_unbiased = TRUE)
  )

  saveRDS(result, file = intermediate_result_filename)
}


gc()
print("post-processing")
# !/usr/bin/env Rscript

verbose <- TRUE



# ran from 10:40 am to 06:05 am the following day.
lasso_result <- readRDS(intermediate_result_filename)
fit <- lasso_result$estimate_unbiased

## Collect coefficients
fx_main <- fit$main_effects %>% mutate(j = NA, k = NA)
fx_int <- fit$pairwise_effects %>% mutate(k = NA)
fx_trip <- fit$triple_effects

print(sprintf("fx_main: %d ", nrow(fx_main)))
print(sprintf("fx_int: %d ", nrow(fx_int)))
print(sprintf("fx_trip: %d ", nrow(fx_trip)))

## Fit main effects only for comparison
# fit_main_only = lm(Y ~ X)

## Statistical test if b_i and b_ij are sig. > 0
#Z <- cbind(X[, fx_main[["i"]]])
#if (nrow(fx_int) > 0) {
#  for (i in 1:nrow(fx_int)) {
#    Z <- cbind(Z, X[, fx_int[i, ][["i"]], drop = FALSE] * X[, fx_int[i, ][["j"]], drop = FALSE])
#  }
#}
#Z <- as.matrix(Z)
#colnames(Z) <- rownames(Z) <- NULL
#Ynum <- as.numeric(Y)
#fit_red <- lm(Ynum ~ Z)

Z <- X[,fx_main$i]
if (nrow(fx_int) > 0) {
   Z <- cbind(Z, X[,fx_int$i] * X[,fx_int$j])
}
if (nrow(fx_trip) > 0) {
   Z <- cbind(Z, X[,fx_trip$i] * X[,fx_trip$j] * X[,fx_trip$k])
}

print("Z dim:")
print(dim(Z))

ols_time <- system.time(fit_red <- lm(Y ~ Z))

pvals <- data.frame(id = 1:ncol(Z), coef = coef(fit_red)[-1]) %>%
  filter(!is.na(coef)) %>%
  data.frame(., pval = summary(fit_red)$coef[-1, 4]) %>%
  tbl_df()

names <- dimnames(X)[[2]]

smry <- left_join(rbind(fx_main[], fx_int[], fx_trip[]) %>% data.frame(id = 1:nrow(.), .), pvals, by = "id") %>%
  mutate(pval = ifelse(is.na(pval), 1, pval)) %>%
  mutate(i_name = names[i]) %>%
  mutate(j_name = ifelse(is.na(j), NA, names[j])) %>%
  mutate(k_name = ifelse(is.na(k), NA, names[k])) %>%
  # mutate(j_name = names[j]) %>%
  rename(coef.est = coef)

sig <- smry %>% filter(pval < 0.05)
sig_main <- sig[is.na(sig$j), ]
sig_int <- sig[!is.na(sig$j), ]
Z_sig <- cbind(X[, sig_main[["i"]]])
if (nrow(sig_int) > 0) {
  for (i in 1:nrow(sig_int)) {
    Z_sig <- cbind(Z_sig, X[, sig_int[i, ][["i"]], drop = FALSE] * X[, sig_int[i, ][["j"]], drop = FALSE])
  }
}
Z_sig <- as.matrix(Z_sig)
colnames(Z_sig) <- rownames(Z_sig) <- NULL
Ynum <- as.numeric(Y)
fit_red_sig <- lm(Ynum ~ Z_sig)

## Write out
if (verbose) cat("Saving\n")
saveRDS(list(
  fit = fit,
  fx_int = fx_int,
  fx_main = fx_main,
  fit_red = fit_red,
  fit_red_sig = fit_red_sig,
  time = time,
  smry = smry
),
file = result_filename
)
