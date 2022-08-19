#!/usr/bin/Rscript
require(Matrix)
require(dplyr)

verbose <- TRUE
lambda_min_ratio = 5e-2
if (verbose) cat("Reading data set\n")
Q <- readRDS("./Q1_binary.rds") # possibly not large enough for p > 1000


args <- commandArgs(trailingOnly = TRUE)
#args <- c("1000", "100", "10", "50", "20", "0")

print(args)
n <- args[1] %>% as.numeric
p <- args[2] %>% as.numeric
SNR <- args[3] %>% as.numeric
num_bi <- args[4] %>% as.numeric
num_bij <- args[5] %>% as.numeric
num_bijk <- args[6] %>% as.numeric
num_lethals <- args[7] %>% as.numeric

#' Perturbation matrix
#'
#' Sample a n \times p perturbation matrix from a frequency vector
#' using the Bernoulli distribution
#'
#' @author Fabian Schmich
#' @param n Number of observations
#' @param p Number of genes
#' @return Perturbation matrix
perturbation_matrix <- function(n, p) {
  colnames(Q) <- rownames(Q) <- NULL
  X <- Q[sample(1:nrow(Q), n), sample(1:ncol(Q), p)] %>% as.matrix
  dup <- which(duplicated(t(X)) == TRUE)
  while(length(dup) > 0) {
    repl <- Q[sample(1:nrow(Q), n), sample(1:ncol(Q), length(dup)), drop = FALSE] %>% as.matrix
    for (i in 1:ncol(repl)) {
      X[,dup[i]] <- repl[,i]
    }
    dup <- which(duplicated(t(X)) == TRUE)
  }
  # X <- matrix(0, nrow = n, ncol = p)
  # for (i in 1:nrow(X)) {
  #   X[i,] <- rbinom(n = p, size = 1, prob = sample(x = freqs, size = 1))
  # }
  return(X)
}

## Simulate
# Perturbation matrix
if (verbose) cat("Make perturbation matrix\n")
X <- perturbation_matrix(n = n, p = p)

if (verbose) cat("Finding interactions present\n")

single_obs <- lapply(1:p, function(i) {
  data.frame(o1 = sum(X[,i])) %>%
  mutate(gene_i = i) %>% tibble::as_tibble()
}) %>% do.call("rbind", .) %>% arrange(gene_i) %>% select(gene_i, o1) %>% group_by(gene_i, o1) %>% ungroup


# obs gives us every possible interaction we can use to generate Y.
pairwise_obs <- lapply(2:p, function(i) {
  apply(X[,1:(i-1), drop = FALSE], 2, function(col) {
    data.frame(o00 = (!col & !X[,i]) %>% sum,
               o10 = (col & !X[,i]) %>% sum,
               o01 = (!col & X[,i]) %>% sum,
               o11 = (col & X[,i]) %>% sum)
    }) %>% do.call("rbind", .) %>%
    tibble::as_tibble() %>%
    mutate(gene_j = i, gene_i = 1:(i-1))
}) %>% do.call("rbind", .) %>%
  tibble::as_tibble() %>%
  arrange(gene_i, gene_j) %>%
  select(gene_i, gene_j, o00, o01, o10, o11) %>%
  group_by(gene_i, gene_j, o00, o01, o10, o11) %>% 
  summarise(omin = min(o00, o10, o01, o11)) %>%
  ungroup


threeway_obs <- apply((combn(1:p, 3)), 2, function(a) {
  data.frame(
    gene_i=a[1],
    gene_j=a[2],
    gene_k=a[3],
    o111=sum(X[,a[1]] * X[,a[2]] * X[,a[3]])
  )
}) %>% do.call("rbind", .) %>% tibble::as_tibble()

if (verbose) cat("Choosing effect coefficients\n")
# main effects
nz_single_obs = single_obs %>% filter(o1 > 0)
bi_ind <- nz_single_obs[sample(nrow(nz_single_obs), num_bi),] %>%
  rowwise %>% mutate(coef = rnorm(1, mean = 0, sd = 2)) # double amplitude sd

# pairwise
nz_pairwise_obs = pairwise_obs %>% filter(o11 > 0)
bij_ind <- nz_pairwise_obs[sample(nrow(nz_pairwise_obs), num_bij),] %>%
  rowwise %>% mutate(coef = rnorm(1, mean = 0, sd = 2)) # double amplitude sd

# pairwise lethal
lethal_ind <- nz_pairwise_obs[sample(nrow(nz_pairwise_obs), num_lethals),] %>%
  mutate(coef = -1000) # double amplitude sd

# three-way
nz_threeway_obs = threeway_obs %>% filter(o111 > 0)
bijk_ind <- nz_threeway_obs[sample(nrow(nz_threeway_obs), num_bijk),] %>%
  rowwise %>% mutate(coef = rnorm(1, mean = 0, sd = 2)) # double amplitude sd

## Fitness
if (verbose) cat("Sampling fitness\n")
Y <- X[,bi_ind[["gene_i"]], drop = FALSE] %*% bi_ind[["coef"]]
for (i in 1:nrow(bij_ind)) {
  Y <- Y + (X[,bij_ind[i,][["gene_i"]], drop = FALSE] * X[,bij_ind[i,][["gene_j"]], drop = FALSE]) %*% bij_ind[i,][["coef"]]
}
for (i in 1:nrow(bijk_ind)) {
  Y <- Y + (X[,bijk_ind[i,][["gene_i"]], drop = FALSE] * X[,bijk_ind[i,][["gene_j"]], drop=FALSE] * X[,bijk_ind[i,][["gene_k"]], drop = FALSE] %*% bijk_ind[i,][["coef"]])
}
if (num_lethals > 0) {
    for (i in 1:nrow(lethal_ind)) {
      Y <- Y + (X[,lethal_ind[i,][["gene_i"]], drop = FALSE] * X[,lethal_ind[i,][["gene_j"]], drop = FALSE]) %*% lethal_ind[i,][["coef"]]
    }
}
## add noise
noise <- (rnorm(n = nrow(Y), mean = 0, sd = 1))
Y <- Y + sqrt(var(Y[,1])/(SNR * var(noise))) * noise

## Write out
if (verbose) cat("Saving\n")
saveRDS(list(X=X, Y=Y,
             bi_ind = bi_ind,
             bij_ind = bij_ind,
             bijk_ind = bijk_ind,
             lethal_ind = lethal_ind,
             single_obs=single_obs,
             pairwise_obs=pairwise_obs,
             threeway_obs=threeway_obs),
        file = sprintf("./simulated_data/n%d_p%d_SNR%d_nbi%d_nbij%d_nbijk%d_nlethals%d_%d.rds",
                       n, p, SNR, num_bi, num_bij, num_bijk, num_lethals, (runif(1) * 1e5) %>% floor))
