#!/usr/bin/Rscript
require(Matrix)
require(dplyr)
library(foreach)
library(doMC)

registerDoMC(cores=detectCores())

verbose <- TRUE
lambda_min_ratio = 5e-2
if (verbose) cat("Reading data set\n")
Q <- readRDS("./Q1_binary.rds") # possibly not large enough for p > 1000

generate_set <- function(n, p, SNR, num_bi, num_bij, num_bijk, num_lethals) {
    # args <- commandArgs(trailingOnly = TRUE)
    #args <- c("1000", "100", "10", "50", "20", "0")

    # print(args)
    # n <- args[1] %>% as.numeric
    # p <- args[2] %>% as.numeric
    # SNR <- args[3] %>% as.numeric
    # num_bi <- args[4] %>% as.numeric
    # num_bij <- args[5] %>% as.numeric
    # num_bijk <- args[6] %>% as.numeric
    # num_lethals <- args[7] %>% as.numeric

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

    ## single_obs <- lapply(1:p, function(i) {
    ##   data.frame(o1 = sum(X[,i])) %>%
    ##   mutate(gene_i = i) %>% tibble::as_tibble()
    ## }) %>% do.call("rbind", .) %>% arrange(gene_i) %>% select(gene_i, o1) %>% group_by(gene_i, o1) %>% ungroup

    ## single_obs <- data.frame(gene_i = 1:p, o1 = colSums(X))

    # single
    ## tmp <- apply(X, 2, sum)

    # pair
    ## cX <- crossprod(X, X)
    ## cX2 <- einsum("ik,jk->ik", X, X)
    ## tmp2 <- lapply(1:p, function(i) {
    ##     lapply(i:p, function(j) {
    ##         return(cX[i,j])
    ##     })
    ## })
    ## tmp2 <- tmp2 %>% do.call("rbind", .) %>% tibble::as_tibble()

    ## pairwise_obst <- apply((combn(1:p, 2)), 2, function(a) {
    ##   data.frame(
    ##     gene_i=a[1],
    ##     gene_j=a[2],
    ##     o11=cX[a[1], a[2]]
    ##   )
    ## }) %>% do.call("rbind", .) %>% tibble::as_tibble()

    # torch is /much/ faster than R for this sort of thing
    # single
    library(torch)
    X_tensor <- torch_tensor(X)
    test1 <- torch_einsum("ij->j", list(X_tensor))
    single_obs = data.frame(gene_i = 1:p, o1 = as_array(test1))

    # pairwise
    test2 <- torch_einsum("xi,xj->ij", list(X_tensor, X_tensor))
    t2_arr = as_array(test2)
    ## all_ij <- combn(1:p, 2)
    all_ij <- as_array(torch_combinations(torch_tensor(1:p), 2))
    pairwise_obs = data.frame(gene_i = all_ij[,1], gene_j = all_ij[,2], o11 = t2_arr[all_ij])

    # triple (subsample)
    if (num_bijk > 0) {
        subsample_size <- min(num_bijk, p)
        ijk_subsample <- sample(1:p, subsample_size)
        small_X_tensor <- torch_tensor(X[,ijk_subsample])
        test3 <- torch_einsum("xi,xj,xk->ijk", list(small_X_tensor, small_X_tensor, small_X_tensor))
        test3 <- as_array(test3)
        available_ijk <- combn(ijk_subsample, 3)
        ijk_sub_tensor <- torch_tensor(ijk_subsample)
        ijk_sub_indices <- as_array(torch_combinations(torch_tensor(1:subsample_size), 3))
        available_ijk <- torch_combinations(ijk_sub_tensor, r=3, with_replacement = FALSE)
        available_ijk <- as_array(available_ijk)
        threeway_obs <- data.frame(gene_i = available_ijk[,1], gene_j = available_ijk[,2],  gene_k = available_ijk[,3], o111 = test3[ijk_sub_indices])
    } else {
        threeway_obs <- NA
        bijk_ind <- NA
    }

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
    if (num_bijk > 0) {
        nz_threeway_obs = threeway_obs %>% filter(o111 > 0)
        cat(sprintf("nz_3way: %d, sampling: %d\n", nrow(nz_threeway_obs), num_bijk))
        bijk_ind <- nz_threeway_obs[sample(nrow(nz_threeway_obs), num_bijk),] %>%
        rowwise %>% mutate(coef = rnorm(1, mean = 0, sd = 2)) # double amplitude sd
    }

    ## Fitness
    if (verbose) cat("Sampling fitness\n")
    Y <- X[,bi_ind[["gene_i"]], drop = FALSE] %*% bi_ind[["coef"]]
    for (i in 1:nrow(bij_ind)) {
      Y <- Y + (X[,bij_ind[i,][["gene_i"]], drop = FALSE] * X[,bij_ind[i,][["gene_j"]], drop = FALSE]) %*% bij_ind[i,][["coef"]]
    }
    if (num_bijk > 0) {
        for (i in 1:nrow(bijk_ind)) {
        Y <- Y + (X[,bijk_ind[i,][["gene_i"]], drop = FALSE] * X[,bijk_ind[i,][["gene_j"]], drop=FALSE] * X[,bijk_ind[i,][["gene_k"]], drop = FALSE] %*% bijk_ind[i,][["coef"]])
        }
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
    # if (verbose) cat("Saving\n")
    return(list(X=X, Y=Y,
                 bi_ind = bi_ind,
                 bij_ind = bij_ind,
                 bijk_ind = bijk_ind,
                 lethal_ind = lethal_ind,
                 single_obs=single_obs,
                 pairwise_obs=pairwise_obs,
                 threeway_obs=threeway_obs))
}
