#!/usr/bin/env Rscript
library(Pint)
library(dplyr)
library(ggplot2)
library(digest)

# f = readRDS("simulated_data/n1000_p100_SNR5_nbi20_nbij30_nbijk5_nlethals0_94723.rds")
# f = readRDS("simulated_data/n1000_p100_SNR5_nbi20_nbij30_nbijk50_nlethals0_12086.rds")
# f = readRDS("simulated_data/n1000_p100_SNR10_nbi10_nbij100_nbijk1000_nlethals0_85437.rds")

# Single value for now
find_matching_tp <- function(vals, X, true_hashes) {
    cols = X[,vals]
    if (length(vals) > 1) {
        actual_col = apply(X = cols, FUN = prod, MARGIN = 1)
    } else {
        actual_col = cols
    }
    hash = digest(actual_col)
    return(hash %in% true_hashes)
}

get_true_effect_hashes <- function(f) {
    X = f$X
    # main effects
    main_hashes = apply(X[,f$bi_ind$gene_i], 2, digest)
    #pairwise
    pairwise_inds = rbind(f$bij_ind$gene_i, f$bij_ind$gene_j)
    pairwise_hashes = apply(pairwise_inds, 2, function(x) {
        pairwise_cols = X[,x]
        pairwise_actual_col = apply(pairwise_cols, 1, prod)
        return(digest(pairwise_actual_col))
    })
    #threeway
    threeway_inds = rbind(f$bijk_ind$gene_i, f$bijk_ind$gene_j, f$bijk_ind$gene_k)
    threeway_hashes = apply(threeway_inds, 2, function(x) {
        threeway_cols = X[,x]
        threeway_actual_col = apply(threeway_cols, 1, prod)
        return(digest(threeway_actual_col))
    })
    all_hashes = c(main_hashes, pairwise_hashes, threeway_hashes)
}

process <- function(ret, f) {
  all_hashes = get_true_effect_hashes(f)
  me <- ret$main
  fx_main <- data.frame(gene_i = as.numeric(me$effects$i), effect = me$effects$strength) %>%
      arrange(gene_i) %>%
    left_join(., f$single_obs, by = c("gene_i")) %>%
    mutate(type = "main") %>%
    rowwise() %>%
    left_join(., f$bi_ind, by = c("gene_i")) %>%
    mutate(TP = !is.na(coef), gene_j = NA, gene_k = NA, coef.est = effect) %>%
    select(gene_i, gene_j, gene_k, TP, coef, coef.est, type) %>%
    arrange(desc(TP))
  fx_main["equiv_tp"] = sapply(fx_main$gene_i, function(x) {
      return(find_matching_tp(x, f$X, all_hashes))
  })

  pe <- ret$pairwise$effects
  fx_pairwise <- data.frame(gene_i = as.numeric(pe$i),
                            gene_j = as.numeric(pe$j),
                            effect = pe$strength) %>%
    arrange(gene_i, gene_j) %>%
    left_join(., f$pairwise_obs, by = c("gene_i", "gene_j")) %>%
    mutate(type = "pairwise") %>%
    rowwise() %>%
    left_join(., f$bij_ind, by = c("gene_i", "gene_j")) %>%
    mutate(TP = !is.na(coef), gene_k = NA, coef.est = effect) %>%
    select(gene_i, gene_j, gene_k, TP, coef, coef.est, type) %>%
    arrange(desc(TP))
  fx_pairwise["equiv_tp"] = apply(fx_pairwise, 1, function(x) {
      cols = c(as.numeric(x["gene_i"]), as.numeric(x["gene_j"]))
      return(find_matching_tp(cols, f$X, all_hashes))
  })

  tre <- ret$triple$effects
  fx_threeway <- data.frame(gene_i = as.numeric(tre$i), gene_j = as.numeric(tre$j), gene_k = as.numeric(tre$k), effect = tre$strength) %>%
    arrange(gene_i, gene_j, gene_k) %>%
    left_join(., f$threeway_obs, by = c("gene_i", "gene_j", "gene_k")) %>%
    mutate(type = "three-way") %>%
    rowwise() %>%
    left_join(., f$bijk_ind, by = c("gene_i", "gene_j", "gene_k")) %>%
    mutate(TP = !is.na(coef), coef.est = effect) %>%
    select(gene_i, gene_j, gene_k, TP, coef, coef.est, type) %>%
    arrange(desc(TP))
  fx_threeway["equiv_tp"] = apply(fx_threeway, 1, function(x) {
      cols = c(as.numeric(x["gene_i"]), as.numeric(x["gene_j"]), as.numeric(x["gene_k"]))
      return(find_matching_tp(cols, f$X, all_hashes))
  })

  X <- f$X
  # Z <- X[,me$i]
  # if (nrow(pe) > 0) {
  #    Z <- cbind(Z, X[,pe$i] * X[,pe$j])
  # }
  # if (nrow(tre) > 0) {
  #    Z <- cbind(Z, X[,tre$a] * X[,tre$b] * X[,tre$c])
  # }

  # ols_time = system.time(fit_red <- lm(f$Y ~ Z))

  # pvals <- data.frame(id = 1:ncol(Z), coef = coef(fit_red)[-1]) %>%
  #  filter(!is.na(coef)) %>%
  #  data.frame(., pval = summary(fit_red)$coef[-1,4]) %>%
  #  mutate(coef.lm_est=coef) %>%
  #  select(id, coef.lm_est, pval) %>%
  #  tbl_df

  # coefs = summary(fit_red)$coef
  # coefs = fit_red$coefficients
  # tmp = data.frame(rowname=rownames(coefs), coefs)
  # pv_tmp = data.frame(rowname=rownames(pvals), pvals)
  # smry = merge(tmp, pv_tmp, all=TRUE, by="rowname") %>%
  #    arrange(as.numeric(id)) %>% filter(!is.na(id)) %>%
  #    select(id, coef.lm_est, pval)


  # smry <- merge(rbind(fx_main, fx_pairwise, fx_threeway) %>% data.frame(id = 1:nrow(.), .), smry, all=TRUE, by = "id") #%>%
  smry <- data.frame(rbind(fx_main, fx_pairwise, fx_threeway))
  #   mutate(pval = ifelse(is.na(pval), 1, pval)) %>%
  #   rename(coef.est = coef) %>%
  #   left_join(., obs, by = c("gene_i", "gene_j"))

  # Z <- cbind(X[,pint_fx_main[["gene_i"]]])
  # if (nrow(interactions) > 0) {
  #  for (i in 1:nrow(pint_fx_int)) {
  #    Z <- cbind(Z, X[,pint_fx_int[i,][["gene_i"]], drop = FALSE] * X[,pint_fx_int[i,][["gene_j"]], drop = FALSE])
  #  }
  # }
  # Z <- as.matrix(Z)
  # print(summary(fx_main$TP))
  # print(summary(fx_pairwise$TP))
  # print(summary(fx_threeway$TP))
  # return(list(fx_main=fx_main, fx_pairwise=fx_pairwise, fx_threeway=fx_threeway))
  return(smry)
}

run_23way <- function(fname) {
  f <- readRDS(fname)
  # running_time <- system.time(
  ret3 <- interaction_lasso(X = f$X, Y = f$Y, depth = 3, check_duplicates = TRUE, max_nz_beta = 100)
  # )
  ret2 <- interaction_lasso(X = f$X, Y = f$Y, depth = 2, check_duplicates = TRUE, max_nz_beta = 100)

  # saveRDS(list(ret2=ret2, ret3=ret3), "lasso_both.rds")


  depth3_results <- process(ret3, f)
  depth2_results <- process(ret2, f)

  cutoff <- sd(depth3_results$coef.est)
  # depth3_results %>% filter(abs(coef.est) > 0.5)
  # depth3_results %>%
  #   select(TP) %>%
  #   summary()
  # depth3_results %>% filter(abs(coef.est) > 0.5) %>% select(TP) %>% summary
  # depth3_results %>% filter(abs(coef.est) > sd(coef.est)) %>% select(TP) %>% summary
  # depth3_results %>% filter(abs(coef.est) > 2*sd(coef.est)) %>% select(TP) %>% summary
  # depth3_results %>% filter(abs(coef.est) > 4*sd(coef.est)) %>% select(TP) %>% summary

  #depth3_results %>% filter(abs(coef.est) > 2 * sd(coef.est))
  #depth2_results %>% filter(abs(coef.est) > 2 * sd(coef.est))

  #depth2_results %>%
  #  filter(abs(coef.est) > 2 * sd(coef.est)) %>%
  #  select(TP) %>%
  #  summary()

  # rbind(depth2_results, depth3_results) %>% filter(abs(coef.est) > 2*sd(coef.est))

  all_real_effects <- rbind(
    f$bi_ind %>% mutate(gene_j = NA, gene_k = NA) %>% select(gene_i, gene_j, gene_k, coef),
    f$bij_ind %>% mutate(gene_k = NA) %>% select(gene_i, gene_j, gene_k, coef),
    f$bijk_ind %>% select(gene_i, gene_j, gene_k, coef)
  )

  big3_results <- depth3_results %>% filter(abs(coef.est) > cutoff)
  big2_results <- depth2_results %>% filter(abs(coef.est) > cutoff)

  big_exist <- sum(abs(all_real_effects$coef) > cutoff)
  big_real <- all_real_effects[abs(all_real_effects$coef) > cutoff, ]
  big_main <- big_real[is.na(big_real$gene_j), ]
  big_pair <- big_real[is.na(big_real$gene_k) & !is.na(big_real$gene_j), ]
  big_triple <- big_real[!is.na(big_real$gene_k) & !is.na(big_real$gene_j), ]

  all_main <- all_real_effects[is.na(all_real_effects$gene_j), ]
  all_pair <- all_real_effects[is.na(all_real_effects$gene_k) & !is.na(all_real_effects$gene_j), ]
  all_triple <- all_real_effects[!is.na(all_real_effects$gene_k) & !is.na(all_real_effects$gene_j), ]

  # big_found <- nrow(depth3_results %>% filter(abs(coef.est) > cutoff) %>% filter(TP==TRUE))

  # depth3_results %>% ggplot(aes(x=0, y=coef.est)) +  geom_point()
  #ggplot(depth3_results %>% filter(TP == FALSE), aes(x = coef.est)) +
  #  geom_histogram(binwidth = 0.1) +
  #  geom_density() +
  #  # geom_density(data=data.frame(sample=c(rnorm(nrow(depth3_results)))), aes(x=sample), color="red")
  #  geom_histogram(data = (depth3_results %>% filter(TP == TRUE)), binwidth = 0.1, color = "red") +
  #  geom_vline(xintercept = cutoff) +
  #  geom_vline(xintercept = -cutoff)


  get_prec_rec <- function(rs) {
    # pr <- data.frame(type=NA, prec=NA, rec=NA)

    main <- rs %>% filter(type == "main")
    pair <- rs %>% filter(type == "pairwise")
    triple <- rs %>% filter(type == "three-way")

    main_prec <- sum(main$TP) / nrow(main)
    pair_prec <- sum(pair$TP) / nrow(pair)
    triple_prec <- sum(triple$TP) / nrow(triple)

    main_rec <- sum(main$TP) / nrow(all_main)
    pair_rec <- sum(pair$TP) / nrow(all_pair)
    triple_rec <- sum(triple$TP) / nrow(all_triple)

    total_effects = nrow(all_main) + nrow(all_pair) + nrow(all_triple)
    prec_overall <- sum(rs$TP) / nrow(rs)
    rec_overall <- sum(rs$TP) / total_effects

    pr <- data.frame(
      prec_main = main_prec, rec_main = main_rec,
      prec_pair = pair_prec, rec_pair = pair_rec,
      prec_trip = triple_prec, rec_trip = triple_rec,
      prec_overall = prec_overall, rec_overall = rec_overall
    )

    return(pr)
  }

  big3_pr <- get_prec_rec(big3_results)
  big2_pr <- get_prec_rec(big2_results)

  all3_pr <- get_prec_rec(depth3_results)
  all2_pr <- get_prec_rec(depth2_results)

  return(list(depth2_results = depth2_results, depth3_results = depth3_results, big2_results = big2_results, big3_results = big3_results, big2_pr = big2_pr, big3_pr = big3_pr, all2_pr = all2_pr, all3_pr = all3_pr))
}

# saveRDS(run_23way(f), "results.rds")
# saveRDS(list(depth2_results=depth2_results, depth3_results=depth3_results, big2_results=big2_results, big3_results=big3_results, big2_pr=big2_pr, big3_pr=big3_pr), "results.rds")
