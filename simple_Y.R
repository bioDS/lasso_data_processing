#!/usr/bin/Rscript

library(LassoTesting)
library(dplyr)
library(Matrix)

yfile = "infx_simple_simulated_y_full_risearch.rds"
result_filename = "infx_simulated_full_int_smry_risearch.rds"

main_effects_file = "sim_y_full_main_effects.rds"
int_effects_file = "sim_y_full_int_effects.rds"

#Y = readRDS("infx_Y_vaccinia.rds")
#X = readRDS("infx_X_binary_cutoff0_vaccinia.rds")
#X = readRDS("infx_X_risearch_10k_cutoff.rds")
X = readRDS("infx_X_binary_risearch_vaccinia.rds")

n = nrow(X)
p = ncol(X)

if (file.exists(yfile)) {
    Y = readRDS(yfile)
} else {
    if (file.exists(main_effects_file)) {
        main_effects = readRDS(main_effects_file)
    } else {
        print("simulating main effects")
        main_effects = vector(mode="numeric", p)
        for (col in 1:p) {
            if (runif(1) * 10 < 1) {
                main_effects[col] = rnorm(1, mean = 0, sd = 2)
            }
        }
        saveRDS(main_effects, main_effects_file)
    }

    if (file.exists(int_effects_file)) {
        int_effects = readRDS(int_effects_file)
    } else {
        print("simulating interactions")
        int_effects = Matrix(0,p,p, sparse=TRUE)
        for (col in 1:p) {
            for (row in (col+1):p) {
                if (runif(1) * 1000 < 1) {
                    int_effects[row,col] = rnorm(1, mean = 0, sd = 2)
                }
            }
        }
        saveRDS(int_effects, int_effects_file)
    }

    if (file.exists(yfile)) {
        Y = readRDS(yfile)
    } else {
        print("adding noise and finding Y")
        Y = vector(mod="numeric", n)
        # main effects
        for (row in 1:n) {
            Y[row] = sum(unlist(main_effects[X[row,] == 1])) # main effects
            + sum((int_effects[X[row,] == 1, X[row,] == 1])) # interactions
            + rnorm(1, mean = 0, sd = 10) # noise
        }
        saveRDS(Y, yfile)
    }
}


fit = overlap_lasso(as.matrix(X), as.numeric(Y), max_interaction_distance=1000)

## Collect coefficients
fx_main <- fit$main_effects %>% mutate(j=NA)
fx_int <- fit$interaction_effects

## Fit main effects only for comparison
#fit_main_only = lm(Y ~ X)

## Statistical test if b_i and b_ij are sig. > 0
Z <- cbind(X[,fx_main[["i"]]])
if (nrow(fx_int) > 0) {
    for (i in 1:nrow(fx_int)) {
      Z <- cbind(Z, X[,fx_int[i,][["i"]], drop = FALSE] * X[,fx_int[i,][["j"]], drop = FALSE])
    }
}
Z <- as.matrix(Z)
colnames(Z) <- rownames(Z) <- NULL
Ynum <- as.numeric(Y)
fit_red <- lm(Ynum ~ Z)


pvals <- data.frame(id = 1:ncol(Z), coef = coef(fit_red)[-1]) %>%
  filter(!is.na(coef)) %>%
  data.frame(., pval = summary(fit_red)$coef[-1,4]) %>%
  tbl_df

names = dimnames(X)[[2]]

smry <- left_join(rbind(fx_main[], fx_int[]) %>% data.frame(id = 1:nrow(.), .), pvals, by = "id") %>%
  mutate(pval = ifelse(is.na(pval), 1, pval)) %>%
  mutate(i_name = names[i]) %>%
  mutate(j_name = ifelse(is.na(j), NA, names[j])) %>%
  #mutate(j_name = names[j]) %>%
  rename(coef.est = coef)

sig = smry %>% filter(pval < 0.05)
sig_main = sig[is.na(sig$j),]
sig_int = sig[!is.na(sig$j),]
Z_sig <- cbind(X[,sig_main[["i"]]])
if (nrow(sig_int) > 0) {
    for (i in 1:nrow(sig_int)) {
      Z_sig <- cbind(Z_sig, X[,sig_int[i,][["i"]], drop = FALSE] * X[,sig_int[i,][["j"]], drop = FALSE])
    }
}
Z_sig <- as.matrix(Z_sig)
colnames(Z_sig) <- rownames(Z_sig) <- NULL
Ynum <- as.numeric(Y)
fit_red_sig <- lm(Ynum ~ Z_sig)

## Write out
saveRDS(list(fit = fit,
             fx_int = fx_int,
             fx_main = fx_main,
             fit_red = fit_red,
             fit_red_sig = fit_red_sig,
             smry = smry),
        file = result_filename)
