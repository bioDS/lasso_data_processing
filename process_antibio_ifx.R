#!/usr/bin/env Rscript
require(Matrix)
require(dplyr)


verbose <- TRUE

adcal = TRUE
max_dist = 10


# ran from 10:40 am to 06:05 am the following day.
fit = readRDS(sprintf("./bacteria_lasso_output_adcal%s_maxdist%d.rds", adcal, max_dist))
X = readRDS("./dedup2_X.rds")
Y = readRDS("./antibiotic_Y.rds")

if (verbose) cat("Collecting stats\n")

## Collect coefficients
fx_main <- fit$main_effects %>% mutate(j=NA)
fx_int <- fit$interaction_effects

## Fit main effects only for comparison
#fit_main_only = lm(Y ~ X)

## Statistical test if b_i and b_ij are sig. > 0
Z <- cbind(X[,fx_main[["i"]]])
for (i in 1:nrow(fx_int)) {
  Z <- cbind(Z, X[,fx_int[i,][["i"]], drop = FALSE] * X[,fx_int[i,][["j"]], drop = FALSE])
}
Z <- as.matrix(Z)
colnames(Z) <- rownames(Z) <- NULL
Ynum <- as.numeric(Y)
fit_red <- lm(Ynum ~ Z)


pvals <- data.frame(id = 1:ncol(Z), coef = coef(fit_red)[-1]) %>%
  filter(!is.na(coef)) %>%
  #mutate(coef = ifelse(is.na(coef), 0, coef)) %>%
  data.frame(., pval = summary(fit_red)$coef[-1,4]) %>%
  tbl_df

names = dimnames(X)[[2]]

#keep = !is.na(fit_red$coefficients[-1])

#smry <- left_join(rbind(fx_main[keep[0:length(fx_main$i)],], fx_int[keep[-(0:length(fx_main$i))],]) %>% data.frame(id = 1:nrow(.), .), pvals, by = "id") %>%
smry <- left_join(rbind(fx_main[], fx_int[]) %>% data.frame(id = 1:nrow(.), .), pvals, by = "id") %>%
  mutate(pval = ifelse(is.na(pval), 1, pval)) %>%
  #left_join(., id, by = c("gene_i", "gene_j")) %>%
  mutate(i_name = names[i]) %>%
  mutate(j_name = names[j]) %>%
  rename(coef.est = coef)

## Write out
if (verbose) cat("Saving\n")
saveRDS(list(fit = fit,
             fx_int = fx_int,
             fx_main = fx_main,
             fit_red = fit_red,
             smry = smry),
        file = sprintf("ifx_full_stats_adcal%s_maxdist%d.rds", adcal, max_dist))

write.table(smry, row.names=F, sep=',',
			file=sprintf("ifx_full_stats_adcal%s_maxdist%d.csv", adcal, max_dist))
