#!/usr/bin/Rscript

library(data.table)
library(Pint)
require(Matrix)
require(dplyr)
verbose <- TRUE

X_file = "dedup_X.rds"
Y_file = "antibiotic_Y.rds"
max_nz = 50
column_limit <- 0 # the input already has <5 removed
dedup2_file <- sprintf("dedup_X_colLim%d.rds", column_limit)
depth <- 2

print("checking for existing processed data")
if (file.exists(X_file) && file.exists(Y_file)) {
	print("using existing data")
	X = readRDS(X_file)
	Y = readRDS(Y_file)
} else {
	print("data not found, reading table")
	ad = data.table::fread("binary_table_T_anonymised.csv")

	print("formatting data")
	Y = as.numeric(unlist(ad[,"MIC Log2"]))
	tX = ad[,!c("V1", "MIC Log2", "MIC (Cipro)", "MIC_Class")]
	X = as.matrix(tX)
	print("original dimensions")
	dim(X)
	X = t(unique(t(X))) # remove duplicate columns
	print("reduced dimensions")
	dim(X)
	saveRDS(X, file=X_file)
	saveRDS(Y, file=Y_file)
}

#dedup = readRDS("dedup_X.rds")
dedup = X
print(dim(X))

print("removing small columns")
where = apply(dedup, 2, sum) < column_limit
where2 = apply(dedup, 2, sum) >= nrow(X) - column_limit
where_all = (where | where2)
dedup2 = dedup[,!where_all]
#print("NOT removing small columns")
#dedup2 = dedup
print("final dimensions")
dim(dedup2)

saveRDS(dedup2, file=dedup2_file)

results_file <- sprintf("bacteria_lasso_output_2_depth%d_colLim%d.rds", depth, column_limit)

if (!file.exists(results_file)) {
    print("beginning lasso")
    time <- system.time(result <- interaction_lasso(dedup2, Y, depth=depth, max_nz_beta=max_nz))

    print("saving results")
    saveRDS(result, file=results_file)
}
result <- readRDS(results_file)

fit <- readRDS(results_file)
X = readRDS(dedup2_file)
Y = readRDS("./antibiotic_Y.rds")

if (verbose) cat("Collecting stats\n")

## Collect coefficients
fx_main <- fit$main_effects %>% mutate(j=NA)
fx_int <- fit$pairwise_effects

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
			 time = time,
             smry = smry),
        file = sprintf("ifx_full_stats_2.rds"))

write.table(smry, row.names=F, sep=',',
			file=sprintf("ifx_full_stats_2.csv"))
