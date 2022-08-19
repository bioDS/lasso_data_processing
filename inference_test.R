#!/usr/bin/env Rscript
library(selectiveInference)
#detach("package:selectiveInference", unload=TRUE)
#library(selectiveInference)
library(dplyr)
library(Pint)

f <- readRDS("simulated_data/n100_p20_SNR5_nbi10_nbij20_nbijk30_nlethals0_15743.rds")
# f <- readRDS("../simulated_data/simulated_data_small_repeat/n1000_p100_SNR5_nbi20_nbij20_nlethals0_viol0_9284.rds")
X <- f$X
Y <- f$Y

n = nrow(X)
p = ncol(X)

p_int = p*(p+1)/2
pairwise_col_ind = matrix(nrow = p, ncol = p)
inters_at_val = matrix(nrow=p_int, ncol=2)

X2 <- matrix(nrow = n, ncol = p_int)

count = 1;
for (main in 1:p) {
    for (inter in main:p) {
        X2[,count] <- X[,main] * X[,inter]
        pairwise_col_ind[main, inter] = count
        inters_at_val[count,1] <- main
        inters_at_val[count,2] <- inter
        count <- count+1
    }
}

pr = pairwise_lasso(X, Y, depth=2)

beta = rep(0.0, p_int)

for (ind in 1:length(pr$main_effects$i)) {
    i = pr$main_effects$i[ind]
    col = pairwise_col_ind[i, i]
    beta[col] <- pr$main_effects$strength[ind]
}
for (ind in 1:length(pr$pairwise_effects$i)) {
    i = pr$pairwise_effects$i[ind]
    j = pr$pairwise_effects$j[ind]
    col = pairwise_col_ind[i, j]
    beta[col] <- pr$pairwise_effects$strength[ind]
}

ndX2 = X2[,!duplicated(t(X2))]
ndX2 = X2

glm = cv.glmnet(x=ndX2, Y, familly="guassian", intercept=FALSE, thresh=1e-10)
glm_beta = glm$beta
lambdas = length(glm$lambda)
# ll = glm$lambda.1se*n
ll = glm$lambda.min*n

glm_beta <- coef(glm, x=ndX2, y=Y, s=ll/n, exact=TRUE)[-1]

# fixedLassoInf(x=ndX2, y=Y, lambda=ll, beta=coef(glm, x=ndX2, y=Y, s=ll/n, exact=TRUE)[-1], intercept=FALSE)
# glm_ci = ROSI(X=ndX2, y=Y, lambda=ll, soln=glm_beta)
#detach("package:selectiveInference", unload=TRUE)
#library(selectiveInference)
ci = ROSI(X=X2, y=Y, lambda=0.05, soln=beta, verbose=TRUE, solver="glmnet")
# fixedLassoInf(x=X2, y=Y, lambda=0.051342, beta=beta)

use_beta=beta

nz_glm_beta = glm_beta[glm_beta != 0]

# summary = data.frame(col=1:p_int, gene_i=NA, gene_j=NA, beta=glm_beta, pval=NA, lo=NA, hi=NA, rosi_coef=NA, actual_coef=NA)
summary = data.frame(col=1:p_int, gene_i=NA, gene_j=NA, beta=use_beta, pval=NA, lo=NA, hi=NA, rosi_coef=NA, actual_coef=NA, equivalent=NA)

count <- 1
for (i in 1:length(use_beta)) {
    if (glm_beta[i] != 0) {
        summary$rosi_coef[i] <- ci$estimate[count]
        summary$pval[i] <- ci$pvalues[count]
        summary$lo[i] <- ci$intervals[count,1]
        summary$hi[i] <- ci$intervals[count,2]
        summary$gene_i[i]=inters_at_val[summary$col[i],1]
        summary$gene_j[i]=inters_at_val[summary$col[i],2]
        #if (summary$gj[i] == summary$gi[i] && f$bi_ind$coef) {
        #    summary$actual_coef[i]<-(f$bi_ind %>% filter(gene_i == summary$gi[i]))$coef
        #}

        count <- count+1
    }
}

# summary <- left_join(summary, f$bi_ind, by="gene_i")
smry_main = summary %>% filter(gene_i == gene_j) %>% left_join(f$bi_ind)
smry_int  = summary %>% filter(gene_i != gene_j) %>% left_join(f$bij_ind)

sig_main <- smry_main %>% filter(sign(lo) == sign(hi))
sig_int  <- smry_int  %>% filter(sign(lo) == sign(hi))
# smry_int %>% filter(sign(lo) == sign(hi)) %>% filter(sign(hi) == sign(beta))

for (sigrow in 1:nrow(sig_int)) {
    p_i <- sig_int[sigrow,]$gene_i
    p_j <- sig_int[sigrow,]$gene_j

    comp_x <- X[,p_i]*X[,p_j]

    equal <- data.frame(i=NA, j=NA)
    count=1
    for (main in 1:p) {
        for (int in main:p) {
            if (sum(X[,main]*X[,int] == comp_x) == 100) {
                equal[count,] = c(main, int)
                count <- count+1
            }
        }
    }
    equal

    true_effects <- NULL
    for (ind in 1:nrow(equal)) {
        i = equal[ind,]$i
        j = equal[ind,]$j
        true_effects = rbind(true_effects, f$bij_ind %>% filter(gene_i == i) %>% filter(gene_j == j))
    }

    if (nrow(true_effects) > 0) {
        # print(p_i)
        # print(p_j)
        sig_int$equivalent[sigrow] <- TRUE
    }
}

sig_int

# summary %>% filter(sign(lo) == sign(hi))

