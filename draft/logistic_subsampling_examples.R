source('R/family_expand.R')
source('R/glm_subsampling_functions.R')
source('draft/generate_logit_data.R')
rm(list = setdiff(ls(), lsf.str()))
# library(RcppArmadillo)
#
# library(Rcpp)
# sourceCpp("R/do_it.cpp")
N <-  1e4
beta0 <- rep(0.5, 7)
# beta0 <- c(1, rep(0.5, 6))
d <- length(beta0)
seed <- 1
set.seed(seed)
#### JASA2018 simulation########
family <- binomial.expand()
# FULL.DATA <- simu_mzNormal(seed = seed, N = N, beta0 = beta0, corr = 0.5)
FULL.DATA <- simu_mzNormal(seed = seed, N = N, beta0 = beta0, corr = 0.5)
# FULL.DATA <- simu_ueNormal(seed = seed, N = N, beta0 = beta0, corr = 0.5)
# FULL.DATA <- simu_nzNormal(seed = seed, N = N, beta0 = beta0, corr = 0.5)
# FULL.DATA <- simu_EXP(seed = seed, N = N, beta0 = beta0, rate = 2) # simu 6, EXP
# FULL.DATA <- simu_ueNormal(seed = j, N = N, beta0 = beta0, corr = 0.5) # simu 3, multinormal, Hetero
# FULL.DATA <- simu_T3(seed = seed, N = N, beta0 = beta0, corr = 0.5, df = 3)  # simu 5, T3
X.full <- FULL.DATA$X
Y.full <- FULL.DATA$Y
print(paste('sum: ', sum(Y.full)))
data <- as.data.frame(cbind(Y.full, X.full))
formula <- Y.full ~ .  # if your design matrix already has an intercept column, you should add '-1' to the formula to avoid adding another intercept column by 'model.matrix()'.
n.plt <- 200
n.ssp <- 500
#
# FULL <- logistic.coef.estimate(X = cbind(1, X.full), Y = Y.full, family=family)$beta
# MSE_FULL <- sum((FULL - beta0)^2)
# Uni <-  logistic.optimal.subsampling(formula, data, n.plt, n.ssp,
#                                               estimate.method = 'Uni')

# subsampling.summary(Uni)
###############################################################################
#
Wet_WR_optA <-  glm.optimal.subsampling(formula, data, n.plt, n.ssp,
                                             family = family$family.name,
                                             criterion = 'OptA',
                                             sampling.method = 'Poisson', # WithReplacement, Poisson
                                             estimate.method = 'Weighted',
                                             alpha = 0,
                                             b=2)


Wet_WR_optA_plt <- sum((Wet_WR_optA$beta.plt - beta0)^2)
Wet_WR_optA_ssp <- sum((Wet_WR_optA$beta.ssp - beta0)^2)
Wet_WR_optA_cmb <- sum((Wet_WR_optA$beta.cmb - beta0)^2)
Wet_WR_optA_var_plt <- sum(diag(Wet_WR_optA$var.plt))
Wet_WR_optA_var_ssp <- sum(diag(Wet_WR_optA$var.ssp))
Wet_WR_optA_var_ssp_true <- sum(diag(Wet_WR_optA$var.ssp.true))
Wet_WR_optA_var_cmb <- sum(diag(Wet_WR_optA$var.cmb))
Wet_WR_optA_var_cmb_true <- sum(diag(Wet_WR_optA$var.cmb.true))
###############################################################################
Wet_Poi_optA <-  glm.optimal.subsampling(formula, data, n.plt, n.ssp,
                                              family = family$family.name,
                                              criterion = 'LCC',
                                              sampling.method = 'Poisson', # WithReplacement, Poisson
                                              estimate.method = 'Weighted',
                                              alpha = 0,
                                              b=2)
subsampling.summary(Wet_Poi_optA)
Wet_Poi_optA_ssp <- sum((Wet_Poi_optA$beta.ssp - beta0)^2)
Wet_Poi_optA_cmb <- sum((Wet_Poi_optA$beta.cmb - beta0)^2)
Wet_Poi_optA_var_plt <- sum(diag(Wet_Poi_optA$var.plt))
Wet_Poi_optA_var_ssp <- sum(diag(Wet_Poi_optA$var.ssp))
Wet_Poi_optA_var_ssp_true <- sum(diag(Wet_Poi_optA$var.ssp.true))
Wet_Poi_optA_var_cmb <- sum(diag(Wet_Poi_optA$var.cmb))
Wet_Poi_optA_var_cmb_true <- sum(diag(Wet_Poi_optA$var.cmb.true))

#
# Odds
Odds_Poi_optA <-  logistic.optimal.subsampling(formula, data, n.plt, n.ssp,
                                               family = family$family.name,
                                                   criterion = 'LCC',
                                                   sampling.method = 'Poisson', # WithReplacement, Poisson
                                                   estimate.method = 'LogOddsCorrection',
                                                   alpha = 0,
                                                   b=2)
Odds_Poi_optA_ssp <- sum((Odds_Poi_optA$beta.ssp - beta0)^2)
Odds_Poi_optA_cmb <- sum((Odds_Poi_optA$beta.cmb - beta0)^2)
#



#
Odds_WR_optA <-  logistic.optimal.subsampling(formula, data, n.plt, n.ssp,
                                              family = family$family.name,
                                               criterion = 'optA',
                                               sampling.method = 'WithReplacement', # WithReplacement, Poisson
                                               estimate.method = 'LogOddsCorrection',
                                               b=2)
Odds_WR_optA_ssp <- sum((Odds_WR_optA$beta.ssp - beta0)^2)
Odds_WR_optA_cmb <- sum((Odds_WR_optA$beta.cmb - beta0)^2)
#


Odds_LCC_test <-  logistic.optimal.subsampling(formula, data, n.plt, n.ssp,
                                               family = family$family.name,
                                                     criterion = 'LCC',
                                                     sampling.method = 'Poisson', # WithReplacement
                                                     estimate.method = 'LogOddsCorrection',
                                                     b=2)
# Odds_LCC_test$beta.cmb
# MSE_Odds_LCC_test <- sum((Odds_LCC_test$beta.cmb - beta0)^2)
#
#
#
# Wet_LCC_test <-  logistic.optimal.subsampling(X.full, Y.full, n.plt, n.ssp,
#                                                criterion = 'LCC',
#                                                sampling.method = 'Poisson', # WithReplacement
#                                                estimate.method = 'Weighted',
#                                                b=2)
# Wet_LCC_test$beta.cmb
# MSE_Wet_LCC_test <- sum((Wet_LCC_test$beta.cmb - beta0)^2)












# source('R/rare_event_functions.R')
# We_WRep_optA <- RareLogistic(X.full, Y.full, n.plt, n.ssp, method = 'Weighted', criterion = 'optA')
# We_WRep_optA$beta.cmb
# MSE <- sum((We_WRep_optA$beta.cmb - beta0)^2)

We_WRep_optL_test <-  logistic.optimal.subsampling(X.full, Y.full, n.plt, n.ssp,
                                              criterion = 'optL',
                                              method = 'Weighted',
                                              b = 2)
We_WRep_optL_test$beta.cmb
MSE_We_WRep_optL_test <- sum((We_WRep_optL_test$beta.cmb - beta0)^2)

We_WRep_LCC_test <-  logistic.optimal.subsampling(X.full, Y.full, n.plt, n.ssp,
                                              criterion = 'LCC',
                                              method = 'Weighted',
                                              b = 2)
We_WRep_LCC_test$beta.cmb
MSE_We_WRep_LCC_test <- sum((We_WRep_LCC_test$beta.cmb - beta0)^2)

############
Unwe_Poi_optA_test <-  logistic.optimal.subsampling(X.full, Y.full, n.plt, n.ssp,
                                              criterion = 'optA',
                                              method = 'LogOddsCorrection',
                                              b = 2)
Unwe_Poi_optA_test$beta.cmb
MSE_Unwe_Poi_optA_test <- sum((Unwe_Poi_optA_test$beta.cmb - beta0)^2)

Unwe_Poi_optL_test <-  logistic.optimal.subsampling(X.full, Y.full, n.plt, n.ssp,
                                               criterion = 'optL',
                                               method = 'LogOddsCorrection',
                                               b = 2)
Unwe_Poi_optL_test$beta.cmb
MSE_Unwe_Poi_optL_test <- sum((Unwe_Poi_optL_test$beta.cmb - beta0)^2)

###############
Unwe_Poi_LCC <-  logistic.optimal.subsampling(X.full, Y.full, n.plt, n.ssp,
                                               criterion = 'LCC',
                                               method = 'Poisson',
                                               b = 2)
# LogOddsCorrection.optL.results <- RareLogistic(X.full, Y.full, n.plt, n.ssp, method = 'LogOddsCorrection', criterion = 'optL')
# Weighted.optA.results <- RareLogistic(X.full, Y.full, n.plt, n.ssp, method = 'Weighted', criterion = 'optA')
# LogOddsCorrection.optA.results <- RareLogistic(X.full, Y.full, n.plt, n.ssp, method = 'LogOddsCorrection', criterion = 'optA')
#
# Betas[,1] <- logistic.coef.estimate(X=X.full, Y=Y.full)$beta #beta.full
# Betas[,2] <- RareLogistic(X.full, Y.full, n.plt, n.ssp, method = 'UniW')$beta.ssp
# Betas[,3] <- RareLogistic(X.full, Y.full, n.plt, n.ssp, method = 'Uni')$beta.ssp
# Betas[,4] <- RareLogistic(X.full, Y.full, n.plt, n.ssp, criterion = 'LCC')$beta.ssp
#
# Betas[,5] <- Weighted.optL.results$beta.ssp
# Betas[,6] <- Weighted.optA.results$beta.ssp
# Betas[,7] <- LogOddsCorrection.optL.results$beta.ssp
# Betas[,8] <- LogOddsCorrection.optA.results$beta.ssp
#
# Betas[,9] <- Weighted.optL.results$beta.cmb
# Betas[,10] <- Weighted.optA.results$beta.cmb
# Betas[,11] <- LogOddsCorrection.optL.results$beta.cmb
# Betas[,12] <- LogOddsCorrection.optA.results$beta.cmb


MSE <- as.data.frame(apply((Betas - beta0)^2 ,2,sum))
colnames(MSE) <- 'MSE'
rownames(MSE) <- c('full', 'UniW', 'Uni', 'LCC', 'Weighted+optL', 'Weighted+optA', 'LOC+optL', 'LOC+optA',
                   'Weighted+optL+cmb', 'Weighted+optA+cmb', 'LOC+optL+cmb', 'LOC+optA+cmb')
BETAS <- as.data.frame(cbind(beta0, Betas), row.names =c(1:d))
colnames(BETAS) <- c('true', 'full', 'UniW', 'Uni', 'LCC', 'Weighted+optL', 'Weighted+optA', 'LOC+optL', 'LOC+optA',
                     'Weighted+optL+cmb', 'Weighted+optA+cmb', 'LOC+optL+cmb', 'LOC+optA+cmb')
rownames(BETAS) <- c(1:d)
MSE
log(MSE)
BETAS

### test
criterion = 'optA'
method = 'SwR'
weighted.estimator = T
b = 2
N <- nrow(X.full)
d <- ncol(X.full)
N1 <- sum(Y.full)
n.plt <- N - N1
# pilot step
if (weighted.estimator == T){
  index.plt <- uni.indx(N, Y.full, n.plt) #first n.plt/2 comes from Y==0, last n.plt/2 comes from y==1,
  pinv.plt <- 1/c(rep(1/(2 * n.plt), n.plt/2),
                  rep(1/(2 * N1), n.plt/2)) #pinv means pi's inv
  beta.plt <- logistic.coef.estimate(X.full[index.plt, ], Y.full[index.plt], weights = pinv.plt)$beta
  P.plt <- pbeta(X.full, beta.plt)
  MN.plt <- MN(X.full[index.plt, ], P.plt[index.plt], pinv.plt)
  Psi.plt <- Psi(X.full[index.plt, ], Y.full[index.plt], P.plt[index.plt], pinv.plt)
  cc.plt <- n.plt / pinv.plt[index.plt]
  return(list(
    beta.plt = beta.plt,
    MN.plt = MN.plt,
    Psi.plt = Psi.plt,
    P.plt = P.plt,
    index.plt = index.plt,
    cc.plt = cc.plt
  ))
}
#
if (criterion == "optA"){
  dm <- sqrt((Y.full - P.plt)^2 * rowSums((X.full %*% solve(MN.plt)) ^ 2))
} else if (criterion == "optL"){
  dm <- sqrt((Y.full - P.plt)^2 * rowSums(X.full^2))
} else if (criterion == "LCC"){
  dm <- abs(Y - pbeta_pilot)
}
# # subsample estimate step
# ssp.estimate.reults <- subsample.estimate(x.ssp, y.ssp, w.ssp, beta.plt, criterion, method)
# beta.ssp <- ssp.estimate.reults$beta.ssp
# MN.ssp <- ssp.estimate.reults$MN.ssp
# Psi.ssp <- ssp.estimate.reults$Psi.ssp
# #
# if (criterion == 'LCC'){
#   beta.ssp <-
#     logistic.coef.estimate(X = x.ssp, Y = y.ssp, weights = w.ssp)$beta + beta.plt
# } else if (method == "Weighted"){
#   beta.ssp <-
#     logistic.coef.estimate(X = x.ssp, Y = y.ssp, weights = 1 / w.ssp)$beta
# } else if (method == "LogOddsCorrection"){
#   p.ssp <- subsampling.results$p.ssp
#   beta.ssp <-
#     logistic.coef.estimate(
#       X = x.ssp,
#       Y = y.ssp, #start balue
#       weights = rep(1, length(y.ssp)),
#       offset = -log(p.ssp)
#     )$beta
#   # beta.ssp <- getMSLE(x,y,p,beta.plt)$beta
# }
# P.ssp  <- pbeta(x.ssp, beta.ssp)
# MN.ssp <- MN(x.ssp, P.ssp, 1) #1/w.ssp
# Psi.ssp <- Psi(x.ssp, y.ssp, P.ssp, 1)
#logistic regression
set.seed(1)
N <- 1e4
beta0 <- rep(-0.5, 7)
d <- length(beta0) - 1
X <- matrix(0, N, d)
generate_rexp <- function(x) x <- rexp(N, rate = 2)
X <- apply(X, 2, generate_rexp)
Y <- rbinom(N, 1, 1 - 1 / (1 + exp(beta0[1] + X %*% beta0[-1])))
print(paste('N: ', N))
print(paste('sum(Y): ', sum(Y)))
data <- as.data.frame(cbind(Y, X))
formula <- Y ~ .
n.plt <- 200
n.ssp <- 600
subsampling.results <- glm.optimal.subsampling(formula, data, n.plt, n.ssp, family = 'binomial', criterion = "OptL", sampling.method = 'Poisson', estimate.method = "LogOddsCorrection")
subsampling.summary(subsampling.results)
subsampling.results <- glm.optimal.subsampling(formula, data, n.plt, n.ssp, family = 'binomial', criterion = "OptL", sampling.method = 'WithReplacement', estimate.method = "Weighted")
subsampling.summary(subsampling.results)
Uni.subsampling.results <- glm.optimal.subsampling(formula, data, n.plt, n.ssp, family = 'binomial', estimate.method = 'Uni')
subsampling.summary(Uni.subsampling.results)

