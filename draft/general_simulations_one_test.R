library(parallel)
library(ggplot2)
library(reshape2)
rm(list = ls())
source('R/logistic_subsampling_functions.R')
# source('draft/combine_sample_test.R')
source('draft/generate_logit_data.R')
rpt <- 1
N <-  1e4
beta0 <- c(rep(0.5, 7))
d <- length(beta0)
n.plt <- 200
n.ssp_min <- 200
n.ssp_max <- 1000
n.ssp <- c(100, seq(n.ssp_min, n.ssp_max, 200))
method.all <- c('Wet.optA', 'Wet.optL', 'Wet.LCC', 'Odds.optA', 'Odds.optL', 'Odds.LCC', 'Uni')
num.method <- length(method.all)
Betas <- array(data = NA, dim = c(rpt, length(n.ssp), num.method, d))
beta.full <- matrix(NA, rpt, d)
var.cmb <- var.cmb.true <- array(data = NA, dim = c(rpt, length(n.ssp), num.method-1))
SubsampleSize <- array(data = NA, dim = c(rpt, length(n.ssp), num.method-1))
j <- 1
beta.1rep <- array(data = NA, dim = c(length(n.ssp), num.method, d))
var.cmb <- var.cmb.true <- SubsampleSize <- array(data = NA, dim = c(length(n.ssp), num.method-1))
# FULL.DATA <- simu_mzNormal(seed = j, N = N, beta0 = beta0, corr = 0.5)      # simu 1, multinormal, balanced
FULL.DATA <- simu_nzNormal(seed = j, N = N, beta0 = beta0, corr = 0.5)    # simu 2, multinormal, imbalanced
# FULL.DATA <- simu_ueNormal(seed = j, N = N, beta0 = beta0, corr = 0.5)    # simu 3, multinormal, Hetero
# FULL.DATA <- simu_mixNormal(seed = j, N = N, beta0 = beta0, corr = 0.5)   # simu 4, mix normal
# FULL.DATA <- simu_T3(seed = j, N = N, beta0 = beta0, corr = 0.5, df = 3)  # simu 5, T3
# FULL.DATA <- simu_EXP(seed = j, N = N, beta0 = beta0, rate = 2)           # simu 6, #EXP
X.full <- FULL.DATA$X
Y.full <- FULL.DATA$Y
data <- as.data.frame(cbind(Y.full, X.full))
formula <- Y.full ~ . -1
beta.full <- logistic.coef.estimate(X=X.full, Y=Y.full)$beta #beta.full
for(i in seq_along(n.ssp)) {
  Wet.optA <- logistic.optimal.subsampling(formula, data, n.plt, n.ssp[i],
                                           criterion = 'optA',
                                           sampling.method = 'WithReplacement', # Poisson; WithReplacement
                                           estimate.method = 'Weighted',
                                           b = 2)
  beta.1rep[i, 1, ] <- Wet.optA$beta.cmb
  var.cmb[i, 1] <- sum(diag(Wet.optA$var.cmb))
  var.cmb.true[i, 1] <- sum(diag(Wet.optA$var.cmb.true))
  SubsampleSize[i, 1] <- length(Wet.optA$index.ssp)
  #
  Wet.optL <- logistic.optimal.subsampling(formula, data, n.plt, n.ssp[i],
                                           criterion = 'optL',
                                           sampling.method = 'WithReplacement',
                                           estimate.method = 'Weighted',
                                           b = 2)
  beta.1rep[i, 2, ] <- Wet.optL$beta.cmb
  var.cmb[i, 2] <- sum(diag(Wet.optL$var.cmb))
  var.cmb.true[i, 2] <- sum(diag(Wet.optL$var.cmb.true))
  SubsampleSize[i, 2] <- length(Wet.optL$index.ssp)
  #
  Wet.LCC <- logistic.optimal.subsampling(formula, data, n.plt, n.ssp[i],
                                          criterion = 'LCC',
                                          sampling.method = 'WithReplacement',
                                          estimate.method = 'Weighted',
                                          b = 2)
  beta.1rep[i, 3, ] <- Wet.LCC$beta.cmb
  var.cmb[i, 3] <- sum(diag(Wet.LCC$var.cmb))
  var.cmb.true[i, 3] <- sum(diag(Wet.LCC$var.cmb.true))
  SubsampleSize[i, 3] <- length(Wet.LCC$index.ssp)
  #
  Odds.optA <- logistic.optimal.subsampling(formula, data, n.plt, n.ssp[i],
                                            criterion = 'optA',
                                            sampling.method = 'WithReplacement',
                                            estimate.method = 'LogOddsCorrection',
                                            b = 2)
  beta.1rep[i, 4, ] <- Odds.optA$beta.cmb
  var.cmb[i, 4] <- sum(diag(Wet.optA$var.cmb))
  var.cmb.true[i, 4] <- sum(diag(Wet.optA$var.cmb.true))
  SubsampleSize[i, 4] <- length(Odds.optA$index.ssp)
  #
  Odds.optL <- logistic.optimal.subsampling(formula, data, n.plt, n.ssp[i],
                                            criterion = 'optL',
                                            sampling.method = 'WithReplacement',
                                            estimate.method = 'LogOddsCorrection',
                                            b = 2)
  beta.1rep[i, 5, ] <- Odds.optL$beta.cmb
  var.cmb[i, 5] <- sum(diag(Wet.optL$var.cmb))
  var.cmb.true[i, 5] <- sum(diag(Wet.optL$var.cmb.true))
  SubsampleSize[i, 5] <- length(Odds.optL$index.ssp)
  #
  Odds.LCC <- logistic.optimal.subsampling(formula, data, n.plt, n.ssp[i],
                                           criterion = 'LCC',
                                           sampling.method = 'WithReplacement',
                                           estimate.method = 'LogOddsCorrection',
                                           b = 2)
  beta.1rep[i, 6, ] <- Odds.LCC$beta.cmb
  var.cmb[i, 6] <- sum(diag(Wet.LCC$var.cmb))
  var.cmb.true[i, 6] <- sum(diag(Wet.LCC$var.cmb.true))
  SubsampleSize[i, 6] <- length(Odds.LCC$index.ssp)
  #
  Uni <- logistic.optimal.subsampling(formula, data, n.plt, n.ssp[i],
                                      estimate.method = 'Uni',
                                      b = 2)
  beta.1rep[i, 7, ] <- Uni$beta.ssp
}
