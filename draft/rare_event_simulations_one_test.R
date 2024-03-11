library(parallel)
library(ggplot2)
library(reshape2)
source('R/rare_event_functions.R')
source('draft/generate_rare_logit_data.R')
rpt <- 1
N <-  5 * 1e5
beta0 <- c(NA, -rep(1, 6))
case <- 1 # case 3 might get error
beta0 <- generate.rare.data(N, case = case, beta0)$beta0
d <- length(beta0)
n.plt <- 200
n.ssp <- c(1000, 2000, 5000, 10000)
# n.ssp_min <- 200
# n.ssp_max <- 1000
# n.ssp <- c(100, seq(n.ssp_min, n.ssp_max, 200))
method.all <- c('Wet_optA', 'Wet_optL', 'Wet_LCC', 'Odds_optA', 'Odds_optL','Odds_LCC', 'Uni')
num.method <- length(method.all)
Betas <- array(data = NA, dim = c(rpt, length(n.ssp), num.method, d))
beta.full <- matrix(NA, rpt, d)
SubsampleSize <- var.ssp <- var.cmb <- array(data = NA, dim = c(rpt, length(n.ssp), num.method-1))
j <- 1
set.seed(j)
beta.1rep <- array(data = NA, dim = c(length(n.ssp), num.method, d))
SubsampleSize <- var.ssp <- var.cmb <- array(data = NA, dim = c(length(n.ssp), (num.method-1)))
FULL.DATA <- generate.rare.data(N, case = case, beta0)
X.full <- FULL.DATA$X
Y.full <- FULL.DATA$Y
N1 <- sum(Y.full == 1)
beta0 <- FULL.DATA$beta0
beta.full <- logistic.coef.estimate(X=X.full, Y=Y.full)$beta #beta.full
alpha <- 0.1
for(i in seq_along(n.ssp)){
  n1 <- n.ssp[i]

  Wet.optA <- RareLogistic(X.full, Y.full, n.plt, n1, estimate.method = 'Weighted', criterion = 'optA', alpha = alpha)
  beta.1rep[i, 1, ] <- Wet.optA$beta.cmb
  SubsampleSize[i, 1] <- length(Wet.optA$index.ssp) - N1
  var.cmb[i, 1] <- sum(diag(Wet.optA$var.cmb))
  var.ssp[i, 1] <- sum(diag(Wet.optA$var.ssp))

  Wet.optL <- RareLogistic(X.full, Y.full, n.plt, n1, estimate.method = 'Weighted', criterion = 'optL', alpha = alpha)
  beta.1rep[i, 2, ] <- Wet.optL$beta.cmb
  SubsampleSize[i, 2] <- length(Wet.optL$index.ssp) - N1
  var.cmb[i, 2] <- sum(diag(Wet.optL$var.cmb))
  var.ssp[i, 2] <- sum(diag(Wet.optL$var.ssp))

  Wet.LCC <- RareLogistic(X.full, Y.full, n.plt, n1, estimate.method = 'Weighted', criterion = 'LCC', alpha = alpha)
  beta.1rep[i, 3, ] <- Wet.LCC$beta.cmb
  SubsampleSize[i, 3] <- length(Wet.LCC$index.ssp) - N1
  var.cmb[i, 3] <- sum(diag(Wet.LCC$var.cmb))
  var.ssp[i, 3] <- sum(diag(Wet.LCC$var.ssp))

  Odds.optA <- RareLogistic(X.full, Y.full, n.plt, n1, estimate.method = 'LogOddsCorrection', criterion = 'optA', alpha = alpha)
  beta.1rep[i, 4, ] <- Odds.optA$beta.cmb
  SubsampleSize[i, 4] <- length(Odds.optA$index.ssp) - N1
  var.cmb[i, 4] <- sum(diag(Odds.optA$var.cmb))
  var.ssp[i, 4] <- sum(diag(Odds.optA$var.ssp))

  Odds.optL <- RareLogistic(X.full, Y.full, n.plt, n1, estimate.method = 'LogOddsCorrection', criterion = 'optL', alpha = alpha)
  beta.1rep[i, 5, ] <- Odds.optL$beta.cmb
  SubsampleSize[i, 5] <- length(Odds.optL$index.ssp) - N1
  var.cmb[i, 5] <- sum(diag(Odds.optL$var.cmb))
  var.ssp[i, 5] <- sum(diag(Odds.optL$var.ssp))

  Odds.LCC <- RareLogistic(X.full, Y.full, n.plt, n1, estimate.method = 'LogOddsCorrection', criterion = 'LCC', alpha = alpha)
  beta.1rep[i, 6, ] <- Odds.LCC$beta.cmb
  SubsampleSize[i, 6] <- length(Odds.LCC$index.ssp) - N1
  var.cmb[i, 6] <- sum(diag(Odds.LCC$var.cmb))
  var.ssp[i, 6] <- sum(diag(Odds.LCC$var.ssp))

  beta.1rep[i, 7, ] <- RareLogistic(X.full, Y.full, n.plt, n1, estimate.method = 'Uni')$beta.ssp
}
