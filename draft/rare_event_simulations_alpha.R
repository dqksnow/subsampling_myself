library(parallel)
library(ggplot2)
library(reshape2)
rm(list = ls())
source('R/rare_event_functions.R')
source('draft/generate_rare_logit_data.R')
title <- paste('Rare, Different alpha, LogOddsCorrection,Case 1') # Poisson
mse_var_title <- paste('Mse vs Var, Case 1, Poisson, Combined estimator')

logistic_simu <- function(replication) { #this function is used for parLapply
  source('R/rare_event_functions.R')
  source('draft/generate_rare_logit_data.R')
  beta.ssp <- beta.cmb <- array(data = NA, dim = c(length(alpha.all), length(n.ssp.all), d))
  var.ssp <- var.ssp.true <- var.cmb <- var.cmb.true <- SubsampleSize <- array(data = NA, dim = c(length(alpha.all), length(n.ssp.all)))
  FULL.DATA <- generate.rare.data(N, case = case, beta0)
  X.full <- FULL.DATA$X
  Y.full <- FULL.DATA$Y
  N1 <- sum(Y.full == 1)
  data <- as.data.frame(cbind(Y.full, X.full))
  formula <- Y.full ~ .
  beta.full <- logistic.coef.estimate(X=cbind(1, X.full), Y=Y.full)$beta #beta.full
  for (i in seq_along(alpha.all)) {
    for(j in seq_along(n.ssp.all)) {
      Wet.optA <- RareLogistic(formula, data, n.plt, n.ssp.all[j], estimate.method = 'LogOddsCorrection', criterion = 'optA', alpha = alpha.all[i]) # Weighted
      beta.ssp[i, j, ] <- Wet.optA$beta.ssp
      var.ssp.true[i, j] <- sum(diag(Wet.optA$var.ssp))
      beta.cmb[i, j, ] <- Wet.optA$beta.cmb
      var.cmb.true[i, j] <- sum(diag(Wet.optA$var.cmb))
      #
    }
  }
  return(list(beta.ssp = beta.ssp,
              beta.cmb = beta.cmb,
              beta.full = beta.full,
              var.ssp.true = var.ssp.true,
              var.cmb.true = var.cmb.true))
}

rpt <- 100
N <-  5 * 1e5
beta0 <- c(NA, -rep(1, 6))
case <- 1
beta0 <- generate.rare.data(N, case = case, beta0)$beta0
d <- length(beta0)
n.plt <- 200
n.ssp.all <- c(1e3, 2e3, 5e3, 1e4)

alpha.all <- seq(0, 0.3, 0.1)
num.alpha <- length(alpha.all)
method.all <- c('Wet.optA', 'Wet.optL', 'Wet.LCC', 'Odds.optA', 'Odds.optL', 'Odds.LCC', 'Uni')
num.method <- length(method.all)
Betas.ssp <- array(data = NA, dim = c(rpt, num.alpha, length(n.ssp.all), d))
Betas.cmb <- array(data = NA, dim = c(rpt, num.alpha, length(n.ssp.all), d))

beta.full <- matrix(NA, rpt, d)
var.ssp <- var.ssp.true <- var.cmb <- var.cmb.true <- array(data = NA, dim = c(rpt, num.alpha, length(n.ssp.all)))
# SubsampleSize <- array(data = NA, dim = c(rpt, length(n.ssp.all), num.method-1))

cl <- makeCluster(8) # how many CPU cores are called
t1 <- proc.time()
clusterExport(cl=cl,
              varlist=c('n.plt', 'n.ssp.all', 'beta0', 'alpha.all', 'd', 'N', 'case'),
              envir=environment())
results <- parLapply(cl, 1:rpt, logistic_simu)
t2 <- proc.time()
time.cost <- t2 - t1
for(i in 1:rpt){
  Betas.ssp[i,,,] <- results[[i]][[1]]
  Betas.cmb[i,,,] <- results[[i]][[2]]
  beta.full[i,] <- results[[i]][[3]]
  var.ssp.true[i,,] <- results[[i]][[4]]
  var.cmb.true[i,,] <- results[[i]][[5]]
}

MeanVar.ssp.true <- apply(var.ssp.true, c(2, 3), mean)
MeanVar.cmb.true <- apply(var.cmb.true, c(2, 3), mean)
# MeanSubsampleSize <- apply(SubsampleSize, c(2, 3), mean)
# MeanSubsampleSize <- cbind(n.ssp.all, MeanSubsampleSize)
mse.ssp <- mse.cmb <- array(data = NA, dim = c(num.alpha, length(n.ssp.all)))
for (idn in 1:num.alpha){
  for (m in seq_along(n.ssp.all)) {
    mse.ssp[idn, m] <- mean(apply(sweep(Betas.ssp[, idn, m,], 2, beta0) ^ 2 , 1, sum))
    mse.cmb[idn, m] <- mean(apply(sweep(Betas.cmb[, idn, m,], 2, beta0) ^ 2 , 1, sum))
    # fname <- paste0("rare_event_simu_output/csv/case",case,"N",N,"method",m,"n",n,".csv")
    # write.csv(Betas[,idn,,m],file=fname,row.names=FALSE)
  }
}

mse.full <- mean(apply(sweep(beta.full, 2, beta0) ^ 2 , 1, sum))
# mse <- log(mse)
# mse.full <- log(mse.full)
# plot: mse ###################################################################
mse_plot <- data.frame(alpha = alpha.all, rou = rep(n.ssp.all / N, each = num.alpha), mse = as.vector(mse.cmb))
mse_figure <- ggplot(mse_plot, aes(x = rou, y = mse, group = alpha, color = alpha)) +
  geom_hline(aes(yintercept = mse.full), linetype = "dashed") +
  geom_line() + geom_point() + ggtitle(title) +
  scale_x_continuous(breaks = n.ssp.all / N) +
  labs(
    x = "n.ssp",
    y = "mse",
    color = "alpha"
  ) +
  scale_color_gradient(low = "blue", high = "red") +
  theme_set(theme_bw())
mse_figure

mse_vs_var_plot <- data.frame(rou = n.ssp.all / N, mse = as.vector(mse.cmb[1,]), var = MeanVar.cmb.true[1,])
mse_vs_var_plot <- melt(mse_vs_var_plot, id.vars = 'rou', variable.name="method", value.name="mse") # convert data structure to fit ggplot
mse_vs_var_figure <- ggplot(mse_vs_var_plot, aes(x = rou, y = mse)) +
  geom_hline(aes(yintercept = mse.full), linetype = "dashed") +
  geom_line(aes(colour = method, linetype = method)) + geom_point() + ggtitle(title) +
  scale_x_continuous(breaks = n.ssp.all / N) +
  labs(
    x = "n.ssp",
    y = "mse",
    color = "alpha"
  ) +
  scale_color_manual(
    values = c(
      'mse' = "red",
      'var' = "blue"
    ))  +
  theme_set(theme_bw())
mse_vs_var_figure
print(paste0('total time: ',time.cost[3][[1]],'s'))
