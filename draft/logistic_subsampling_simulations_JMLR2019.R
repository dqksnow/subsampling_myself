library(parallel)
library(ggplot2)
library(reshape2)
rm(list = ls())
source('R/glm_subsampling_functions.R')
# source('draft/logistic_subsampling_functions_glm.R')
# source('R/cmb_sample_logistic_subsampling_functions.R')
source('draft/generate_logit_data.R')
title <- paste('MSE, Case 1, WithReplacement') # Poisson
var_title <- paste('Var of combined estimator with MLE vs with True, Case 1, SampleMethod WithReplacement')
mse_var_title <- paste('Mse vs Var, Case 1, WithReplacement, Combined estimator')

logistic_simu <- function(j) { #this function is used for parLapply
  source('R/glm_subsampling_functions.R')
  source('R/family_expand.R')
  # source('draft/logistic_subsampling_functions_glm.R')
  # source('R/cmb_sample_logistic_subsampling_functions.R')
  # source('draft/combine_sample_test.R')
  source('draft/generate_logit_data.R')
  beta.ssp <- beta.cmb <- array(data = NA, dim = c(length(n.ssp), num.method, d))
  var.ssp <- var.ssp.true <- var.cmb <- var.cmb.true <- SubsampleSize <- array(data = NA, dim = c(length(n.ssp), num.method-1))
  FULL.DATA <- simu_mzNormal(seed = j, N = N, beta0 = beta0, corr = 0.5)      # simu 1, multinormal, balanced
  # FULL.DATA <- simu_nzNormal(seed = j, N = N, beta0 = beta0, corr = 0.5)    # simu 2, multinormal, imbalanced
  # FULL.DATA <- simu_ueNormal(seed = j, N = N, beta0 = beta0, corr = 0.5)    # simu 3, multinormal, Hetero
  # FULL.DATA <- simu_mixNormal(seed = j, N = N, beta0 = beta0, corr = 0.5)   # simu 4, mix normal
  # FULL.DATA <- simu_T3(seed = j, N = N, beta0 = beta0, corr = 0.5, df = 3)  # simu 5, T3
  # FULL.DATA <- simu_EXP(seed = j, N = N, beta0 = beta0, rate = 2)           # simu 6, #EXP
  family <- binomial.expand()
  X.full <- FULL.DATA$X
  Y.full <- FULL.DATA$Y
  data <- as.data.frame(cbind(Y.full, X.full))
  # colnames(data) <- c("response", "intercept", paste0("v", 1:6))
  # formula <- response ~ . -1
  formula <- Y.full ~ .
  beta.full <- glm.coef.estimate(X=cbind(1, X.full), Y=Y.full, family = family)$beta #beta.full
  for(i in seq_along(n.ssp)) {
    Wet.optA <- glm.optimal.subsampling(formula, data, n.plt, n.ssp[i],
                                             criterion = 'optA',
                                             sampling.method = 'Poisson', # Poisson; WithReplacement
                                             estimate.method = 'Weighted',
                                             alpha = 0,
                                             b = 2,
                                             family = family$family.name)
    beta.ssp[i, 1, ] <- Wet.optA$beta.ssp
    var.ssp[i, 1] <- sum(diag(Wet.optA$var.ssp))
    var.ssp.true[i, 1] <- sum(diag(Wet.optA$var.ssp.true))
    beta.cmb[i, 1, ] <- Wet.optA$beta.cmb
    var.cmb[i, 1] <- sum(diag(Wet.optA$var.cmb))
    var.cmb.true[i, 1] <- sum(diag(Wet.optA$var.cmb.true))
    # SubsampleSize[i, 1] <- length(Wet.optA$index.ssp)
    #
    Wet.optL <- glm.optimal.subsampling(formula, data, n.plt, n.ssp[i],
                                                         criterion = 'optL',
                                                    sampling.method = 'Poisson',
                                                    estimate.method = 'Weighted',
                                             alpha = 0,
                                                         b = 2,
                                             family = family$family.name)
    beta.ssp[i, 2, ] <- Wet.optL$beta.ssp
    var.ssp[i, 2] <- sum(diag(Wet.optL$var.ssp))
    var.ssp.true[i, 2] <- sum(diag(Wet.optL$var.ssp.true))
    beta.cmb[i, 2, ] <- Wet.optL$beta.cmb
    var.cmb[i, 2] <- sum(diag(Wet.optL$var.cmb))
    var.cmb.true[i, 2] <- sum(diag(Wet.optL$var.cmb.true))
    #
    Wet.LCC <- glm.optimal.subsampling(formula, data, n.plt, n.ssp[i],
                                                         criterion = 'LCC',
                                                   sampling.method = 'Poisson',
                                                   estimate.method = 'Weighted',
                                            alpha = 0,
                                                         b = 2,
                                            family = family$family.name)
    beta.ssp[i, 3, ] <- Wet.LCC$beta.ssp
    var.ssp[i, 3] <- sum(diag(Wet.LCC$var.ssp))
    var.ssp.true[i, 3] <- sum(diag(Wet.LCC$var.ssp.true))
    beta.cmb[i, 3, ] <- Wet.LCC$beta.cmb
    var.cmb[i, 3] <- sum(diag(Wet.LCC$var.cmb))
    var.cmb.true[i, 3] <- sum(diag(Wet.LCC$var.cmb.true))
    #
    Odds.optA <- glm.optimal.subsampling(formula, data, n.plt, n.ssp[i],
                                                         criterion = 'optA',
                                                     sampling.method = 'Poisson',
                                                     estimate.method = 'LogOddsCorrection',
                                                         b = 2,
                                              family = family$family.name)
    beta.ssp[i, 4, ] <- Odds.optA$beta.ssp
    var.ssp[i, 4] <- sum(diag(Odds.optA$var.ssp))
    var.ssp.true[i, 4] <- sum(diag(Odds.optA$var.ssp.true))
    beta.cmb[i, 4, ] <- Odds.optA$beta.cmb
    var.cmb[i, 4] <- sum(diag(Odds.optA$var.cmb))
    var.cmb.true[i, 4] <- sum(diag(Odds.optA$var.cmb.true))
    #
    Odds.optL <- glm.optimal.subsampling(formula, data, n.plt, n.ssp[i],
                                                         criterion = 'optL',
                                                     sampling.method = 'Poisson',
                                                     estimate.method = 'LogOddsCorrection',
                                                         b = 2,
                                              family = family$family.name)
    beta.ssp[i, 5, ] <- Odds.optL$beta.ssp
    var.ssp[i, 5] <- sum(diag(Odds.optL$var.ssp))
    var.ssp.true[i, 5] <- sum(diag(Odds.optL$var.ssp.true))
    beta.cmb[i, 5, ] <- Odds.optL$beta.cmb
    var.cmb[i, 5] <- sum(diag(Odds.optL$var.cmb))
    var.cmb.true[i, 5] <- sum(diag(Odds.optL$var.cmb.true))
    #
    Odds.LCC <- glm.optimal.subsampling(formula, data, n.plt, n.ssp[i],
                                                         criterion = 'LCC',
                                                    sampling.method = 'Poisson',
                                                    estimate.method = 'LogOddsCorrection',
                                                         b = 2,
                                             family = family$family.name)
    beta.ssp[i, 6, ] <- Odds.LCC$beta.ssp
    var.ssp[i, 6] <- sum(diag(Odds.LCC$var.ssp))
    var.ssp.true[i, 6] <- sum(diag(Odds.LCC$var.ssp.true))
    beta.cmb[i, 6, ] <- Odds.LCC$beta.cmb
    var.cmb[i, 6] <- sum(diag(Odds.LCC$var.cmb))
    var.cmb.true[i, 6] <- sum(diag(Odds.LCC$var.cmb.true))
    #
    Uni <- glm.optimal.subsampling(formula, data, n.plt, n.ssp[i],
                                             estimate.method = 'Uni',
                                             b = 2,
                                        family = family$family.name)
    beta.ssp[i, 7, ] <- beta.cmb[i, 7, ] <- Uni$beta
  }
  return(list(beta.ssp = beta.ssp,
              beta.cmb = beta.cmb,
              beta.full = beta.full,
              var.ssp = var.ssp,
              var.ssp.true = var.ssp.true,
              var.cmb = var.cmb,
              var.cmb.true = var.cmb.true,
              SubsampleSize = SubsampleSize))
}

rpt <- 100
N <-  1e4
beta0 <- c(rep(0.5, 7))
d <- length(beta0)
n.plt <- 200
# n.ssp_min <- 200
# n.ssp_max <- 3000
# n.ssp <- c(100, seq(n.ssp_min, n.ssp_max, 400))
n.ssp <- c(100, 200, 500, 1000)
method.all <- c('Wet.optA', 'Wet.optL', 'Wet.LCC', 'Odds.optA', 'Odds.optL', 'Odds.LCC', 'Uni')
num.method <- length(method.all)
Betas.ssp <- array(data = NA, dim = c(rpt, length(n.ssp), num.method, d))
Betas.cmb <- array(data = NA, dim = c(rpt, length(n.ssp), num.method, d))

beta.full <- matrix(NA, rpt, d)
var.ssp <- var.ssp.true <- var.cmb <- var.cmb.true <- array(data = NA, dim = c(rpt, length(n.ssp), num.method-1))
SubsampleSize <- array(data = NA, dim = c(rpt, length(n.ssp), num.method-1))


cl <- makeCluster(8) # how many CPU cores are called
t1 <- proc.time()
clusterExport(cl=cl,
              varlist=c('n.plt', 'n.ssp', 'beta0', 'num.method', 'd', 'N'),
              envir=environment())
results <- parLapply(cl, 1:rpt, logistic_simu)
t2 <- proc.time()
time.cost <- t2 - t1
for(i in 1:rpt){
  Betas.ssp[i,,,] <- results[[i]][[1]]
  Betas.cmb[i,,,] <- results[[i]][[2]]
  beta.full[i,] <- results[[i]][[3]]
  var.ssp[i,,] <- results[[i]][[4]]
  var.ssp.true[i,,] <- results[[i]][[5]]
  var.cmb[i,,] <- results[[i]][[6]]
  var.cmb.true[i,,] <- results[[i]][[7]]
  SubsampleSize[i,,] <- results[[i]][[8]]
}

MeanVar.ssp <- apply(var.ssp, c(2, 3), mean)
MeanVar.ssp.true <- apply(var.ssp.true, c(2, 3), mean)
MeanVar.cmb <- apply(var.cmb, c(2, 3), mean)
MeanVar.cmb.true <- apply(var.cmb.true, c(2, 3), mean)
# MeanSubsampleSize <- apply(SubsampleSize, c(2, 3), mean)
# MeanSubsampleSize <- cbind(n.ssp, MeanSubsampleSize)
mse.ssp <- mse.cmb <- array(data = NA, dim = c(length(n.ssp), num.method))
for (idn in seq_along(n.ssp)){
  for (m in 1:num.method) {
    mse.ssp[idn, m] <- mean(apply(sweep(Betas.ssp[, idn, m,], 2, beta0) ^ 2 , 1, sum))
    mse.cmb[idn, m] <- mean(apply(sweep(Betas.cmb[, idn, m,], 2, beta0) ^ 2 , 1, sum))
    # fname <- paste0("rare_event_simu_output/csv/case",case,"N",N,"method",m,"n",n,".csv")
    # write.csv(Betas[,idn,,m],file=fname,row.names=FALSE)
  }
}

# alpha_beta_diff_ss <- array(NA, dim = c(rpt, length(n.ssp), num.method))
#
# for (method in 1:num.method) {
#   for (i in 1:length(n.ssp)) {
#     alpha_beta_diff_ss[, i, method] <- rowSums((beta.full - Betas.cmb[, i, method, ])^2)
#   }
# }
#
# mean_alpha_beta_diff_ss <- apply(alpha_beta_diff_ss, c(2, 3), mean)  # Mean for each setting and method



mse.full <- mean(apply(sweep(beta.full, 2, beta0) ^ 2 , 1, sum))
# mse <- log(mse)
# mse.full <- log(mse.full)

# plot: mse ###################################################################
mse_plot <- data.frame(rou = n.ssp / N,
                       Wet.optA = mse.cmb[, 1],
                       Wet.optL = mse.cmb[, 2],
                       Wet.LCC = mse.cmb[, 3],
                       Odds.optA = mse.cmb[, 4],
                       Odds.optL = mse.cmb[, 5],
                       Odds.LCC = mse.cmb[, 6],
                       Uni = mse.ssp[, 7]
)
mse_plot<-melt(mse_plot, id.vars = 'rou', variable.name = "method",
               value.name="mse") # convert data structure to fit ggplot
mse_figure <- ggplot(mse_plot, aes(x = rou, y = mse)) +
  geom_hline(aes(yintercept = mse.full), linetype = "dashed") +
  geom_line(aes(colour = method, linetype = method)) + geom_point() + ggtitle(title) +
  scale_x_continuous(breaks = n.ssp / N) +
  theme_set(theme_bw()) +
  scale_color_manual(
    values = c(
      'Wet.optA' = "red",
      'Wet.optL' = "blue",
      'Wet.LCC' = "green",
      'Odds.optA' = "red",
      'Odds.optL' = "blue",
      'Odds.LCC' = "green",
      'Uni' = "gray"
    )) +
  scale_linetype_manual(
    values = c(
      'Wet.optA' = "solid",
      'Wet.optL' = "solid",
      'Wet.LCC' = "solid",
      'Odds.optA' = "dashed",
      'Odds.optL' = "dashed",
      'Odds.LCC' = "dashed",
      'Uni' = "solid"
    )
  )
# mse_figure
# print(paste0('total time: ',time.cost[3][[1]],'s'))

# plot: estimate of var, for full estimator and for true parameter ############
# var_plot <- data.frame(rou = n.ssp / N,
#                        Wet.optA.MLE = MeanVar.cmb[, 1],
#                        Wet.optL.MLE = MeanVar.cmb[, 2],
#                        Wet.LCC.MLE = MeanVar.cmb[, 3],
#                        Wet.optA.true = MeanVar.cmb.true[, 1],
#                        Wet.optL.true = MeanVar.cmb.true[, 2],
#                        Wet.LCC.true = MeanVar.cmb.true[, 3]
# )
# var_plot <- melt(var_plot, id.vars = 'rou', variable.name="method",
#                value.name="var") # convert data structure to fit ggplot
# var_figure <- ggplot(var_plot, aes(x = rou, y = var)) +
#   geom_line(aes(colour = method, linetype = method)) + geom_point() + ggtitle(var_title) +
#   scale_x_continuous(breaks = n.ssp / N) +
#   theme_set(theme_bw()) +
#   scale_color_manual(
#     values = c(
#       'Wet.optA.MLE' = "red",
#       'Wet.optL.MLE' = "blue",
#       'Wet.LCC.MLE' = "green",
#       'Wet.optA.true' = "red",
#       'Wet.optL.true' = "blue",
#       'Wet.LCC.true' = "green"
#     )) +
#   scale_linetype_manual(
#     values = c(
#       'Wet.optA.MLE' = "solid",
#       'Wet.optL.MLE' = "solid",
#       'Wet.LCC.MLE' = "solid",
#       'Wet.optA.true' = "dashed",
#       'Wet.optL.true' = "dashed",
#       'Wet.LCC.true' = "dashed"
#     )
#   )

# plot: compare of mse and var.true ###########################################
mse_vs_2var_plot <- data.frame(rou = n.ssp / N,
                              Odds.optA.MSE = mse.cmb[, 5],
                              Odds.optA.Var.MLE = MeanVar.cmb[, 5],
                              Odds.optA.Var.true = MeanVar.cmb.true[, 5]) # MeanVar.ssp.true[, 5]
mse_vs_2var_plot<-melt(mse_vs_2var_plot, id.vars = 'rou', variable.name="method",
                      value.name="var") # convert data structure to fit ggplot
mse_vs_2var_figure <- ggplot(mse_vs_2var_plot, aes(x = rou, y = var)) +
  geom_line(aes(colour = method, linetype = method)) + geom_point() + ggtitle(mse_var_title) +
  geom_hline(aes(yintercept = mse.full), linetype = "dashed") +
  scale_x_continuous(breaks = n.ssp / N) +
  theme_set(theme_bw()) +
  scale_color_manual(
    values = c(
      'Odds.optA.MSE' = "red",
      'Odds.optA.Var.MLE' = "blue",
      'Odds.optA.Var.true' = "green"
      )) +
  scale_linetype_manual(
    values = c(
      'Odds.optA.MSE' = "solid",
      'Odds.optA.Var.MLE' = "dashed",
      'Odds.optA.Var.true' = "dashed"
    )
  )

print(paste0('total time: ',time.cost[3][[1]],'s'))
# mse
mse_figure
# var_figure
# mse_vs_var_figure
mse_vs_2var_figure
# MeanSubsampleSize
# ggsave(paste('draft/0510/', title, '.png'), width = 6, height = 3, dpi = 300)



###############################################################################
###################
# single core test
# set.seed(1)

#
# for(i in seq_along(n.ssp)){
#   We.SwR.optA <- logistic.optimal.subsampling(X.full, Y.full, n.plt, n.ssp[i],
#                                                        criterion = 'optA',
#                                                        method = 'SwR',
#                                                        weighted.estimator = T,
#                                                        b = 2)
#   beta.1rep[i, 1, ] <- We.SwR.optA$beta.cmb
#   var[i, 1] <- sum(diag(We.SwR.optA$var.beta))
#   #
#   We.SwR.optL <- logistic.optimal.subsampling(X.full, Y.full, n.plt, n.ssp[i],
#                                                        criterion = 'optL',
#                                                        method = 'SwR',
#                                                        weighted.estimator = T,
#                                                        b = 2)
#   beta.1rep[i, 2, ] <- We.SwR.optL$beta.cmb
#   var[i, 2] <- sum(diag(We.SwR.optL$var.beta))
#   #
#   We.SwR.LCC <- logistic.optimal.subsampling(X.full, Y.full, n.plt, n.ssp[i],
#                                                       criterion = 'LCC',
#                                                       method = 'SwR',
#                                                       weighted.estimator = T,
#                                                       b = 2)
#   beta.1rep[i, 3, ] <- We.SwR.LCC$beta.cmb
#   var[i, 3] <- sum(diag(We.SwR.LCC$var.beta))
#   #
#   LOC.optA <- logistic.optimal.subsampling(X.full, Y.full, n.plt, n.ssp[i],
#                                                         criterion = 'optA',
#                                                         method = 'Poisson',
#                                                         weighted.estimator = F,
#                                                         b = 2)
#   beta.1rep[i, 4, ] <- LOC.optA$beta.cmb
#   var[i, 4] <- sum(diag(LOC.optA$var.beta))
#   #
#   We.SwR.optL <- logistic.optimal.subsampling(X.full, Y.full, n.plt, n.ssp[i],
#                                                        criterion = 'optL',
#                                                        method = 'Poisson',
#                                                        weighted.estimator = F,
#                                                        b = 2)
#   beta.1rep[i, 5, ] <- We.SwR.optL$beta.cmb
#   var[i, 5] <- sum(diag(We.SwR.optL$var.beta))
#   #
#   We.SwR.LCC <- logistic.optimal.subsampling(X.full, Y.full, n.plt, n.ssp[i],
#                                                       criterion = 'LCC',
#                                                       method = 'Poisson',
#                                                       weighted.estimator = F,
#                                                       b = 2)
#   beta.1rep[i, 6, ] <- We.SwR.LCC$beta.cmb
#   var[i, 6] <- sum(diag(We.SwR.LCC$var.beta))
# }
##
