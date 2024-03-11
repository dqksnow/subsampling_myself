library(parallel)
library(ggplot2)
library(reshape2)
rm(list = ls())
source('R/glm_subsampling_functions.R')
# source('draft/logistic_subsampling_functions_glm.R')
# source('R/cmb_sample_logistic_subsampling_functions.R')
source('draft/generate_nb_data.R')
title <- paste('NB Regression, MSE, Case 4') # Poisson
var_title <- paste('Var of combined estimator with MLE vs with True, Case 1, SampleMethod WithReplacement')
mse_var_title <- paste('Mse vs Var, Case 1, WithReplacement, Combined estimator')

logistic_simu <- function(j) { #this function is used for parLapply
  source('R/glm_subsampling_functions.R')
  source('R/family_expand.R')
  # source('draft/logistic_subsampling_functions_glm.R')
  # source('R/cmb_sample_logistic_subsampling_functions.R')
  # source('draft/combine_sample_test.R')
  source('draft/generate_nb_data.R')
  beta.ssp <- beta.cmb <- array(data = NA, dim = c(length(n.ssp), num.method, d))
  var.ssp <- var.ssp.true <- var.cmb <- var.cmb.true <- SubsampleSize <- array(data = NA, dim = c(length(n.ssp), num.method-1))
  FULL.DATA <- generate.nb.data(N, case=4, beta0, seed = j)
  family <- negative.binomial.expand()
  X.full <- FULL.DATA$X
  Y.full <- FULL.DATA$Y
  data <- as.data.frame(cbind(Y.full, X.full))
  # colnames(data) <- c("response", "intercept", paste0("v", 1:6))
  # formula <- response ~ . -1
  formula <- Y.full ~ .
  beta.full <- glm.coef.estimate(X=cbind(1, X.full), Y=Y.full, family = family)$beta #beta.full
  for(i in seq_along(n.ssp)) {
    optA.poi <- glm.optimal.subsampling(formula, data, n.plt, n.ssp[i],
                                        criterion = 'optA',
                                        sampling.method = 'Poisson', # Poisson; WithReplacement
                                        estimate.method = 'Weighted',
                                        alpha = 0,
                                        b = 2,
                                        family = family$family.name)
    beta.ssp[i, 1, ] <- optA.poi$beta.ssp
    var.ssp[i, 1] <- sum(diag(optA.poi$var.ssp))
    var.ssp.true[i, 1] <- sum(diag(optA.poi$var.ssp.true))
    beta.cmb[i, 1, ] <- optA.poi$beta.cmb
    var.cmb[i, 1] <- sum(diag(optA.poi$var.cmb))
    var.cmb.true[i, 1] <- sum(diag(optA.poi$var.cmb.true))
    # SubsampleSize[i, 1] <- length(optA.poi$index.ssp)
    #
    optL.poi <- glm.optimal.subsampling(formula, data, n.plt, n.ssp[i],
                                        criterion = 'optL',
                                        sampling.method = 'Poisson',
                                        estimate.method = 'Weighted',
                                        alpha = 0,
                                        b = 2,
                                        family = family$family.name)
    beta.ssp[i, 2, ] <- optL.poi$beta.ssp
    var.ssp[i, 2] <- sum(diag(optL.poi$var.ssp))
    var.ssp.true[i, 2] <- sum(diag(optL.poi$var.ssp.true))
    beta.cmb[i, 2, ] <- optL.poi$beta.cmb
    var.cmb[i, 2] <- sum(diag(optL.poi$var.cmb))
    var.cmb.true[i, 2] <- sum(diag(optL.poi$var.cmb.true))
    #
    LCC.poi <- glm.optimal.subsampling(formula, data, n.plt, n.ssp[i],
                                       criterion = 'LCC',
                                       sampling.method = 'Poisson',
                                       estimate.method = 'Weighted',
                                       alpha = 0,
                                       b = 2,
                                       family = family$family.name)
    beta.ssp[i, 3, ] <- LCC.poi$beta.ssp
    var.ssp[i, 3] <- sum(diag(LCC.poi$var.ssp))
    var.ssp.true[i, 3] <- sum(diag(LCC.poi$var.ssp.true))
    beta.cmb[i, 3, ] <- LCC.poi$beta.cmb
    var.cmb[i, 3] <- sum(diag(LCC.poi$var.cmb))
    var.cmb.true[i, 3] <- sum(diag(LCC.poi$var.cmb.true))
    #
    optA.rep <- glm.optimal.subsampling(formula, data, n.plt, n.ssp[i],
                                        criterion = 'optA',
                                        sampling.method = 'WithReplacement',
                                        estimate.method = 'Weighted',
                                        b = 2,
                                        family = family$family.name)
    beta.ssp[i, 4, ] <- optA.rep$beta.ssp
    var.ssp[i, 4] <- sum(diag(optA.rep$var.ssp))
    var.ssp.true[i, 4] <- sum(diag(optA.rep$var.ssp.true))
    beta.cmb[i, 4, ] <- optA.rep$beta.cmb
    var.cmb[i, 4] <- sum(diag(optA.rep$var.cmb))
    var.cmb.true[i, 4] <- sum(diag(optA.rep$var.cmb.true))
    #
    optL.rep <- glm.optimal.subsampling(formula, data, n.plt, n.ssp[i],
                                        criterion = 'optL',
                                        sampling.method = 'WithReplacement',
                                        estimate.method = 'Weighted',
                                        b = 2,
                                        family = family$family.name)
    beta.ssp[i, 5, ] <- optL.rep$beta.ssp
    var.ssp[i, 5] <- sum(diag(optL.rep$var.ssp))
    var.ssp.true[i, 5] <- sum(diag(optL.rep$var.ssp.true))
    beta.cmb[i, 5, ] <- optL.rep$beta.cmb
    var.cmb[i, 5] <- sum(diag(optL.rep$var.cmb))
    var.cmb.true[i, 5] <- sum(diag(optL.rep$var.cmb.true))
    #
    LCC.rep <- glm.optimal.subsampling(formula, data, n.plt, n.ssp[i],
                                       criterion = 'LCC',
                                       sampling.method = 'WithReplacement',
                                       estimate.method = 'Weighted',
                                       b = 2,
                                       family = family$family.name)
    beta.ssp[i, 6, ] <- LCC.rep$beta.ssp
    var.ssp[i, 6] <- sum(diag(LCC.rep$var.ssp))
    var.ssp.true[i, 6] <- sum(diag(LCC.rep$var.ssp.true))
    beta.cmb[i, 6, ] <- LCC.rep$beta.cmb
    var.cmb[i, 6] <- sum(diag(LCC.rep$var.cmb))
    var.cmb.true[i, 6] <- sum(diag(LCC.rep$var.cmb.true))
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

rpt <- 1000
N <-  1e5
beta0 <- rep(0.5, 7)
d <- length(beta0)
n.plt <- 200
n.ssp_min <- 250
n.ssp_max <- 1500
n.ssp <- c(seq(n.ssp_min, n.ssp_max, 250))
# n.ssp <- c(100, 200, 500, 1000)
method.all <- c('optA.poi', 'optL.poi', 'LCC.poi', 'optA.rep', 'optL.rep', 'LCC.rep', 'Uni')
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
                       optA.poi = mse.cmb[, 1],
                       optL.poi = mse.cmb[, 2],
                       LCC.poi = mse.cmb[, 3],
                       optA.rep = mse.cmb[, 4],
                       optL.rep = mse.cmb[, 5],
                       LCC.rep = mse.cmb[, 6],
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
      'optA.poi' = "red",
      'optL.poi' = "blue",
      'LCC.poi' = "green",
      'optA.rep' = "red",
      'optL.rep' = "blue",
      'LCC.rep' = "green",
      'Uni' = "gray"
    )) +
  scale_linetype_manual(
    values = c(
      'optA.poi' = "solid",
      'optL.poi' = "solid",
      'LCC.poi' = "solid",
      'optA.rep' = "dashed",
      'optL.rep' = "dashed",
      'LCC.rep' = "dashed",
      'Uni' = "solid"
    )
  )
# mse_figure
# print(paste0('total time: ',time.cost[3][[1]],'s'))

# plot: estimate of var, for full estimator and for true parameter ############
# var_plot <- data.frame(rou = n.ssp / N,
#                        optA.poi.MLE = MeanVar.cmb[, 1],
#                        optL.poi.MLE = MeanVar.cmb[, 2],
#                        LCC.poi.MLE = MeanVar.cmb[, 3],
#                        optA.poi.true = MeanVar.cmb.true[, 1],
#                        optL.poi.true = MeanVar.cmb.true[, 2],
#                        LCC.poi.true = MeanVar.cmb.true[, 3]
# )
# var_plot <- melt(var_plot, id.vars = 'rou', variable.name="method",
#                value.name="var") # convert data structure to fit ggplot
# var_figure <- ggplot(var_plot, aes(x = rou, y = var)) +
#   geom_line(aes(colour = method, linetype = method)) + geom_point() + ggtitle(var_title) +
#   scale_x_continuous(breaks = n.ssp / N) +
#   theme_set(theme_bw()) +
#   scale_color_manual(
#     values = c(
#       'optA.poi.MLE' = "red",
#       'optL.poi.MLE' = "blue",
#       'LCC.poi.MLE' = "green",
#       'optA.poi.true' = "red",
#       'optL.poi.true' = "blue",
#       'LCC.poi.true' = "green"
#     )) +
#   scale_linetype_manual(
#     values = c(
#       'optA.poi.MLE' = "solid",
#       'optL.poi.MLE' = "solid",
#       'LCC.poi.MLE' = "solid",
#       'optA.poi.true' = "dashed",
#       'optL.poi.true' = "dashed",
#       'LCC.poi.true' = "dashed"
#     )
#   )

# plot: compare of mse and var.true ###########################################
mse_vs_2var_plot <- data.frame(rou = n.ssp / N,
                               optA.poi.MSE = mse.cmb[, 1],
                               optA.poi.Var.MLE = MeanVar.cmb[, 1],
                               optA.poi.Var.true = MeanVar.cmb.true[, 1])
mse_vs_2var_plot<-melt(mse_vs_2var_plot, id.vars = 'rou', variable.name="method",
                       value.name="var") # convert data structure to fit ggplot
mse_vs_2var_figure <- ggplot(mse_vs_2var_plot, aes(x = rou, y = var)) +
  geom_line(aes(colour = method, linetype = method)) + geom_point() + ggtitle(mse_var_title) +
  geom_hline(aes(yintercept = mse.full), linetype = "dashed") +
  scale_x_continuous(breaks = n.ssp / N) +
  theme_set(theme_bw()) +
  scale_color_manual(
    values = c(
      'optA.poi.MSE' = "red",
      'optA.poi.Var.MLE' = "blue",
      'optA.poi.Var.true' = "green"
    )) +
  scale_linetype_manual(
    values = c(
      'optA.poi.MSE' = "solid",
      'optA.poi.Var.MLE' = "dashed",
      'optA.poi.Var.true' = "dashed"
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
