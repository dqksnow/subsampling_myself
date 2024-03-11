library(parallel)
library(ggplot2)
library(reshape2)
source('R/rare_event_functions.R')
source('draft/generate_rare_logit_data.R')
rm(list = setdiff(ls(), lsf.str())) # remove all variables except functions

rare_event_simu <- function(j){ #this function is used for parLapply
  source('R/rare_event_functions.R')
  source('draft/generate_rare_logit_data.R')
  set.seed(j)
  beta.ssp <- beta.cmb <- array(data = NA, dim = c(length(n.ssp), num.method, d))
  SubsampleSize <- var.ssp <- var.cmb <- array(data = NA, dim = c(length(n.ssp), num.method))
  FULL.DATA <- generate.rare.data(N, case = case, beta0)
  X.full <- FULL.DATA$X
  Y.full <- FULL.DATA$Y
  data <- as.data.frame(cbind(Y.full, X.full))
  formula <- Y.full ~ .
  N1 <- sum(Y.full == 1)
  beta0 <- FULL.DATA$beta0
  beta.full <- glm.coef.estimate(X=cbind(1, X.full), Y=Y.full)$beta #beta.full
  alpha <- 0
  for(i in seq_along(n.ssp)){
    n1 <- n.ssp[i]

    Wet.OptA <- rare.logistic.subsampling(formula, data, n.plt, n1, estimate.method = 'Weighted', criterion = 'OptA', alpha = alpha)
    beta.ssp[i, 1, ] <- Wet.OptA$beta.ssp
    beta.cmb[i, 1, ] <- Wet.OptA$beta.cmb
    SubsampleSize[i, 1] <- length(Wet.OptA$index.ssp) - N1
    var.cmb[i, 1] <- sum(diag(Wet.OptA$var.cmb))
    var.ssp[i, 1] <- sum(diag(Wet.OptA$var.ssp))

    Wet.OptL <- rare.logistic.subsampling(formula, data, n.plt, n1, estimate.method = 'Weighted', criterion = 'OptL', alpha = alpha)
    beta.ssp[i, 2, ] <- Wet.OptL$beta.ssp
    beta.cmb[i, 2, ] <- Wet.OptL$beta.cmb
    SubsampleSize[i, 2] <- length(Wet.OptL$index.ssp) - N1
    var.cmb[i, 2] <- sum(diag(Wet.OptL$var.cmb))
    var.ssp[i, 2] <- sum(diag(Wet.OptL$var.ssp))

    Wet.LCC <- rare.logistic.subsampling(formula, data, n.plt, n1, estimate.method = 'Weighted', criterion = 'LCC', alpha = alpha)
    beta.ssp[i, 3, ] <- Wet.LCC$beta.ssp
    beta.cmb[i, 3, ] <- Wet.LCC$beta.cmb
    SubsampleSize[i, 3] <- length(Wet.LCC$index.ssp)
    var.cmb[i, 3] <- sum(diag(Wet.LCC$var.cmb))
    var.ssp[i, 3] <- sum(diag(Wet.LCC$var.ssp))

    Odds.OptA <- rare.logistic.subsampling(formula, data, n.plt, n1, estimate.method = 'LogOddsCorrection', criterion = 'OptA', alpha = alpha)
    beta.ssp[i, 4, ] <- Odds.OptA$beta.ssp
    beta.cmb[i, 4, ] <- Odds.OptA$beta.cmb
    SubsampleSize[i, 4] <- length(Odds.OptA$index.ssp) - N1
    var.cmb[i, 4] <- sum(diag(Odds.OptA$var.cmb))
    var.ssp[i, 4] <- sum(diag(Odds.OptA$var.ssp))

    Odds.OptL <- rare.logistic.subsampling(formula, data, n.plt, n1, estimate.method = 'LogOddsCorrection', criterion = 'OptL', alpha = alpha)
    beta.ssp[i, 5, ] <- Odds.OptL$beta.ssp
    beta.cmb[i, 5, ] <- Odds.OptL$beta.cmb
    SubsampleSize[i, 5] <- length(Odds.OptL$index.ssp) - N1
    var.cmb[i, 5] <- sum(diag(Odds.OptL$var.cmb))
    var.ssp[i, 5] <- sum(diag(Odds.OptL$var.ssp))

    Odds.LCC <- rare.logistic.subsampling(formula, data, n.plt, n1, estimate.method = 'LogOddsCorrection', criterion = 'LCC', alpha = alpha)
    beta.ssp[i, 6, ] <- Odds.LCC$beta.ssp
    beta.cmb[i, 6, ] <- Odds.LCC$beta.cmb
    SubsampleSize[i, 6] <- length(Odds.LCC$index.ssp) - N1 # should I minus N1?
    var.cmb[i, 6] <- sum(diag(Odds.LCC$var.cmb))
    var.ssp[i, 6] <- sum(diag(Odds.LCC$var.ssp))

    Uni <- rare.logistic.subsampling(formula, data, n.plt, n1, estimate.method = 'Uni')
    beta.ssp[i, 7, ] <- beta.cmb[i, 7, ] <- Uni$beta
    var.cmb[i, 7] <-  var.ssp[i, 7] <- sum(diag(Uni$var))
    SubsampleSize[i, 7] <- length(Uni$index) - N1

    # UniW <- rare.logistic.subsampling(formula, data, n.plt, n1, estimate.method = 'UniW')
    # beta.ssp[i, 8, ] <- beta.cmb[i, 8, ] <- UniW$beta
    # var.cmb[i, 8] <- var.ssp[i, 8] <- sum(diag(UniW$var))
    # SubsampleSize[i, 8] <- length(UniW$index) - N1
  }
  return(list(beta.ssp = beta.ssp,
              beta.cmb = beta.cmb,
              beta.full = beta.full,
              var.cmb = var.cmb,
              var.ssp = var.ssp,
              SubsampleSize = SubsampleSize))
}
rpt <- 320
N <-  5 * 1e5
beta0 <- c(NA, -rep(1, 6))
case <- 2
beta0 <- generate.rare.data(N, case = case, beta0)$beta0
d <- length(beta0)
n.plt <- 200
n.ssp <- c(1e3, 2e3, 5e3, 1e4)
# n.ssp_min <- 200
# n.ssp_max <- 1000
# n.ssp <- c(100, seq(n.ssp_min, n.ssp_max, 200))
method.all <- c('Wet_OptA', 'Wet_OptL', 'Wet_LCC', 'Odds_OptA', 'Odds_OptL','Odds_LCC', 'Uni')
num.method <- length(method.all)
Betas.ssp <- Betas.cmb <- array(data = NA, dim = c(rpt, length(n.ssp), num.method, d))
beta.full <- matrix(NA, rpt, d)
SubsampleSize <- var.ssp <- var.cmb <- array(data = NA, dim = c(rpt, length(n.ssp), num.method))

cl <- makeCluster(8) # how many CPU cores are called
t1 <- proc.time()
clusterExport(cl=cl,
              varlist=c('n.plt', 'case', 'n.ssp', 'beta0', 'num.method', 'd', 'N'),
              envir=environment()) # import environment variables into the function 'subsampling_simu'
results <- parLapply(cl, 1:rpt, rare_event_simu)
t2 <- proc.time()
time.cost <- t2 - t1
for(i in 1:rpt){
  Betas.ssp[i,,,] <- results[[i]][[1]]
  Betas.cmb[i,,,] <- results[[i]][[2]]
  beta.full[i,] <- results[[i]][[3]]
  var.cmb[i,,] <- results[[i]][[4]]
  var.ssp[i,,] <- results[[i]][[5]]
  SubsampleSize[i,,] <- results[[i]][[6]]
}

# calculate results and plot
Mean.var.cmb <- apply(var.cmb, c(2, 3), mean)
Mean.var.ssp <- apply(var.ssp, c(2, 3), mean)
# Mean.var.ssp = log(Mean.var.ssp)
# Mean.var.cmb = log(Mean.var.cmb)
MeanSubsampleSize <- apply(SubsampleSize, c(2, 3), mean)
MeanSubsampleSize <- cbind(n.ssp, MeanSubsampleSize)

mse.ssp <- mse.cmb <- array(data = NA, dim = c(length(n.ssp), num.method))
for (idn in seq_along(n.ssp)){
  for (m in 1:num.method) {
    mse.ssp[idn, m] <- mean(apply(sweep(Betas.ssp[, idn, m,], 2, beta0) ^ 2 , 1, sum))
    mse.cmb[idn, m] <- mean(apply(sweep(Betas.cmb[, idn, m,], 2, beta0) ^ 2 , 1, sum))
    # fname <- paste0("rare_event_simu_output/csv/case",case,"N",N,"method",m,"n",n,".csv")
    # write.csv(Betas[,idn,,m],file=fname,row.names=FALSE)
  }
}
# mse.ssp <- log(mse.ssp)
# mse.cmb <- log(mse.cmb)
mse.full <- mean(apply(sweep(beta.full, 2, beta0) ^ 2 , 1, sum))
# mse.full <- log(mse.full)
# plot of mse ########################################################################
mse_cmb_plot <- data.frame(rou = n.ssp/N,
                       Wet_OptA_MSE_cmb = mse.cmb[, 1],
                       Wet_OptL_MSE_cmb = mse.cmb[, 2],
                       Wet_LCC_MSE_cmb = mse.cmb[, 3],
                       Odds_OptA_MSE_cmb = mse.cmb[, 4],
                       Odds_OptL_MSE_cmb = mse.cmb[, 5],
                       Odds_LCC_MSE_cmb = mse.cmb[, 6],
                       Uni_MSE = mse.cmb[, 7]
                       # UniW_MSE = mse.cmb[, 8]
                       )
mse_cmb_plot<-melt(mse_cmb_plot, id.vars = 'rou', variable.name="method",
               value.name="mse") # convert data structure to fit ggplot
mse_cmb_figure <- ggplot(mse_cmb_plot, aes(x = rou, y = mse)) +
  geom_hline(aes(yintercept = mse.full), linetype = "dashed") +
  geom_line(aes(colour = method, linetype=method)) + geom_point() +
  scale_x_continuous(breaks = n.ssp / N) +
  theme_set(theme_bw()) +
  scale_color_manual(
    values = c(
      "Wet_OptA_MSE_cmb" = "red",
      "Wet_OptL_MSE_cmb" = "blue",
      "Wet_LCC_MSE_cmb" = "green",
      "Odds_OptA_MSE_cmb" = "red",
      "Odds_OptL_MSE_cmb" = "blue",
      "Odds_LCC_MSE_cmb" = "green",
      'Uni_MSE' = 'black'
      # 'UniW_MSE' = 'gray'
    )) +
  scale_linetype_manual(
    values = c(
      "Wet_OptA_MSE_cmb" = "solid",
      "Wet_OptL_MSE_cmb" = "solid",
      "Wet_LCC_MSE_cmb" = "solid",
      "Odds_OptA_MSE_cmb" = "dashed",
      "Odds_OptL_MSE_cmb" = "dashed",
      "Odds_LCC_MSE_cmb" = "dashed",
      'Uni_MSE' = 'solid'
      # 'UniW_MSE' = 'solid'
    )
  )

# plot: compare of mse and var ###########################################
mse_var_title <- paste('Rare, Mse vs Var, Case 1, Poisson, Wet, Combined estimator')
mse_vs_var_plot <- data.frame(rou = n.ssp / N,
                              # Wet.OptA.MSE.ssp = mse.ssp[, 1],
                              # Wet.OptA.Var.ssp = Mean.var.ssp[, 1],
                              Wet.OptA.MSE.cmb = mse.cmb[, 4],
                              Wet.OptA.Var.cmb = Mean.var.cmb[, 4])
mse_vs_var_plot<-melt(mse_vs_var_plot, id.vars = 'rou', variable.name="method",
                       value.name="var") # convert data structure to fit ggplot
mse_vs_var_figure <- ggplot(mse_vs_var_plot, aes(x = rou, y = var)) +
  geom_line(aes(colour = method, linetype = method)) + geom_point() + ggtitle(mse_var_title) +
  geom_hline(aes(yintercept = mse.full), linetype = "dashed") +
  scale_x_continuous(breaks = n.ssp / N) +
  theme_set(theme_bw()) +
  scale_color_manual(
    values = c(
      # 'Wet.OptA.MSE.ssp' = "red",
      # 'Wet.OptA.Var.ssp' = "blue",
      'Wet.OptA.MSE.cmb' = "red",
      'Wet.OptA.Var.cmb' = 'blue'
    )) +
  scale_linetype_manual(
    values = c(
      # 'Wet.OptA.MSE.ssp' = "solid",
      # 'Wet.OptA.Var.ssp' = "solid",
      'Wet.OptA.MSE.cmb' = "solid",
      'Wet.OptA.Var.cmb' = 'dashed'
    )
  )
mse_cmb_figure
mse_vs_var_figure
mse.ssp
mse.cmb
MeanSubsampleSize
print(paste0('total time: ',time.cost[3][[1]],'s'))
# ggsave(paste('R/rare_event_simu_output/fig/mse_case_', case, 'n0_', n0, '.png'), width = 6, height = 3, dpi = 300)
# ggsave(paste('draft/0510/rare/case_', case, '_N_', N, '.png'), width = 6, height = 3, dpi = 300)
###############################################################################
# ## single time test
# set.seed(1)
# beta.cmb <- array(data = NA, dim = c(length(n.ssp), num.method, d))
# FULL.DATA <- generate.rare.data(N, case = case, beta0)
# X.full <- FULL.DATA$X
# Y.full <- FULL.DATA$Y
# print(mean(Y.full))
# print(sum(Y.full))
# beta0 <- FULL.DATA$beta0
# for(i in seq_along(n.ssp)){
#   n1 <- n.ssp[i]
#   beta.cmb[i, 1, ] <- logistic.coef.estimate(X=X.full, Y=Y.full)$beta #beta.full
#   beta.cmb[i, 2, ] <- rare.logistic.subsampling(X.full, Y.full, n.plt, n1, method = 'UniW')$beta.est
#   beta.cmb[i, 3, ] <- rare.logistic.subsampling(X.full, Y.full, n.plt, n1, method = 'Uni')$beta.est
#   beta.cmb[i, 4, ] <- rare.logistic.subsampling(X.full, Y.full, n.plt, n1, method = 'Weighted')$beta.est
#   # beta.pilot <- rare.logistic.subsampling(X.full, Y.full, n.plt, n.ssp, method = 'Weighted')$beta.pilot
#   beta.cmb[i, 5, ] <- rare.logistic.subsampling(X.full, Y.full, n.plt, n1, method = 'LogOddsCorrection')$beta.est
#   beta.cmb[i, 6, ] <- rare.logistic.subsampling(X.full, Y.full, n.plt, n1, method = 'LCC')$beta.est
#   # beta.cmb[i, 7, ] <- rare.logistic.subsampling(X.full, Y.full, n.plt, n.ssp, method = 'Weighted', criteria = 'OptA')$beta.est
#   # beta.cmb[i, 8, ] <- rare.logistic.subsampling(X.full, Y.full, n.plt, n.ssp, method = 'LogOddsCorrection', criteria = 'OptA')$beta.est
# }
