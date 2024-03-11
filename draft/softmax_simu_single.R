rm(list = ls())
library(ggplot2)
library(reshape2)
source('R/softmax_main_function.R')
source('R/softmax_functions.R')
source('draft/generate_logit_data.R')

rpt <- 500
N <-  5e4
d <- 2 # dim of covariates
K <- 2 # K + 1 classes
G <- rbind(rep(-1/(K+1), K), diag(K) - 1/(K+1)) %x% diag(d)
beta.true <- 0.5 * matrix(-1, d, K)
beta.true.sum <- cbind(rep(1, d), beta.true)
# beta.true <- matrix((0.05 * c(1:8)), d, K)
# beta.true.sum <- matrix(G %*% c(beta.true), d, (K+1))
n.plt <- 3000
n.ssp <- c(500, 1500, 3000)
method.all <- c('MSPE', 'OptA.s', 'Uni')
num.method <- length(method.all)
Betas.plt <- Betas.ssp <- Betas.cmb <-
  array(data = NA, dim = c(rpt, length(n.ssp), num.method, (K+1)*d))
cov.plt <- cov.ssp <- cov.cmb <- cov.cmb.full <- cov.ssp.full <- SubsampleSize <-
  array(data = NA, dim = c(rpt, length(n.ssp), num.method))
Beta.full <- matrix(NA, rpt, (K+1)*d)
###############################################################################
t1 <- proc.time()
for (j in 1:rpt) {
  print(j)
  FULL.DATA <- simu_softmax(seed = j, case = 1, N = N, d = d, K = K, beta = beta.true.sum)
  X <- FULL.DATA$X
  Y <- FULL.DATA$Y
  fit.full <- nnet::multinom(Y ~ X - 1, trace = FALSE)
  beta.full.b <- as.vector(t(coef(fit.full)))
  beta.full.s <- G %*% beta.full.b
  Beta.full[j, ] <- beta.full.s
  for(i in seq_along(n.ssp)) {
    MSPE <- softmax.optimal.subsampling(X, Y, n.plt, n.ssp[i],
                                                criterion = 'MSPE',
                                                sampling.method = 'WithReplacement',
                                                estimate.method = 'Weighted',
                                        constraint = 'summation')
    Betas.plt[j, i, 1, ] <- MSPE$beta.plt
    Betas.ssp[j, i, 1, ] <- MSPE$beta.ssp
    Betas.cmb[j, i, 1, ] <- MSPE$beta.cmb
    cov.plt[j, i, 1] <- sum(diag(MSPE$cov.plt))
    cov.ssp[j, i, 1] <- sum(diag(MSPE$cov.ssp))
    cov.cmb[j, i, 1] <- sum(diag(MSPE$cov.cmb))
    #
    OptA <- softmax.optimal.subsampling(X, Y, n.plt, n.ssp[i],
                                                  criterion = 'OptL',
                                                  sampling.method = 'WithReplacement',
                                                  estimate.method = 'Weighted',
                                                  constraint = 'summation')
    Betas.plt[j, i, 2, ] <- OptA$beta.plt
    Betas.ssp[j, i, 2, ] <- OptA$beta.ssp
    Betas.cmb[j, i, 2, ] <- OptA$beta.cmb
    cov.plt[j, i, 2] <- sum(diag(OptA$cov.plt))
    cov.ssp[j, i, 2] <- sum(diag(OptA$cov.ssp))
    cov.cmb[j, i, 2] <- sum(diag(OptA$cov.cmb))
    #
    Uniform <- softmax.optimal.subsampling(X, Y, n.plt, n.ssp[i],
                                           estimate.method = 'Uniform',
                                           sampling.method = 'WithReplacement'
                                           )
    Betas.plt[j, i, 3, ] <- Betas.ssp[j, i, 3, ] <- Betas.cmb[j, i, 3, ] <-
      Uniform$beta
    cov.plt[j, i, 3] <- cov.ssp[j, i, 3] <- cov.cmb[j, i, 3] <- sum(diag(Uniform$cov))
  }
}
t2 <- proc.time()
time.cost <- t2 - t1
##############################################################################
mse.plt <- mse.ssp <- mse.cmb <- array(data = NA, dim = c(length(n.ssp), num.method))
beta.true.sum.matrix <- matrix(replicate(rpt, as.vector(beta.true.sum)), nrow = rpt, byrow = T)
for (idn in seq_along(n.ssp)){
  for (m in 1:num.method) {
    mse.plt[idn, m] <- mean(rowSums((Betas.plt[, idn, m, ] - beta.true.sum.matrix) ^ 2))
    mse.ssp[idn, m] <- mean(rowSums((Betas.ssp[, idn, m, ] - beta.true.sum.matrix) ^ 2))
    mse.cmb[idn, m] <- mean(rowSums((Betas.cmb[, idn, m, ] - beta.true.sum.matrix) ^ 2))
  }
}
mse.full <- mean(rowSums((Beta.full - beta.true.sum.matrix) ^ 2))
Meancov.plt <- apply(cov.plt, c(2, 3), mean)
Meancov.ssp <- apply(cov.ssp, c(2, 3), mean)
Meancov.cmb <- apply(cov.cmb, c(2, 3), mean)
mse.plt/Meancov.plt
mse.ssp/Meancov.ssp
mse.cmb/Meancov.cmb

##############################################################################
title <- paste('N=', N, ', n.plt=', n.plt, ', rou=n.ssp/N, ', 'Pilot is sampled with rep') # WithReplacement
plot_subtitle <- paste('MSPE: WithRep + Weighted + cmb, OptL: WithRep + Weighted + cmb, ',
                       'Uniform: WithRep + Weighted')
mse_vs_var_plot <- data.frame(rou = n.ssp / N,
                              MSPE.MSE = mse.cmb[, 1],
                              OptL.MSE = mse.cmb[, 2],
                              Uni.MSE = mse.cmb[, 3],
                              MSPE.Var.true = Meancov.cmb[, 1],
                              OptL.Var.true = Meancov.cmb[, 2],
                              Uni.Var.true = Meancov.cmb[, 3]
)
mse_vs_var_plot<-melt(mse_vs_var_plot, id.vars = 'rou', variable.name="method",
                      value.name="var") # convert data structure to fit ggplot
mse_vs_var_figure <- ggplot(mse_vs_var_plot, aes(x = rou, y = var)) +
  geom_line(aes(colour = method, linetype = method)) + geom_point() +
  ggtitle(title) +
  labs(subtitle = plot_subtitle) +
  ylab("mse") +
  # geom_hline(aes(yintercept = mse.full), linetype = "dashed") +
  scale_x_continuous(breaks = n.ssp / N) +
  theme_set(theme_bw()) +
  scale_color_manual(
    values = c(
      MSPE.MSE = "green",
      OptL.MSE = "blue",
      Uni.MSE = 'gray',
      MSPE.Var.true = "green",
      OptL.Var.true = "blue",
      Uni.Var.true = 'gray'
    )) +
  scale_linetype_manual(
    values = c(
      MSPE.MSE = "solid",
      OptL.MSE = "solid",
      Uni.MSE = "solid",
      MSPE.Var.true = "dashed",
      OptL.Var.true = "dashed",
      Uni.Var.true = "dashed"
    )
  )
# mse
# mse_figure
# var_figure
mse_vs_var_figure
print(paste0('total time: ',time.cost[3][[1]],'s'))
