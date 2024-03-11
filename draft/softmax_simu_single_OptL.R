rm(list = ls())
library(ggplot2)
library(reshape2)
source('R/softmax_main_function.R')
source('R/softmax_functions.R')
source('draft/generate_logit_data.R')
rpt <- 1
N <-  1e4
d <- 3 # dim of covariates
K <- 5 # K + 1 classes
G <- rbind(rep(-1/(K+1), K), diag(K) - 1/(K+1)) %x% diag(d)
beta.true <- 0.2 * matrix(-1, d, K)
beta.true.sum <- cbind(rep(1, d), beta.true)

# beta.true <- matrix((0.05 * c(1:8)), d, K)
# beta.true.sum <- matrix(G %*% c(beta.true), d, (K+1))

n.plt <- 1000
n.ssp <- c(500, 1500, 3000)
method.all <- c('Uni')
num.method <- length(method.all)
Betas.cmb <- Betas.ssp <- Betas.cmb <-
  array(data = NA, dim = c(rpt, length(n.ssp), num.method, (K+1)*d))
cov.plt <- cov.plt.nnet <- unique.plt <- SubsampleSize <-
  array(data = NA, dim = c(rpt, length(n.ssp), num.method))
Beta.full <- matrix(NA, rpt, (K+1)*d)
cov.full <- cov.full.nnet <- rep(NA, rpt)
###############################################################################
t1 <- proc.time()
for (j in 1:rpt) {
  print(j)
  FULL.DATA <- simu_softmax(seed = j, case = 1, N = N, d = d, K = K, beta = beta.true.sum)
  X <- FULL.DATA$X
  Y <- FULL.DATA$Y
  Y.matrix <- matrix(0, nrow = N, ncol = K)
  Y.matrix[cbind(seq_along(Y), Y)] <- 1
  # fit.full <- nnet::multinom(Y ~ X - 1, trace = FALSE)
  # P1.full <- fit.full$fitted.values
  # ddL.plt <- softmax.ddL(X = X, P = P1.full[, -1], p = 1, K, d, scale = 1)
  # dL.sq.plt <- softmax.dL.sq(X = X, Y.matrix = Y.matrix,
  #                            P = P1.full[, -1], p = 1, K = K, d = d, scale = 1)
  # cov.full.nnet[j] <- sum((summary(fit.full)$standard.errors)^2)
  # cov.full[j] <- sum(diag(solve(ddL.plt) %*% dL.sq.plt %*% solve(ddL.plt)))
  #
  # beta.full.b <- as.vector(t(coef(fit.full)))
  # Beta.full[j, ] <- beta.full.b

  for(i in seq_along(n.ssp)) {
    OptA <- softmax.optimal.subsampling(X, Y, n.plt, n.ssp[i],
                                        criterion = 'OptL',
                                        sampling.method = 'WithReplacement',
                                        estimate.method = 'Weighted',
                                        constraint = 'summation')
    Betas.cmb[j, i, 1, ]  <- OptA$beta.cmb
    cov.plt[j, i, 1] <-  sum(diag(OptA$cov.cmb))
    cov.plt.nnet[j, i, 1] <- sum(diag(OptA$cov.cmb.full))
  }
}
t2 <- proc.time()
time.cost <- t2 - t1
##############################################################################
mse.plt <- array(data = NA, dim = c(length(n.ssp), num.method))
beta.true.matrix <- matrix(replicate(rpt, as.vector(beta.true.sum)), nrow = rpt, byrow = T)
for (idn in seq_along(n.ssp)){
  for (m in 1:num.method) {
    mse.plt[idn, m] <- mean(rowSums((Betas.cmb[, idn, m, ] - beta.true.matrix) ^ 2))
  }
}
mse.full <- rowSums((Beta.full - beta.true.matrix) ^ 2)
Meanmse.full <- mean(rowSums((Beta.full - beta.true.matrix) ^ 2))
Meancov.full <- mean(cov.full)
Meancov.full.nnet <- mean(cov.full.nnet)
Meancov.plt <- apply(cov.plt, c(2, 3), mean)
Meancov.plt.nnet <- apply(cov.plt.nnet, c(2, 3), mean)
Meanunique.plt <- apply(unique.plt, c(2, 3), mean)
# boxplot(cbind(mse.full, cov.full, cov.full.nnet))
# boxplot(cbind(rowSums((Betas.cmb[, idn, m, ] - beta.true.matrix) ^ 2),
#               cov.plt[, idn, 1],
#               cov.plt.nnet[, idn, 1]))
mse.plt / Meancov.plt
Meancov.plt / Meancov.plt.nnet
mse.plt / Meancov.plt.nnet
##############################################################################
title <- paste('N=', N, ', fix n.ssp=', n.plt, ' , rou=(n.plt+n.ssp)/N, ', 'Pilot is sampled with rep') # WithReplacement
plot_subtitle <- paste('OptL: Rep + Weighted + cmb sample')
mse_vs_var_plot <- data.frame(rou = (n.plt+n.ssp) / N,
                              OptL.MSE = mse.plt[, 1],
                              OptL.Var.true = Meancov.plt[, 1])
                              # OptL.Var.full = Meancov.plt.nnet[, 1])
mse_vs_var_plot<-melt(mse_vs_var_plot, id.vars = 'rou', variable.name="method",
                      value.name="var") # convert data structure to fit ggplot
mse_vs_var_figure <- ggplot(mse_vs_var_plot, aes(x = rou, y = var)) +
  geom_line(aes(colour = method, linetype = method)) + geom_point() +
  # geom_hline(aes(yintercept = Meanmse.full), linetype = "dashed") +
  scale_x_continuous(breaks = (n.plt+n.ssp) / N) +
  ggtitle(title) +
  labs(subtitle = plot_subtitle) +
  ylab("mse") +
  theme_set(theme_bw()) +
  scale_color_manual(
    values = c(
      OptL.MSE = "green",
      OptL.Var.true = "green"
      # OptL.Var.full = 'gray'
    )) +
  scale_linetype_manual(
    values = c(
      OptL.MSE = "solid",
      OptL.Var.true = "dashed"
      # OptL.Var.full = "dotted"
    )
  )
# mse
# mse_figure
# var_figure
mse_vs_var_figure
print(paste0('total time: ',time.cost[3][[1]],'s'))


