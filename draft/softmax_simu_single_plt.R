rm(list = ls())
library(ggplot2)
library(reshape2)
source('R/softmax_main_function.R')
source('R/softmax_functions.R')
source('draft/generate_logit_data.R')
rpt <- 100
N <-  1e4
d <- 4 # dim of covariates
K <- 2 # K + 1 classes
G <- rbind(rep(-1/(K+1), K), diag(K) - 1/(K+1)) %x% diag(d)
N <- 1e4
# beta.true <- 0.25 * matrix(-1, d, K)
# beta.true.sum <- cbind(rep(1, d), beta.true)

beta.true <- matrix((0.05 * c(1:8)), d, K)
beta.true.sum <- matrix(G %*% c(beta.true), d, (K+1))

n.ssp <- 1000
n.plt <- c(500, 1000, 1500, 3000)
method.all <- c('MSPE', 'OptA.s', 'Uni')
num.method <- length(method.all)
Betas.plt <- Betas.ssp <- Betas.cmb <-
  array(data = NA, dim = c(rpt, length(n.plt), num.method, (K)*d))
cov.plt <- cov.ssp <- cov.cmb <- cov.cmb.full <- cov.ssp.full <- SubsampleSize <-
  array(data = NA, dim = c(rpt, length(n.plt), num.method))
Beta.full <- matrix(NA, rpt, (K)*d)
cov.full <- cov.full.nnet <- rep(NA, rpt)
###############################################################################
t1 <- proc.time()
for (j in 1:rpt) {
  print(j)
  FULL.DATA <- simu_softmax(seed = j, case = 5, N = N, d = d, K = K, beta = beta.true.sum)
  X <- FULL.DATA$X
  Y <- FULL.DATA$Y
  Y.matrix <- matrix(0, nrow = N, ncol = K)
  Y.matrix[cbind(seq_along(Y), Y)] <- 1
  fit.full <- nnet::multinom(Y ~ X - 1, trace = FALSE)
  P1.full <- fit.full$fitted.values
  ddL.plt <- softmax.ddL(X = X, P = P1.full[, -1], p = 1, K, d, scale = 1)
  dL.sq.plt <- softmax.dL.sq(X = X, Y.matrix = Y.matrix,
                             P = P1.full[, -1], p = 1, K = K, d = d, scale = 1)
  cov.full.nnet[j] <- sum((summary(fit.full)$standard.errors)^2)
  cov.full[j] <- sum(diag(solve(ddL.plt) %*% dL.sq.plt %*% solve(ddL.plt)))

  beta.full.b <- as.vector(t(coef(fit.full)))
  beta.full.s <- G %*% beta.full.b
  Beta.full[j, ] <- beta.full.b



  for(i in seq_along(n.plt)) {
    MSPE <- softmax.optimal.subsampling(X, Y, n.plt[i], n.ssp,
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
    OptA <- softmax.optimal.subsampling(X, Y, n.plt[i], n.ssp,
                                        criterion = 'OptA',
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
    Uniform <- softmax.optimal.subsampling(X, Y, n.plt[i], n.ssp,
                                           estimate.method = 'Uniform')
    Betas.plt[j, i, 3, ] <- Betas.ssp[j, i, 3, ] <- Betas.cmb[j, i, 3, ] <-
      Uniform$beta
    cov.plt[j, i, 3] <- cov.ssp[j, i, 3] <- cov.cmb[j, i, 3] <- sum(diag(Uniform$cov))
  }
}
t2 <- proc.time()
time.cost <- t2 - t1
##############################################################################
mse.plt <- mse.ssp <- mse.cmb <- array(data = NA, dim = c(length(n.plt), num.method))
beta.true.matrix <- matrix(replicate(rpt, as.vector(beta.true)), nrow = rpt, byrow = T)
for (idn in seq_along(n.plt)){
  for (m in 1:num.method) {
    mse.plt[idn, m] <- mean(rowSums((Betas.plt[, idn, m, ] - beta.true.matrix) ^ 2))
    mse.ssp[idn, m] <- mean(rowSums((Betas.ssp[, idn, m, ] - beta.true.matrix) ^ 2))
    mse.cmb[idn, m] <- mean(rowSums((Betas.cmb[, idn, m, ] - beta.true.matrix) ^ 2))
  }
}
mse.full <- mean(rowSums((Beta.full - beta.true.matrix) ^ 2))
Meancov.full <- mean(cov.full)
Meancov.full.nnet <- mean(cov.full.nnet)
Meancov.plt <- apply(cov.plt, c(2, 3), mean)
Meancov.ssp <- apply(cov.ssp, c(2, 3), mean)
Meancov.cmb <- apply(cov.cmb, c(2, 3), mean)
boxplot(cbind(rowSums((Beta.full - beta.true.matrix) ^ 2), cov.full, cov.full.nnet))
boxplot(cbind(as.vector(rowSums((Betas.plt[, idn, m, ] - beta.true.matrix) ^ 2)), as.vector(cov.plt[, idn, m])))
mse.plt / Meancov.plt
##############################################################################
mse_vs_var_plot <- data.frame(rou = n.plt / N,
                              MSPE.MSE = mse.plt[, 1],
                              OptA.MSE = mse.plt[, 2],
                              Uni.MSE = mse.plt[, 3],
                              MSPE.Var.true = Meancov.plt[, 1],
                              OptA.Var.true = Meancov.plt[, 2],
                              Uni.Var.true = Meancov.plt[, 3]
)
mse_vs_var_plot<-melt(mse_vs_var_plot, id.vars = 'rou', variable.name="method",
                      value.name="var") # convert data structure to fit ggplot
mse_vs_var_figure <- ggplot(mse_vs_var_plot, aes(x = rou, y = var)) +
  geom_line(aes(colour = method, linetype = method)) + geom_point() +
  geom_hline(aes(yintercept = mse.full), linetype = "dashed") +
  scale_x_continuous(breaks = n.ssp / N) +
  theme_set(theme_bw()) +
  scale_color_manual(
    values = c(
      MSPE.MSE = "green",
      OptA.MSE = "blue",
      Uni.MSE = 'gray',
      MSPE.Var.true = "green",
      OptA.Var.true = "blue",
      Uni.Var.true = 'gray'
    )) +
  scale_linetype_manual(
    values = c(
      MSPE.MSE = "solid",
      OptA.MSE = "solid",
      Uni.MSE = "solid",
      MSPE.Var.true = "dashed",
      OptA.Var.true = "dashed",
      Uni.Var.true = "dashed"
    )
  )
# mse
# mse_figure
# var_figure
mse_vs_var_figure
print(paste0('total time: ',time.cost[3][[1]],'s'))
a <- sum((summary(fit.full)$standard.errors)^2)


