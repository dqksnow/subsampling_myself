rm(list = ls())
library(ggplot2)
library(reshape2)
source('R/glm_subsampling_functions.R')
source('R/family_expand.R')
source('draft/generate_logit_data.R')
rpt <- 100
N <-  1e4
beta.true <- c(rep(0.5, 7))
d <- length(beta.true)


n.plt <- 500
n.ssp <- c(500, 1000, 1500, 3000)
method.all <- c('Uni', 'OptA')
num.method <- length(method.all)
Betas.plt <- Betas.ssp <- Betas.cmb <-
  array(data = NA, dim = c(rpt, length(n.ssp), num.method, d))
cov.plt <- cov.plt.nnet <- unique.plt <- SubsampleSize <-
  array(data = NA, dim = c(rpt, length(n.ssp), num.method))
Beta.full <- matrix(NA, rpt, d)
cov.full <- cov.full.nnet <- rep(NA, rpt)
###############################################################################
t1 <- proc.time()
for (j in 1:rpt) {
  print(j)
  family <- binomial.expand()
  # FULL.DATA <- simu_mzNormal(seed = j, N = N, beta0 = beta.true, corr = 0.5)
  FULL.DATA <- simu_nzNormal(seed = j, N = N, beta0 = beta.true, corr = 0.5)
  X.full <- FULL.DATA$X
  Y.full <- FULL.DATA$Y
  data <- as.data.frame(cbind(Y.full, X.full))
  formula <- Y.full ~ .
  beta.full <- glm.coef.estimate(X=cbind(1, X.full), Y=Y.full, family = family)$beta #beta.full
  Beta.full[j, ] <- beta.full
  for(i in seq_along(n.ssp)) {
    Uniform <- glm.optimal.subsampling(formula, data, n.plt, n.ssp[i],
                                   estimate.method = 'Uni',
                                   family = family$family.name)
    Betas.plt[j, i, 1, ]  <- Uniform$beta
    cov.plt[j, i, 1] <-  sum(diag(Uniform$var))
    unique.plt[j, i, 1] <- (length(unique(Uniform$index)))
    OptA <- glm.optimal.subsampling(formula, data, n.plt, n.ssp[i],
                                       estimate.method = 'Weighted',
                                    sampling.method = 'WithReplacement',
                                    criterion = 'OptA',
                                       family = family$family.name)
    Betas.plt[j, i, 2, ]  <- OptA$beta.cmb
    cov.plt[j, i, 2] <-  sum(diag(OptA$var.cmb.true))
    unique.plt[j, i, 2] <- (length(unique(OptA$index.ssp)))
  }
}
t2 <- proc.time()
time.cost <- t2 - t1
##############################################################################
mse.plt <- array(data = NA, dim = c(length(n.ssp), num.method))
beta.true.matrix <- matrix(replicate(rpt, as.vector(beta.true)), nrow = rpt, byrow = T)
# for (idn in seq_along(n.plt)){
#   for (m in 1:num.method) {
#     mse.plt[idn, m] <- mean(rowSums((Betas.plt[, idn, m, ] - beta.true.matrix) ^ 2))
#   }
# }
for (idn in seq_along(n.ssp)){
  for (m in 1:num.method) {
    mse.plt[idn, m] <- mean(rowSums((Betas.plt[, idn, m, ] - Beta.full) ^ 2))
  }
}
mse.full <- rowSums((Beta.full - beta.true.matrix) ^ 2)
Meanmse.full <- mean(rowSums((Beta.full - beta.true.matrix) ^ 2))
Meancov.full <- mean(cov.full)
Meancov.plt <- apply(cov.plt, c(2, 3), mean)
Meanunique.plt <- apply(unique.plt, c(2, 3), mean)
# boxplot(cbind(rowSums((Betas.plt[, idn, m, ] - beta.true.matrix) ^ 2),
#               cov.plt[, idn, 1]))
# (mse.plt * (Meanunique.plt / (n.plt+n.ssp))) / Meancov.plt
mse.plt / Meancov.plt
##############################################################################
mse_vs_var_plot <- data.frame(rou = (n.plt+n.ssp) / N,
                              Uni.MSE = mse.plt[, 1],
                              Uni.Var = Meancov.plt[, 1])
mse_vs_var_plot<-melt(mse_vs_var_plot, id.vars = 'rou', variable.name="method",
                      value.name="var") # convert data structure to fit ggplot
mse_vs_var_figure <- ggplot(mse_vs_var_plot, aes(x = rou, y = var)) +
  geom_line(aes(colour = method, linetype = method)) + geom_point() +
  geom_hline(aes(yintercept = Meanmse.full), linetype = "dashed") +
  scale_x_continuous(breaks = (n.plt+n.ssp) / N) +
  theme_set(theme_bw()) +
  scale_color_manual(
    values = c(
      Uni.MSE = "green",
      Uni.Var = "green"
    )) +
  scale_linetype_manual(
    values = c(
      Uni.MSE = "solid",
      Uni.Var = "dotted"
    )
  )
# mse
# mse_figure
# var_figure
mse_vs_var_figure
print(paste0('total time: ',time.cost[3][[1]],'s'))


