library(parallel)
library(ggplot2)
library(reshape2)
rm(list = ls())
source('R/softmax_main_function.R')
source('R/softmax_functions.R')
source('draft/generate_logit_data.R')
title <- paste('MSE, Case 2, Combined Estimator') # WithReplacement
MSPE_title <- paste('MSPE, Case 2, Combined Estimator') # WithReplacement
# var_title <- paste('Var of combined estimator with MLE vs with True, Case 1, SampleMethod WithReplacement')
mse_var_title <- paste('MSE vs Var, Case 2, Combined estimator, N=1e4')
mse_var_ssp_title <- paste('MSE vs Var, Case 2, ssp estimator, N=1e4')

softmax_simu <- function(j) { #this function is used for parLapply
  source('R/softmax_main_function.R')
  source('R/softmax_functions.R')
  source('draft/generate_logit_data.R')
  beta.ssp <- beta.cmb <- array(data = NA, dim = c(length(n.ssp), num.method, (K+1)*d))
  cov.ssp <- cov.cmb <- cov.cmb.full <- cov.ssp.full <- SubsampleSize <-
    array(data = NA, dim = c(length(n.ssp), num.method))
  MSE.cmb <- MSPE.cmb <- MSE.ssp <- array(data = NA, dim = c(length(n.ssp), num.method))
  FULL.DATA <- simu_softmax(seed = j, case = 5, N = N, d = d, K = K, beta = beta.true.sum)
  X <- FULL.DATA$X
  Y <- FULL.DATA$Y
  prob.true <- FULL.DATA$P.true
  fit.full <- nnet::multinom(Y ~ X - 1)
  beta.true.sum <- as.vector(beta.true.sum)
  beta.full <- matrix(G %*% as.vector(t(coef(fit.full))), nrow = d)
  prob.full <- exp(X %*% beta.full)
  prob.full <- prob.full / rowSums(prob.full)
  beta.full <- as.vector(beta.full)
  for(i in seq_along(n.ssp)) {
    MSCLE.OptL.s <- softmax.optimal.subsampling(X, Y, n.plt, n.ssp[i],
                                                criterion = 'OptL',
                                                sampling.method = 'Poisson',
                                                estimate.method = 'MSCLE',
                                                constraint = 'summation')
    beta.ssp[i, 1, ] <- MSCLE.OptL.s$beta.ssp
    cov.ssp[i, 1] <- sum(diag(MSCLE.OptL.s$cov.ssp))
    beta.cmb[i, 1, ] <- MSCLE.OptL.s$beta.cmb
    cov.cmb[i, 1] <- sum(diag(MSCLE.OptL.s$cov.cmb))
    SubsampleSize[i, 1] <- length(MSCLE.OptL.s$index.ssp)
    cov.cmb.full[i, 1] <- sum(diag(MSCLE.OptL.s$cov.cmb.full))
    cov.ssp.full[i, 1] <- sum(diag(MSCLE.OptL.s$cov.ssp.full))
    MSCLE.OptL.s.measures <- fitting.measure(beta.pred = MSCLE.OptL.s$beta.cmb,
                                             beta.true = beta.true.sum,
                                             P.pred = MSCLE.OptL.s$P.cmb,
                                             P.true = prob.true)
    MSE.cmb[i, 1] <- MSCLE.OptL.s.measures$MSE
    MSPE.cmb[i, 1] <- MSCLE.OptL.s.measures$MSPE
    MSCLE.OptL.s.measures <- fitting.measure(beta.pred = MSCLE.OptL.s$beta.ssp,
                                             beta.true = beta.true.sum,
                                             P.pred = MSCLE.OptL.s$P.cmb,
                                             P.true = prob.true)
    MSE.ssp[i, 1] <- MSCLE.OptL.s.measures$MSE
    #################
    MSCLE.OptL.b <- softmax.optimal.subsampling(X, Y, n.plt, n.ssp[i],
                                                criterion = 'OptL',
                                                sampling.method = 'Poisson',
                                                estimate.method = 'MSCLE',
                                                constraint = 'baseline')
    beta.ssp[i, 2, ] <- MSCLE.OptL.b$beta.ssp
    cov.ssp[i, 2] <- sum(diag(MSCLE.OptL.b$cov.ssp))
    beta.cmb[i, 2, ] <- MSCLE.OptL.b$beta.cmb
    cov.cmb[i, 2] <- sum(diag(MSCLE.OptL.b$cov.cmb))
    SubsampleSize[i, 2] <- length(MSCLE.OptL.b$index.ssp)
    cov.cmb.full[i, 2] <- sum(diag(MSCLE.OptL.b$cov.cmb.full))
    cov.ssp.full[i, 2] <- sum(diag(MSCLE.OptL.b$cov.ssp.full))
    MSCLE.OptL.b.measures <- fitting.measure(beta.pred = MSCLE.OptL.b$beta.cmb,
                                             beta.true = beta.true.sum,
                                             P.pred = MSCLE.OptL.b$P.cmb,
                                             P.true = prob.true)
    MSE.cmb[i, 2] <- MSCLE.OptL.b.measures$MSE
    MSPE.cmb[i, 2] <- MSCLE.OptL.b.measures$MSPE
    MSCLE.OptL.b.measures <- fitting.measure(beta.pred = MSCLE.OptL.b$beta.ssp,
                                             beta.true = beta.true.sum,
                                             P.pred = MSCLE.OptL.b$P.cmb,
                                             P.true = prob.true)
    MSE.ssp[i, 2] <- MSCLE.OptL.b.measures$MSE
    ################
    WithRep.OptA.b <- softmax.optimal.subsampling(X, Y, n.plt, n.ssp[i],
                                                  criterion = 'OptL',
                                                  sampling.method = 'Poisson',
                                                  estimate.method = 'Weighted',
                                                  constraint = 'summation')
    beta.ssp[i, 3, ] <- WithRep.OptA.b$beta.ssp
    cov.ssp[i, 3] <- sum(diag(WithRep.OptA.b$cov.ssp))
    beta.cmb[i, 3, ] <- WithRep.OptA.b$beta.cmb
    cov.cmb[i, 3] <- sum(diag(WithRep.OptA.b$cov.cmb))
    SubsampleSize[i, 3] <- length(WithRep.OptA.b$index.ssp)
    cov.cmb.full[i, 3] <- sum(diag(WithRep.OptA.b$cov.cmb.full))
    cov.ssp.full[i, 3] <- sum(diag(WithRep.OptA.b$cov.ssp.full))
    WithRep.OptA.b.measures <- fitting.measure(beta.pred = WithRep.OptA.b$beta.cmb,
                                               beta.true = beta.true.sum,
                                               P.pred = WithRep.OptA.b$P.cmb,
                                               P.true = prob.true)
    MSE.cmb[i, 3] <- WithRep.OptA.b.measures$MSE
    MSPE.cmb[i, 3] <- WithRep.OptA.b.measures$MSPE
    WithRep.OptA.b.measures <- fitting.measure(beta.pred = WithRep.OptA.b$beta.ssp,
                                               beta.true = beta.true.sum,
                                               P.pred = WithRep.OptA.b$P.cmb,
                                               P.true = prob.true)
    MSE.ssp[i, 3] <- WithRep.OptA.b.measures$MSE
    ###################
    WithRep.OptA.s <- softmax.optimal.subsampling(X, Y, n.plt, n.ssp[i],
                                                  criterion = 'OptL',
                                                  sampling.method = 'Poisson',
                                                  estimate.method = 'Weighted',
                                                  constraint = 'baseline')
    beta.ssp[i, 4, ] <- WithRep.OptA.s$beta.ssp
    cov.ssp[i, 4] <- sum(diag(WithRep.OptA.s$cov.ssp))
    beta.cmb[i, 4, ] <- WithRep.OptA.s$beta.cmb
    cov.cmb[i, 4] <- sum(diag(WithRep.OptA.s$cov.cmb))
    SubsampleSize[i, 4] <- length(WithRep.OptA.s$index.ssp)
    cov.cmb.full[i, 4] <- sum(diag(WithRep.OptA.s$cov.cmb.full))
    cov.ssp.full[i, 4] <- sum(diag(WithRep.OptA.s$cov.ssp.full))
    WithRep.OptA.s.measures <- fitting.measure(beta.pred = WithRep.OptA.s$beta.cmb,
                                               beta.true = beta.true.sum,
                                               P.pred = WithRep.OptA.s$P.cmb,
                                               P.true = prob.true)
    MSE.cmb[i, 4] <- WithRep.OptA.s.measures$MSE
    MSPE.cmb[i, 4] <- WithRep.OptA.s.measures$MSPE
    WithRep.OptA.s.measures <- fitting.measure(beta.pred = WithRep.OptA.s$beta.ssp,
                                               beta.true = beta.true.sum,
                                               P.pred = WithRep.OptA.s$P.cmb,
                                               P.true = prob.true)
    MSE.ssp[i, 4] <- WithRep.OptA.s.measures$MSE
    ##################

    Uniform <- softmax.optimal.subsampling(X, Y, n.plt, n.ssp[i],
                                           estimate.method = 'Uniform')
    beta.ssp[i, 5, ] <- beta.cmb[i, 5, ] <- Uniform$beta
    cov.ssp[i, 5] <- cov.cmb[i, 5] <-
      cov.cmb.full[i, 5] <- cov.ssp.full[i, 5] <- sum(diag(Uniform$cov))
    Uniform.measures <- fitting.measure(beta.pred = Uniform$beta,
                                        beta.true = beta.true.sum,
                                        P.pred = Uniform$P,
                                        P.true = prob.true)
    MSE.ssp[i, 5] <- MSE.cmb[i, 5] <- Uniform.measures$MSE
    MSPE.cmb[i, 5] <- Uniform.measures$MSPE
  }
  return(list(beta.ssp = beta.ssp,
              beta.cmb = beta.cmb,
              beta.full = beta.full,
              cov.ssp = cov.ssp,
              cov.cmb = cov.cmb,
              MSE.cmb = MSE.cmb,
              MSPE.cmb = MSPE.cmb,
              SubsampleSize = SubsampleSize,
              cov.cmb.full = cov.cmb.full,
              cov.ssp.full = cov.ssp.full,
              MSE.ssp = MSE.ssp))
}

rpt <- 100
d <- 4 # dim of covariates
K <- 2 # K + 1 classes
G <- rbind(rep(-1/(K+1), K), diag(K) - 1/(K+1)) %x% diag(d)
N <- 1e4
# beta.true <- 0.2 * matrix(-1, d, K)
# beta.true.sum <- cbind(rep(1, d), beta.true)

beta.base <- matrix((0.05 * c(1:8)), d, K)
beta.true <- cbind(0, beta.base)
beta.true.sum <- matrix(G %*% c(beta.base), d, (K+1))

n.plt <- 300
sample.rate <- c(0.03, 0.06, 0.09, 0.12, 0.30)
n.ssp <- sample.rate*N

# n.ssp <- c(500, 1000, 2000, 3000, 5000)
method.all <- c('MSCLE.s', 'MSCLE.b', 'OptL.s', 'OptL.b', 'Uni') #'
num.method <- length(method.all)
Betas.ssp <- array(data = NA, dim = c(rpt, length(n.ssp), num.method, (K+1)*d))
Betas.cmb <- array(data = NA, dim = c(rpt, length(n.ssp), num.method, (K+1)*d))

Beta.full <- matrix(NA, rpt, (K+1)*d)
cov.ssp <- cov.cmb <- cov.cmb.full <- cov.ssp.full <-MSE.cmb <-
  MSE.ssp <- MSPE.cmb <- array(data = NA, dim = c(rpt, length(n.ssp), num.method))
SubsampleSize <- array(data = NA, dim = c(rpt, length(n.ssp), num.method))


cl <- makeCluster(8) # how many CPU cores are called
t1 <- proc.time()
clusterExport(cl=cl,
              varlist=c('n.plt', 'n.ssp', 'beta.true.sum', 'num.method', 'K', 'd', 'N', 'G'),
              envir=environment())
results <- parLapply(cl, 1:rpt, softmax_simu)
t2 <- proc.time()
time.cost <- t2 - t1
##############################################################################
for(i in 1:rpt){
  Betas.ssp[i,,,] <- results[[i]][[1]]
  Betas.cmb[i,,,] <- results[[i]][[2]]
  Beta.full[i,] <- results[[i]][[3]]
  cov.ssp[i,,] <- results[[i]][[4]]
  cov.cmb[i,,] <- results[[i]][[5]]
  MSE.cmb[i,,] <- results[[i]][[6]]
  MSPE.cmb[i,,] <- results[[i]][[7]]
  SubsampleSize[i,,] <- results[[i]][[8]]
  cov.cmb.full[i,,] <- results[[i]][[9]]
  cov.ssp.full[i,,] <- results[[i]][[10]]
  MSE.ssp[i,,] <- results[[i]][[11]]
}
#
MeanMSE.ssp <- apply(MSE.ssp, c(2, 3), mean)
MeanMSE.cmb <- apply(MSE.cmb, c(2, 3), mean)
MeanMSPE.cmb <- apply(MSPE.cmb, c(2, 3), mean)
Meancov.ssp <- apply(cov.ssp, c(2, 3), mean)
Meancov.ssp.full <- apply(cov.ssp.full, c(2, 3), mean)
Meancov.cmb <- apply(cov.cmb, c(2, 3), mean)
Meancov.cmb.full <- apply(cov.cmb.full, c(2, 3), mean)
MeanSubsampleSize <- apply(SubsampleSize, c(2, 3), mean)
MeanSubsampleSize <- cbind(n.ssp, MeanSubsampleSize)
##############################################################################
# mse.ssp <- mse.cmb <- array(data = NA, dim = c(length(n.ssp), num.method))
# for (idn in seq_along(n.ssp)){
#   for (m in 1:num.method) {
#     mse.ssp[idn, m] <- mean(apply((Betas.ssp[, idn, m,] - Beta.full) ^ 2 , 1, sum))
#     mse.cmb[idn, m] <- mean(apply((Betas.cmb[, idn, m,] - Beta.full) ^ 2 , 1, sum))
#     # fname <- paste0("rare_event_simu_output/csv/case",case,"N",N,"method",m,"n",n,".csv")
#     # write.csv(Betas[,idn,,m],file=fname,row.names=FALSE)
#   }
# }

mse.full <- mean(rowSums((Beta.full - matrix(replicate(rpt, as.vector(beta.true.sum)), nrow = rpt, byrow = T)) ^ 2))

# mse <- log(mse)
# mse.full <- log(mse.full)

# plot: mse ###################################################################
mse_plot <- data.frame(rou = n.ssp / N,
                       MSCLE.s = MeanMSE.cmb[, 1],
                       MSCLE.b = MeanMSE.cmb[, 2],
                       OptL.s = MeanMSE.cmb[, 3],
                       OptL.b = MeanMSE.cmb[, 4],
                       # OptL.b = MeanMSE.cmb[, 4],
                       # OptL.s = MeanMSE.cmb[, 5],
                       Uni = MeanMSE.cmb[, 5]
                       # MSPE.WithRep = MeanMSE.cmb[, 1],
                       # MSPE.Poi = MeanMSE.cmb[, 2],
                       # OptA.s.WithRep = MeanMSE.cmb[, 3],
                       # OptA.s.Poi = MeanMSE.cmb[, 4],
                       # Uni = MeanMSE.cmb[, 6]
)
mse_plot<-melt(mse_plot, id.vars = 'rou', variable.name = "method",
               value.name="mse") # convert data structure to fit ggplot
mse_figure <- ggplot(mse_plot, aes(x = rou, y = mse)) +
  # geom_hline(aes(yintercept = mse.full), linetype = "dashed") +
  geom_line(aes(colour = method, linetype = method)) + geom_point() + ggtitle(title) +
  scale_x_continuous(breaks = n.ssp / N) +
  theme_set(theme_bw()) +
  scale_color_manual(
    values = c(
      'MSCLE.s' = "green",
      'MSCLE.b' = "green",
      # 'OptA.b' = "blue",
      'OptL.s' = "blue",
      # 'OptL.b' = "red",
      'OptL.b' = "red",
      'Uni' = "gray"
      # 'MSPE.WithRep' = "green",
      # 'MSPE.Poi' = "green",
      # 'OptA.s.WithRep' = "blue",
      # 'OptA.s.Poi' = "blue",
      # 'Uni' = "gray"
    )) +
  scale_linetype_manual(
    values = c(
      'MSCLE.s' = "solid",
      'MSCLE.b' = "dashed",

      # 'OptA.b' = "dashed",
      'OptL.s' = "solid",
      # 'OptL.b' = "dashed",
      'OptL.b' = "solid",
      'Uni' = "dashed"
      # 'MSPE.WithRep' = "solid",
      # 'MSPE.Poi' = "dashed",
      # 'OptA.s.WithRep' = "solid",
      # 'OptA.s.Poi' = "dashed",
      # 'Uni' = "dashed"
    )
  )
mse_figure
print(paste0('total time: ',time.cost[3][[1]],'s'))
# plot: mspe ###################################################################
mspe_plot <- data.frame(rou = n.ssp / N,
                        MSCLE.s = MeanMSPE.cmb[, 1],
                        MSCLE.b = MeanMSPE.cmb[, 2],
                        OptL.s = MeanMSPE.cmb[, 3],
                        OptL.b = MeanMSPE.cmb[, 4],
                        # OptL.b = MeanMSE.cmb[, 4],
                        # OptL.s = MeanMSE.cmb[, 5],
                        Uni = MeanMSPE.cmb[, 5]
                        # MSPE.WithRep = MeanMSPE.cmb[, 1],
                        # MSPE.Poi = MeanMSPE.cmb[, 2],
                        # OptA.s.WithRep = MeanMSPE.cmb[, 3],
                        # OptA.s.Poi = MeanMSPE.cmb[, 4],
                        # Uni = MeanMSPE.cmb[, 6]
)
mspe_plot<-melt(mspe_plot, id.vars = 'rou', variable.name = "method",
                value.name="mspe") # convert data structure to fit ggplot
mspe_figure <- ggplot(mspe_plot, aes(x = rou, y = mspe)) +
  # geom_hline(aes(yintercept = mse.full), linetype = "dashed") +
  geom_line(aes(colour = method, linetype = method)) + geom_point() + ggtitle(MSPE_title) +
  scale_x_continuous(breaks = n.ssp / N) +
  theme_set(theme_bw()) +
  scale_color_manual(
    values = c(
      'MSCLE.s' = "green",
      'MSCLE.b' = "green",
      'OptL.s' = "blue",
      'OptL.b' = "red",
      'Uni' = "gray"
    )) +
  scale_linetype_manual(
    values = c(
      'MSCLE.s' = "solid",
      'MSCLE.b' = "dashed",
      'OptL.s' = "solid",
      'OptL.b' = "solid",
      'Uni' = "dashed"
    )
  )
mspe_figure
# print(paste0('total time: ',time.cost[3][[1]],'s'))
print(1 - Meancov.cmb.full / MeanMSE.cmb)

print(1 - Meancov.cmb.full / Meancov.cmb)

print(1 - Meancov.cmb / MeanMSE.cmb)

# plot: compare of empirical mse and theoretical var ###########################################
mse_vs_var_plot <- data.frame(rou = n.ssp / N,
                              MSCLE.MSE = MeanMSE.ssp[, 1],
                              MSCLE.Var.true = Meancov.ssp[, 1],
                              # OptA.b.MSE = MeanMSE.ssp[, 2],
                              OptA.s.MSE = MeanMSE.cmb[, 3],
                              OptA.s.Var.true = Meancov.cmb[, 3],
                              OptA.s.Var.full = Meancov.cmb.full[, 3]
                              # OptA.b.Var.true = Meancov.cmb[, 2],
                              # MSPE.Var.full = Meancov.cmb.full[, 1],
                              # OptA.b.Var.full = Meancov.cmb.full[, 2],
                              # OptA.s.Var.full = Meancov.cmb.full[, 3]
)
mse_vs_var_plot<-melt(mse_vs_var_plot, id.vars = 'rou', variable.name="method",
                      value.name="var") # convert data structure to fit ggplot
mse_vs_var_figure <- ggplot(mse_vs_var_plot, aes(x = rou, y = var)) +
  geom_line(aes(colour = method, linetype = method)) + geom_point() + ggtitle(mse_var_title) +
  # geom_hline(aes(yintercept = mse.full), linetype = "dashed") +
  scale_x_continuous(breaks = n.ssp / N) +
  theme_set(theme_bw()) +
  scale_color_manual(
    values = c(
      'MSCLE.MSE' = "green",
      'MSCLE.Var.true' = "green",
      'OptA.s.MSE' = "blue",
      'OptA.s.Var.true' = "blue",
      'OptA.s.Var.full' = "blue"
    )) +
  scale_linetype_manual(
    values = c(
      'MSCLE.MSE' = "solid",
      'MSCLE.Var.true' = "dashed",
      'OptA.s.MSE' = "solid",
      'OptA.s.Var.true' = "dashed",
      'OptA.s.Var.full' = "dotted"
    )
  )

print(paste0('total time: ',time.cost[3][[1]],'s'))
# mse
# mse_figure
# var_figure
mse_vs_var_figure
# MeanSubsampleSize
# ggsave(paste('draft/0510/', title, '.png'), width = 6, height = 3, dpi = 300)


# plot: compare of empirical mse and theoretical var,
# cmb estimator and full estimator ###########################################
mse_vs_var_plot <- data.frame(rou = n.ssp / N,
                              MSCLE.MSE = MeanMSE.ssp[, 1],
                              MSCLE.Var.true = Meancov.ssp[, 1],
                              # MSPE.Var.full = Meancov.ssp.full[, 1],
                              # OptA.b.MSE = MeanMSE.ssp[, 2],
                              OptA.s.MSE = MeanMSE.ssp[, 3],
                              OptA.s.Var.true = Meancov.ssp[, 3]
                              # OptA.s.Var.full = Meancov.ssp.full[, 3]
                              # OptA.b.Var.true = Meancov.ssp[, 2],
                              # MSPE.Var.full = Meancov.ssp.full[, 1],
                              # OptA.b.Var.full = Meancov.ssp.full[, 2],
                              # OptA.s.Var.full = Meancov.ssp.full[, 3]
)
mse_vs_var_plot<-melt(mse_vs_var_plot, id.vars = 'rou', variable.name="method",
                      value.name="var") # convert data structure to fit ggplot
mse_vs_var_figure <- ggplot(mse_vs_var_plot, aes(x = rou, y = var)) +
  geom_line(aes(colour = method, linetype = method)) + geom_point() + ggtitle(mse_var_ssp_title) +
  # geom_hline(aes(yintercept = mse.full), linetype = "dashed") +
  scale_x_continuous(breaks = n.ssp / N) +
  theme_set(theme_bw()) +
  scale_color_manual(
    values = c(
      'MSCLE.MSE' = "green",
      'MSCLE.Var.true' = "green",
      'OptA.s.MSE' = "blue",
      'OptA.s.Var.true' = "blue"
    )) +
  scale_linetype_manual(
    values = c(
      'MSCLE.MSE' = "solid",
      'MSCLE.Var.true' = "dashed",
      'OptA.s.MSE' = "solid",
      'OptA.s.Var.true' = "dashed"
    )
  )

print(paste0('total time: ',time.cost[3][[1]],'s'))
# mse
# mse_figure
# var_figure
mse_vs_var_figure
