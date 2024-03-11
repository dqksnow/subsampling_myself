library(parallel)
library(ggplot2)
library(reshape2)
rm(list = ls())
source('R/softmax_main_function.R')
source('R/softmax_functions.R')
source('draft/generate_logit_data.R')
title <- paste('MSE, Case 1, Combined Estimator') # WithReplacement
MSPE_title <- paste('MSPE, Case 1, Combined Estimator') # WithReplacement
# var_title <- paste('Var of combined estimator with MLE vs with True, Case 1, SampleMethod WithReplacement')
mse_var_title <- paste('MSE vs Var, Case 1, Combined estimator, N=1e4')
mse_var_ssp_title <- paste('MSE vs Var, Case 1, ssp estimator, N=1e4')

softmax_simu <- function(j) { #this function is used for parLapply
  source('R/softmax_main_function.R')
  source('R/softmax_functions.R')
  source('draft/generate_logit_data.R')
  beta.ssp <- beta.cmb <- array(data = NA, dim = c(length(n.ssp), num.method, (K+1)*d))
  cov.ssp <- cov.cmb <- cov.cmb.full <- cov.ssp.full <- SubsampleSize <-
    array(data = NA, dim = c(length(n.ssp), num.method))
  MSE.cmb <- MSPE.cmb <- MSE.ssp <- array(data = NA, dim = c(length(n.ssp), num.method))
  FULL.DATA <- simu_softmax(seed = j, case = 1, N = N, d = d, K = K, beta = beta.true.sum)
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
    WithRep.MSPE <- softmax.optimal.subsampling(X, Y, n.ssp[i], n.plt,
                                                criterion = 'MSPE',
                                                sampling.method = 'Poisson',
                                                estimate.method = 'Weighted')
    beta.ssp[i, 1, ] <- WithRep.MSPE$beta.ssp
    cov.ssp[i, 1] <- sum(diag(WithRep.MSPE$cov.plt))
    beta.cmb[i, 1, ] <- WithRep.MSPE$beta.cmb
    cov.cmb[i, 1] <- sum(diag(WithRep.MSPE$cov.cmb))
    cov.cmb.full[i, 1] <- sum(diag(WithRep.MSPE$cov.cmb.full))
    cov.ssp.full[i, 1] <- sum(diag(WithRep.MSPE$cov.ssp.full))
    SubsampleSize[i, 1] <- length(WithRep.MSPE$index.ssp)
    WithRep.MSPE.measures <- fitting.measure(beta.pred = WithRep.MSPE$beta.cmb,
                                             beta.true = beta.true.sum,
                                             P.pred = WithRep.MSPE$P.cmb,
                                             P.true = prob.true)
    MSE.cmb[i, 1] <- WithRep.MSPE.measures$MSE
    MSPE.cmb[i, 1] <- WithRep.MSPE.measures$MSPE
    WithRep.MSPE.measures <- fitting.measure(beta.pred = WithRep.MSPE$beta.plt,
                                             beta.true = beta.true.sum,
                                             P.pred = WithRep.MSPE$P.cmb,
                                             P.true = prob.true)
    MSE.ssp[i, 1] <- WithRep.MSPE.measures$MSE
    #################
    WithRep.OptA.b <- softmax.optimal.subsampling(X, Y, n.ssp[i], n.plt,
                                                  criterion = 'OptA',
                                                  sampling.method = 'Poisson',
                                                  estimate.method = 'Weighted',
                                                  constraint = 'baseline')
    beta.ssp[i, 2, ] <- WithRep.OptA.b$beta.ssp
    cov.ssp[i, 2] <- sum(diag(WithRep.OptA.b$cov.plt))
    beta.cmb[i, 2, ] <- WithRep.OptA.b$beta.cmb
    cov.cmb[i, 2] <- sum(diag(WithRep.OptA.b$cov.cmb))
    cov.cmb.full[i, 2] <- sum(diag(WithRep.OptA.b$cov.cmb.full))
    cov.ssp.full[i, 2] <- sum(diag(WithRep.OptA.b$cov.ssp.full))
    SubsampleSize[i, 2] <- length(WithRep.OptA.b$index.ssp)
    WithRep.OptA.b.measures <- fitting.measure(beta.pred = WithRep.OptA.b$beta.cmb,
                                               beta.true = beta.true.sum,
                                               P.pred = WithRep.OptA.b$P.cmb,
                                               P.true = prob.true)
    MSE.cmb[i, 2] <- WithRep.OptA.b.measures$MSE
    MSPE.cmb[i, 2] <- WithRep.OptA.b.measures$MSPE
    WithRep.OptA.b.measures <- fitting.measure(beta.pred = WithRep.OptA.b$beta.plt,
                                               beta.true = beta.true.sum,
                                               P.pred = WithRep.OptA.b$P.cmb,
                                               P.true = prob.true)
    MSE.ssp[i, 2] <- WithRep.OptA.b.measures$MSE
    ###################
    WithRep.OptA.s <- softmax.optimal.subsampling(X, Y, n.ssp[i], n.plt,
                                                  criterion = 'OptA',
                                                  sampling.method = 'Poisson',
                                                  estimate.method = 'Weighted',
                                                  constraint = 'summation')
    beta.ssp[i, 3, ] <- WithRep.OptA.s$beta.ssp
    cov.ssp[i, 3] <- sum(diag(WithRep.OptA.s$cov.plt))
    beta.cmb[i, 3, ] <- WithRep.OptA.s$beta.cmb
    cov.cmb[i, 3] <- sum(diag(WithRep.OptA.s$cov.cmb))
    cov.cmb.full[i, 3] <- sum(diag(WithRep.OptA.s$cov.cmb.full))
    cov.ssp.full[i, 3] <- sum(diag(WithRep.OptA.s$cov.ssp.full))
    SubsampleSize[i, 3] <- length(WithRep.OptA.s$index.ssp)
    WithRep.OptA.s.measures <- fitting.measure(beta.pred = WithRep.OptA.s$beta.cmb,
                                               beta.true = beta.true.sum,
                                               P.pred = WithRep.OptA.s$P.cmb,
                                               P.true = prob.true)
    MSE.cmb[i, 3] <- WithRep.OptA.s.measures$MSE
    MSPE.cmb[i, 3] <- WithRep.OptA.s.measures$MSPE
    WithRep.OptA.s.measures <- fitting.measure(beta.pred = WithRep.OptA.s$beta.plt,
                                               beta.true = beta.true.sum,
                                               P.pred = WithRep.OptA.s$P.cmb,
                                               P.true = prob.true)
    MSE.ssp[i, 3] <- WithRep.OptA.s.measures$MSE
    ####################
    # WithRep.OptL.b <- softmax.optimal.subsampling(X, Y, n.plt, n.ssp[i],
    #                                                      criterion = 'OptL',
    #                                                      sampling.method = 'WithReplacement',
    #                                                      estimate.method = 'Weighted',
    #                                                      constraint = 'baseline')
    # beta.ssp[i, 4, ] <- WithRep.OptL.b$beta.ssp
    # cov.ssp[i, 4] <- sum(diag(WithRep.OptL.b$cov.ssp))
    # beta.cmb[i, 4, ] <- WithRep.OptL.b$beta.cmb
    # cov.cmb[i, 4] <- sum(diag(WithRep.OptL.b$cov.cmb))
    # cov.cmb.full[i, 4] <- sum(diag(WithRep.OptL.b$cov.cmb.full))
    # cov.ssp.full[i, 4] <- sum(diag(WithRep.OptL.b$cov.ssp.full))
    # SubsampleSize[i, 4] <- length(WithRep.OptL.b$index.ssp)
    # WithRep.OptL.b.measures <- fitting.measure(beta.pred = WithRep.OptL.b$beta.cmb,
    #                                            beta.true = beta.true.sum,
    #                                            P.pred = WithRep.OptL.b$P.cmb,
    #                                            P.true = prob.true)
    # MSE.cmb[i, 4] <- WithRep.OptL.b.measures$MSE
    # MSPE.cmb[i, 4] <- WithRep.OptL.b.measures$MSPE
    # WithRep.OptL.b.measures <- fitting.measure(beta.pred = WithRep.OptL.b$beta.ssp,
    #                                            beta.true = beta.true.sum,
    #                                            P.pred = WithRep.OptL.b$P.cmb,
    #                                            P.true = prob.true)
    # MSE.ssp[i, 4] <- WithRep.OptL.b.measures$MSE
    # ##################
    # WithRep.OptL.s <- softmax.optimal.subsampling(X, Y, n.plt, n.ssp[i],
    #                                                       criterion = 'OptL',
    #                                                       sampling.method = 'WithReplacement',
    #                                                       estimate.method = 'Weighted',
    #                                                       constraint = 'summation')
    # beta.ssp[i, 5, ] <- WithRep.OptL.s$beta.ssp
    # cov.ssp[i, 5] <- sum(diag(WithRep.OptL.s$cov.ssp))
    # beta.cmb[i, 5, ] <- WithRep.OptL.s$beta.cmb
    # cov.cmb[i, 5] <- sum(diag(WithRep.OptL.s$cov.cmb))
    # cov.cmb.full[i, 5] <- sum(diag(WithRep.OptL.s$cov.cmb.full))
    # cov.ssp.full[i, 5] <- sum(diag(WithRep.OptL.s$cov.ssp.full))
    # SubsampleSize[i, 5] <- length(WithRep.OptL.s$index.ssp)
    # WithRep.OptL.s.measures <- fitting.measure(beta.pred = WithRep.OptL.s$beta.cmb,
    #                                            beta.true = beta.true.sum,
    #                                            P.pred = WithRep.OptL.s$P.cmb,
    #                                            P.true = prob.true)
    # MSE.cmb[i, 5] <- WithRep.OptL.s.measures$MSE
    # MSPE.cmb[i, 5] <- WithRep.OptL.s.measures$MSPE
    # WithRep.OptL.s.measures <- fitting.measure(beta.pred = WithRep.OptL.s$beta.ssp,
    #                                            beta.true = beta.true.sum,
    #                                            P.pred = WithRep.OptL.s$P.cmb,
    #                                            P.true = prob.true)
    # MSE.ssp[i, 5] <- WithRep.OptL.s.measures$MSE
    ####################
    Uniform <- softmax.optimal.subsampling(X, Y, n.ssp[i], n.plt,
                                           estimate.method = 'Uniform')
    beta.ssp[i, 4, ] <- beta.cmb[i, 4, ] <- Uniform$beta
    cov.ssp[i, 4] <- cov.cmb[i, 4] <-
      cov.cmb.full[i, 4] <- cov.ssp.full[i, 4] <- sum(diag(Uniform$cov))
    Uniform.measures <- fitting.measure(beta.pred = Uniform$beta,
                                        beta.true = beta.true.sum,
                                        P.pred = Uniform$P,
                                        P.true = prob.true)
    MSE.cmb[i, 4] <- Uniform.measures$MSE
    MSPE.cmb[i, 4] <- Uniform.measures$MSPE
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
d <- 3 # dim of covariates
K <- 5 # K + 1 classes
G <- rbind(rep(-1/(K+1), K), diag(K) - 1/(K+1)) %x% diag(d)
N <- 1e4
beta.true <- 0.2 * matrix(-1, d, K)
# beta.true <- matrix(-c(0.8,0.8,0.2,0.2), d, K)
beta.true.sum <- cbind(rep(1, d), beta.true)
n.plt <- 300
sample.rate <- c(0.03, 0.06, 0.09, 0.12, 0.30)
n.ssp <- sample.rate*N

# n.ssp <- c(500, 1000, 2000, 3000, 5000)
method.all <- c('MSPE', 'OptA.b', 'OptA.s', 'Uni') #'OptL.b', 'OptL.s',
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
                       MSPE = MeanMSE.cmb[, 1],
                       OptA.b = MeanMSE.cmb[, 2],
                       OptA.s = MeanMSE.cmb[, 3],
                       # OptL.b = MeanMSE.cmb[, 4],
                       # OptL.s = MeanMSE.cmb[, 5],
                       Uni = MeanMSE.cmb[, 4]
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
      'MSPE' = "green",
      'OptA.b' = "blue",
      'OptA.s' = "blue",
      # 'OptL.b' = "red",
      # 'OptL.s' = "red",
      'Uni' = "gray"
      # 'MSPE.WithRep' = "green",
      # 'MSPE.Poi' = "green",
      # 'OptA.s.WithRep' = "blue",
      # 'OptA.s.Poi' = "blue",
      # 'Uni' = "gray"
    )) +
  scale_linetype_manual(
    values = c(
      'MSPE' = "solid",
      'OptA.b' = "dashed",
      'OptA.s' = "solid",
      # 'OptL.b' = "dashed",
      # 'OptL.s' = "solid",
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
                        MSPE = MeanMSPE.cmb[, 1],
                        OptA.b = MeanMSPE.cmb[, 2],
                        OptA.s = MeanMSPE.cmb[, 3],
                        # OptL.b = MeanMSPE.cmb[, 4],
                        # OptL.s = MeanMSPE.cmb[, 5],
                        Uni = MeanMSPE.cmb[, 4]
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
      'MSPE' = "green",
      'OptA.b' = "blue",
      'OptA.s' = "blue",
      # 'OptL.b' = "red",
      # 'OptL.s' = "red",
      'Uni' = "gray"
      # 'MSPE.WithRep' = "green",
      # 'MSPE.Poi' = "green",
      # 'OptA.s.WithRep' = "blue",
      # 'OptA.s.Poi' = "blue",
      # 'Uni' = "gray"
    )) +
  scale_linetype_manual(
    values = c(
      'MSPE' = "solid",
      'OptA.b' = "dashed",
      'OptA.s' = "solid",
      # 'OptL.b' = "dashed",
      # 'OptL.s' = "solid",
      'Uni' = "dashed"
      # 'MSPE.WithRep' = "solid",
      # 'MSPE.Poi' = "dashed",
      # 'OptA.s.WithRep' = "solid",
      # 'OptA.s.Poi' = "dashed",
      # 'Uni' = "dashed"
    )
  )
mspe_figure
# print(paste0('total time: ',time.cost[3][[1]],'s'))
print(1 - Meancov.cmb.full / MeanMSE.cmb)

print(1 - Meancov.cmb.full / Meancov.cmb)

print(1 - Meancov.cmb / MeanMSE.cmb)

# plot: compare of empirical mse and theoretical var ###########################################
mse_vs_var_plot <- data.frame(rou = n.ssp / N,
                              MSPE.MSE = MeanMSE.cmb[, 1],
                              MSPE.Var.true = Meancov.cmb[, 1],
                              MSPE.Var.full = Meancov.cmb.full[, 1],
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
      'MSPE.MSE' = "green",
      'OptA.s.MSE' = "blue",
      'MSPE.Var.true' = "green",
      'OptA.s.Var.true' = "blue",
      'MSPE.Var.full' = "green",
      # 'OptA.b.Var.full' = "red"
      'OptA.s.Var.full' = "blue"
    )) +
  scale_linetype_manual(
    values = c(
      'MSPE.MSE' = "solid",
      # 'OptA.b.MSE' = "solid",
      'OptA.s.MSE' = "solid",
      'MSPE.Var.true' = "dashed",
      # 'OptA.b.Var.true' = "dashed",
      'OptA.s.Var.true' = "dashed",
      'MSPE.Var.full' = "dotted",
      # 'OptA.b.Var.full' = "dotted",
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
                              MSPE.MSE = MeanMSE.cmb[, 1],
                              MSPE.Var.true = Meancov.cmb[, 1],
                              MSPE.Var.full = Meancov.cmb.full[, 1],
                              # OptA.b.MSE = MeanMSE.cmb[, 2],
                              OptA.s.MSE = MeanMSE.cmb[, 3],
                              OptA.s.Var.true = Meancov.cmb[, 3],
                              OptA.s.Var.full = Meancov.ssp.full[, 3]
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
      'MSPE.MSE' = "green",
      'OptA.s.MSE' = "blue",
      'MSPE.Var.true' = "green",
      'OptA.s.Var.true' = "blue",
      'MSPE.Var.full' = "green",
      # 'OptA.b.Var.full' = "red"
      'OptA.s.Var.full' = "blue"
    )) +
  scale_linetype_manual(
    values = c(
      'MSPE.MSE' = "solid",
      # 'OptA.b.MSE' = "solid",
      'OptA.s.MSE' = "solid",
      'MSPE.Var.true' = "dashed",
      # 'OptA.b.Var.true' = "dashed",
      'OptA.s.Var.true' = "dashed",
      'MSPE.Var.full' = "dotted",
      # 'OptA.b.Var.full' = "dotted",
      'OptA.s.Var.full' = "dotted"
    )
  )

print(paste0('total time: ',time.cost[3][[1]],'s'))
# mse
# mse_figure
# var_figure
mse_vs_var_figure
