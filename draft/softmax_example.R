source('R/softmax_functions.R')
source('R/softmax_main_function.R')
source('draft/generate_logit_data.R')
# library(MASS)
# library(expm)
rm(list = setdiff(ls(), lsf.str()))

d <- 2 # dim of covariates
K <- 5 # K + 1 classes
G <- rbind(rep(-1/(K+1), K), diag(K) - 1/(K+1)) %x% diag(d)
N <- 1e4
n.plt <- 1000
n.ssp <- 1000
# beta.true <- matrix((0.05 * c(1:8)), d, K)
# beta.true.sum <- matrix(G %*% c(beta.true), d, (K+1))
beta.true <- 0.2 * matrix(-1, d, K)
beta.true.sum <- cbind(rep(1, d), beta.true)

# beta.true <- matrix(c(rep(0,d), rep(-1, K*d)), d, K+1) # baseline constraint
FULL.DATA <- simu_softmax(seed = 2, case = 1, N = N, d = d, K = K, beta = beta.true.sum)
X <- FULL.DATA$X
Y <- FULL.DATA$Y
table(Y)
prob.true <- FULL.DATA$P.true


# Y.matrix <- matrix(0, nrow = N, ncol = K)
# Y.matrix[cbind(seq_along(Y), Y)] <- 1
# fit.full <- nnet::multinom(Y ~ X - 1)
# summary(fit.full)
# # fit.full <- nnet::multinom(Y.matrix ~ X - 1)
# prob.fit.full <- fit.full$fitted.values[, -1]
# prob.full <- pbeta.multi(X, t(coef(fit.full)))
# prob.full <- cbind(1 - rowSums(prob.full), prob.full)
# beta.full <- matrix(G %*% as.vector(t(coef(fit.full))), nrow = d)
#
# beta.full <- G %*% as.vector(t(coef(fit.full)))
# mse.full <- sum((beta.full - beta.true.sum)^2)
# mspe.full <- sum(rowSums(prob.full - prob.Kplus2)^2)

# beta.true <- as.vector(beta.true.sum)

#####
#####
Poi.optL.b <- softmax.optimal.subsampling(X, Y, n.plt, n.ssp,
                                          criterion = 'OptL',
                                          sampling.method = 'WithReplacement',
                                          estimate.method = 'Weighted',
                                          constraint = 'summation')
sum((Poi.optL.b$beta.ssp - as.vector(beta.true.sum))^2)
sum((Poi.optL.b$beta.cmb - as.vector(beta.true.sum))^2)
sum(diag(Poi.optL.b$cov.ssp))
sum(diag(Poi.optL.b$cov.cmb))

WithRep.MSPE.measures <- fitting.measure(beta.pred = Poi.optL.b$beta.cmb,
                                         beta.true = beta.true,
                                         P.pred = Poi.optL.b$P.cmb,
                                         P.true = prob.true)
WithRep.MSPE.measures$MSE
WithRep.MSPE.measures$MSPE
#####
WithRep.MSPE <- softmax.optimal.subsampling(X, Y, n.plt, n.ssp,
                                            criterion = 'MSPE',
                                            sampling.method = 'Poisson',
                                            estimate.method = 'Weighted',
                                            constraint = 'baseline')
sum((WithRep.MSPE$beta.ssp - beta.true.sum)^2)
sum((WithRep.MSPE$beta.cmb - beta.true.sum)^2)
sum(diag(WithRep.MSPE$cov.ssp))
sum(diag(WithRep.MSPE$cov.cmb))
WithRep.MSPE$beta.cmb
beta.plt <- matrix(WithRep.MSPE$beta.plt, nrow = d)
MSE <- sum((beta.pred - beta.true)^2)
prob.plt <- exp(X %*% beta.plt)
prob.plt <- prob.plt / rowSums(prob.plt)
#####
Poi.MSPE <- softmax.optimal.subsampling(X, Y, n.plt, n.ssp,
                                            criterion = 'MSPE',
                                            sampling.method = 'Poisson',
                                            estimate.method = 'Weighted',
                                            constraint = 'summation')
beta.plt <- Poi.MSPE$beta.plt
Poi.MSPE.measures <- fitting.measure(beta.pred = Poi.MSPE$beta.plt,
                                         beta.true = beta.true,
                                         P.pred = Poi.MSPE$P.cmb,
                                         P.true = prob.full)
Poi.MSPE.measures$MSE
sum(diag(Poi.MSPE$var.plt))






measure.WithRep.MSPE <- fitting.measure(beta.pred = beta.plt,
                                        beta.true = beta.true.sum,
                                        P.pred = prob.plt,
                                        P.true = prob.Kplus2)
#
n.plt <- 900
n.ssp <- 0
Uniform <- softmax.optimal.subsampling(X, Y, n.plt, n.ssp,
                                       estimate.method = 'Uniform')
beta.uni <- matrix(Uniform$beta, nrow = d)
prob.uni <- exp(X %*% beta.uni)
prob.uni <- prob.uni / rowSums(prob.uni)
prob.uni.test <- exp(X %*% beta.uni - matrixStats::rowLogSumExps(X %*% beta.uni))
measure.Uniform <- fitting.measure(beta.pred = beta.uni,
                           beta.true = beta.true.sum,
                           P.pred = prob.uni,
                           P.true = prob.Kplus2)


###############################################################################
var.WithRep.MSPE <- sum(diag(WithRep.MSPE$var.cmb))
WithRep.MSPE.mse.plt <- sum((WithRep.MSPE$beta.plt - beta.true)^2)
WithRep.MSPE.mse.ssp <- sum((WithRep.MSPE$beta.ssp - beta.true)^2)
WithRep.MSPE.mse.cmb <- sum((WithRep.MSPE$beta.cmb - beta.true)^2)
#####
Uni <- softmax.optimal.subsampling(X, Y, n.plt, n.ssp, estimate.method = 'Uniform')
Uni.mse <- sum((Uni$beta - beta.true)^2)
#####
WithRep.optA.b <- softmax.optimal.subsampling(X, Y, n.plt, n.ssp,
                                            criterion = 'OptA',
                                            sampling.method = 'WithReplacement',
                                            estimate.method = 'Weighted',
                                            constraint = 'baseline')
WithRep.optA.b.mse.plt <- sum((WithRep.optA.b$beta.plt - beta.true.sum)^2)
WithRep.optA.b.mse.ssp <- sum((WithRep.optA.b$beta.ssp - beta.true.sum)^2)
WithRep.optA.b.mse.cmb <- sum((WithRep.optA.b$beta.cmb - beta.true.sum)^2)
#####
WithRep.optL.b <- softmax.optimal.subsampling(X, Y, n.plt, n.ssp,
                                              criterion = 'OptL',
                                              sampling.method = 'WithReplacement',
                                              estimate.method = 'Weighted',
                                              constraint = 'baseline')
WithRep.optL.b.mse.plt <- sum((WithRep.optL.b$beta.plt - beta.true.sum)^2)
WithRep.optL.b.mse.ssp <- sum((WithRep.optL.b$beta.ssp - beta.true.sum)^2)
WithRep.optL.b.mse.cmb <- sum((WithRep.optL.b$beta.cmb - beta.true.sum)^2)
#####
WithRep.optA.s <- softmax.optimal.subsampling(X, Y, n.plt, n.ssp,
                                              criterion = 'OptA',
                                              sampling.method = 'WithReplacement',
                                              estimate.method = 'Weighted',
                                              constraint = 'summation')
WithRep.optA.s.mse.plt <- sum((WithRep.optA.s$beta.plt - beta.true.sum)^2)
WithRep.optA.s.mse.ssp <- sum((WithRep.optA.s$beta.ssp - beta.true.sum)^2)
WithRep.optA.s.mse.cmb <- sum((WithRep.optA.s$beta.cmb - beta.true.sum)^2)
#####
WithRep.optL.s <- softmax.optimal.subsampling(X, Y, n.plt, n.ssp,
                                              criterion = 'OptL',
                                              sampling.method = 'WithReplacement',
                                              estimate.method = 'Weighted',
                                              constraint = 'summation')
WithRep.optL.s.mse.plt <- sum((WithRep.optL.s$beta.plt - beta.true.sum)^2)
WithRep.optL.s.mse.ssp <- sum((WithRep.optL.s$beta.ssp - beta.true.sum)^2)
WithRep.optL.s.mse.cmb <- sum((WithRep.optL.s$beta.cmb - beta.true.sum)^2)




prob.plt <- exp( X %*% WithRep_optA_b$beta.plt)
prob.plt <- prob.plt / rowSums(prob.plt)
mspe.plt <- sum(rowSums(prob.plt - prob.full)^2) # why so small
