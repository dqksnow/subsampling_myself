source('R/rare_event_functions.R')
# source('draft/rare_event_functions_old.R')
source('draft/generate_rare_logit_data.R')
# source('draft/rare_event_functions_old.R')
rm(list = setdiff(ls(), lsf.str())) # remove all variables except functions
N <-  5 * 1e5
beta0 <- c(NA, -rep(1, 6))
d <- length(beta0)
case <- 1
set.seed(2)
FULL.DATA <- generate.rare.data(N, case=case, beta0)
X.full <- FULL.DATA$X
Y.full <- FULL.DATA$Y
data <- as.data.frame(cbind(Y.full, X.full))
# colnames(data) <- c("response", "intercept", paste0("v", 1:6))
# formula <- response ~ . -1
formula <- Y.full ~ .
beta0 <- FULL.DATA$beta0
print(paste('mean: ', mean(Y.full), 'sum: ', sum(Y.full)))
n.plt <- 1000
n.ssp <- 1000
num.method <- 6
Betas <- matrix(NA, d, num.method)
## nss.star <- array(data = NA, dim = c(lns, rpt, 3))
# FULL <- logistic.coef.estimate(X = X.full, Y = Y.full)$beta
# MSE_FULL <- sum((FULL - beta0)^2)
LogOddsCorrection.OptA.results <- rare.logistic.subsampling(formula,
                                                            data,
                                                            n.plt,
                                                            n.ssp,
                                                            estimate.method = 'LogOddsCorrection',
                                                            criterion = 'OptA')
# subsampling.summary(LogOddsCorrection.LCC)
Weighted.LCC.results <- rare.logistic.subsampling(formula,
                                                  data,
                                                  n.plt,
                                                  n.ssp,
                                                  estimate.method = 'Weighted',
                                                  criterion = 'OptL')
# subsampling.summary(Weighted.LCC.results)

Weighted.OptA.results <- rare.logistic.subsampling(formula,
                                                   data,
                                                   n.plt,
                                                   n.ssp,
                                                   estimate.method = 'Weighted',
                                                   criterion = 'OptA')

subsampling.summary(Weighted.OptA.results)



Weighted.OptL.results <- rare.logistic.subsampling(formula,
                                                   data,
                                                   n.plt,
                                                   n.ssp,
                                                   estimate.method = 'Weighted',
                                                   criterion = 'OptL')


subsampling.summary(LogOddsCorrection.OptA.results)
LogOddsCorrection.OptL.results <- rare.logistic.subsampling(formula, data, n.plt, n.ssp, estimate.method = 'LogOddsCorrection', criterion = 'OptL')



Uni <- rare.logistic.subsampling(formula, data, n0, n.ssp, estimate.method = 'Uni')
subsampling.summary(Uni)
UniW <- rare.logistic.subsampling(formula, data, n0, n.ssp, estimate.method = 'UniW')
subsampling.summary(UniW)
# Betas[,1] <- logistic.coef.estimate(X=X.full, Y=Y.full)$beta #beta.full
# Betas[,2] <- rare.logistic.subsampling(X.full, Y.full, n0, nss, method = 'UniW')$beta.ssp
# Betas[,3] <- rare.logistic.subsampling(X.full, Y.full, n0, nss, method = 'Uni')$beta.ssp
# Betas[,4] <- rare.logistic.subsampling(X.full, Y.full, n0, nss, criterion = 'LCC')$beta.ssp

# Betas[,1] <- Weighted.OptL.results$beta.ssp
# Betas[,2] <- Weighted.OptA.results$beta.ssp
# Betas[,3] <- Weighted.LCC.results$beta.ssp
# Betas[,4] <- Weighted.OptL.results$beta.cmb
# Betas[,5] <- Weighted.OptA.results$beta.cmb
# Betas[,6] <- Weighted.LCC.results$beta.cmb

Betas[,1] <- LogOddsCorrection.OptL.results$beta.cmb
Betas[,2] <- LogOddsCorrection.OptA.results$beta.cmb
Betas[,3] <- LogOddsCorrection.LCC$beta.cmb
Betas[,4] <- LogOddsCorrection.OptL.results$beta.cmb
Betas[,5] <- LogOddsCorrection.OptA.results$beta.cmb
Betas[,6] <- LogOddsCorrection.LCC$beta.cmb


MSE <- as.data.frame(apply((Betas - beta0)^2 ,2,sum))
MSE
# colnames(MSE) <- 'MSE'
# rownames(MSE) <- c('full', 'UniW', 'Uni', 'LCC', 'Weighted+OptL', 'Weighted+OptA', 'LOC+OptL', 'LOC+OptA',
#                    'Weighted+OptL+cmb', 'Weighted+OptA+cmb', 'LOC+OptL+cmb', 'LOC+OptA+cmb')
# BETAS <- as.data.frame(cbind(beta0, Betas), row.names =c(1:d))
# colnames(BETAS) <- c('true', 'full', 'UniW', 'Uni', 'LCC', 'Weighted+OptL', 'Weighted+OptA', 'LOC+OptL', 'LOC+OptA',
#                      'Weighted+OptL+cmb', 'Weighted+OptA+cmb', 'LOC+OptL+cmb', 'LOC+OptA+cmb')
# rownames(BETAS) <- c(1:d)
MSE
# log(MSE)
# BETAS

### test
# N <- nrow(X.full)
# d <- ncol(X.full)
# N1 <- sum(Y.full)
# N0 <- N - N1
# loc0 = Y.full == 0
# # pilot step
# plt.estimate.results <- pilot.estimate(X.full, Y.full, n0)
# beta.plt <- plt.estimate.results$beta.plt
# beta.plt <- plt.estimate.results$beta.plt
# MN.plt <- plt.estimate.results$MN.plt
# Psi.plt <- plt.estimate.results$Psi.plt
# P.plt <- plt.estimate.results$P.plt
# index.plt <- plt.estimate.results$index.plt
# cc.plt <- plt.estimate.results$cc.plt
#
# criterion <- 'OptA'
# method <- "Weighted"
# # subsampling step
# ssp.results <- subsampling(X.full, Y.full, nss, criterion,
#                                    P.plt, MN.plt, index.plt, cc.plt)
# index.ssp <- ssp.results$index.ssp
# w.ssp <- ssp.results$w.ssp
# x.ssp <- X.full[index.ssp, ]
# y.ssp <- Y.full[index.ssp]
# # subsample estimate step
# ssp.estimate.reults <- subsample.estimate(x.ssp, y.ssp, w.ssp, beta.plt, criterion, method)
# beta.ssp <- ssp.estimate.reults$beta.ssp
# MN.ssp <- ssp.estimate.reults$MN.ssp
# Psi.ssp <- ssp.estimate.reults$Psi.ssp
# #
# if (criterion == 'LCC'){
#   beta.ssp <-
#     logistic.coef.estimate(X = x.ssp, Y = y.ssp, weights = w.ssp)$beta + beta.plt
# } else if (method == "Weighted"){
#   beta.ssp <-
#     logistic.coef.estimate(X = x.ssp, Y = y.ssp, weights = 1 / w.ssp)$beta
# } else if (method == "LogOddsCorrection"){
#   p.ssp <- subsampling.results$p.ssp
#   beta.ssp <-
#     logistic.coef.estimate(
#       X = x.ssp,
#       Y = y.ssp, #start balue
#       weights = rep(1, length(y.ssp)),
#       offset = -log(p.ssp)
#     )$beta
#   # beta.ssp <- getMSLE(x,y,p,beta.plt)$beta
# }
# P.ssp  <- pbeta(x.ssp, beta.ssp)
# MN.ssp <- MN(x.ssp, P.ssp, 1) #1/w.ssp
# Psi.ssp <- Psi(x.ssp, y.ssp, P.ssp, 1)
