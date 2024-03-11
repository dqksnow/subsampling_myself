
# source('R/poisson_subsampling_functions.R')
source('R/family_expand.R')
source('R/glm_subsampling_functions.R')
source('draft/generate_nb_data.R')
rm(list = setdiff(ls(), lsf.str()))
N <-  1e4
beta0 <- rep(0.5, 7)
# beta0 <- c(0.5, rep(0.5, 7))

seed <- 1
set.seed(seed)
case <- 1
family <- negative.binomial.expand()
FULL.DATA <- generate.nb.data(N, case=case, beta0, v = 2, seed = seed)
X.full <- FULL.DATA$X
Y.full <- FULL.DATA$Y
d <- length(beta0)
hist(Y.full)
data <- as.data.frame(cbind(Y.full, X.full))
formula <- Y.full ~ . # Y.full ~ . -1
n.plt <- 200
n.ssp <- 500
FULL <- glm.coef.estimate(X = cbind(1, X.full), Y = Y.full,  family=family)$beta
MSE_FULL <- sum((FULL - beta0)^2)

Weighted.optL.results <-  glm.optimal.subsampling(formula,
                                                  data,
                                                  family = family$family.name,
                                                  n.plt,
                                                  n.ssp,
                                                  estimate.method = 'Weighted',
                                                  sampling.method = 'WithReplacement', #WithReplacement Poisson
                                                  criterion = 'optA')
MSE <- sum((Weighted.optL.results$beta.cmb - beta0)^2)


subsampling.summary(Weighted.optL.results)
##
Weighted.optL.results <-  glm.optimal.subsampling(formula,
                                                  data,
                                                  family = family$family.name,
                                                  n.plt,
                                                  n.ssp,
                                                  estimate.method = 'Uni')


subsampling.summary(Weighted.optL.results)
###############################################################################
formula <- paste(colnames(data)[1], "~",
                 paste(colnames(data)[-1], collapse = "+"), "-1") # use '-1' to avoid adding intercept column again.
weights = 1
design <- survey::svydesign(ids =  ~ 1,
                            weights =  ~ weights,
                            data = data)
results <- survey::svyglm(as.formula(formula),
                          design = design,
                          family = poisson(link = "log"))
results <- getMSLE(cbind(1, X.full), Y = Y.full, start = rep(0, d), weights = weights)
results$coefficients
FULL
