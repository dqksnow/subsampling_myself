source('R/family_expand.R')
source('R/glm_subsampling_functions.R')
source('draft/generate_poisson_data.R')
rm(list = setdiff(ls(), lsf.str()))
N <- 1e5
p <- 6
beta0 <- rep(0.5, p + 1)
shape <- 2
seed <- 1
FULL.DATA <- generate.gamma.data(N = N, case = 1, beta0 = beta0, shape = shape, seed = seed)
X.full <- FULL.DATA$X
Y.full <- FULL.DATA$Y
data <- as.data.frame(cbind(Y.full, X.full))
formula <- Y.full ~ .
start <- rep(0.1, p+1)
n.plt <- 200
n.ssp <- 500
# m3 <- glm(formula = formula,
#           data = data,
#           start = start,
#           family = Gamma(link = "inverse"))
# summary(m3)

##########################################
family <- Gamma(link = "inverse")

formula <- as.formula(paste(colnames(data)[1], "~",
                            paste(colnames(data)[-1], collapse = "+"))) # with intercept
design <- survey::svydesign(ids = ~ 1,
                            weights =  ~ 1,
                            data = data)
results <- survey::svyglm(formula = formula,
                          design = design,
                          family = family)
summary(results)
##########################################
family = gamma.expand()
index.plt <- random.index(N, n.plt + n.ssp)
x.plt <- X.full[index.plt,]
x.plt <- cbind(1, x.plt)
y.plt <- Y.full[index.plt]
data <- as.data.frame(cbind(y.plt, x.plt))
# unweighted log likelihood estimation
beta.plt <- glm.coef.estimate(X = x.plt, Y = y.plt, family = family)$beta
beta.plt
######################################################
family <- gamma.expand()
formula <- Y.full ~ .
data <- as.data.frame(cbind(Y.full, X.full))
Wet_WR_optA <-  glm.optimal.subsampling(formula, data, n.plt, n.ssp,
                                        family = family$family.name,
                                         estimate.method = 'Uni', # LogOddsCorrection
                                        alpha = 0,
                                        b=2)
Wet_WR_optA$beta
######################################################
family <- gamma.expand()
LOC_Poi_optA <-  glm.optimal.subsampling(formula, data, n.plt, n.ssp,
                                        family = family$family.name,
                                        criterion = 'OptL',
                                        sampling.method = 'Poisson',
                                        estimate.method = 'Weighted',
                                        alpha = 0,
                                        b=2)
LOC_Poi_optA$beta.cmb
######################################################
