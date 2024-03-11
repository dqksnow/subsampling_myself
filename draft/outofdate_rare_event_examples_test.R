source('R/rare_event_functions.R')
# source('rare_event_functions.R')
source('draft/generate_rare_logit_data.R')
rm(list = setdiff(ls(), lsf.str())) # remove all variables except functions

N <-  5 * 1e5
beta0 <- c(NA, -rep(1, 6))
d <- length(beta0)
case <- 3
set.seed(1)
FULL.DATA <- generate.rare.data(N, case=case, beta0)
X.full <- FULL.DATA$X
Y.full <- FULL.DATA$Y
beta0 <- FULL.DATA$beta0
print(paste('mean: ', mean(Y.full), 'sum: ', sum(Y.full)))
n0 <- 200
nss <- 1000
Betas <- matrix(NA, d, 6)
# Betas[,6] <- RareLogistic(X.full, Y.full, n0, nss, method = 'lcc')$beta.est
# RareLogistic

N1 <- sum(Y.full)
pilot.estimate.results <- pilot.estimate(X.full, Y.full, n0)
beta.pilot <- pilot.estimate.results$beta.pilot
# MN.pilot <- pilot.estimate.results$MN.pilot
pi.P <- pilot.estimate.results$pi.P
pi.Plcc <- pilot.estimate.results$pi.Plcc
pilot.index <- pilot.estimate.results$pilot.indx
idx <- runif(N) <= Y.full + (1 - Y.full) * nss * pi.P  # subsample index
x <- X.full[idx, ]
y <- Y.full[idx]
p <- pmin(nss * pi.P[idx],1)
py <- y + (1 - y) * p

beta.popt <- logistic.coef.estimate(X = x, Y = y, weights = 1 / py)$beta
beta.slik <- logistic.coef.estimate(
    X = x,
    Y = y,
    offset = -log(p),
    start = beta.pilot
  )$beta
# slik.seed <- 1
# fit <- try(
#   beta.slik <- logistic.coef.estimate(X = x, Y = y, offset = -log(p))$beta,
#   silent=TRUE
# )
# while ("try-error" %in% class(fit)) {
#   print(paste("method 'slik', try another seed ", slik.seed))
#   slik.seed <- slik.seed+1
#   set.seed(slik.seed)
#   fit <- try(
#     beta.slik <- logistic.coef.estimate(X = x, Y = y, offset = -log(p))$beta,
#     silent=TRUE
#   )
# }
beta.popt
beta.slik
