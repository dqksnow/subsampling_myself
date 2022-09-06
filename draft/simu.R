
library(mvtnorm)
n <- 1e4
beta0  <- c(rep(0.5, 7))
d <- length(beta0)
corr  <- 0.5
sigmax  <- matrix(corr, d-1, d-1) + diag(1-corr, d-1)

set.seed(123)
X <- rmvnorm(n, rep(0, d-1), sigmax)
X <- cbind(1, X)
P <- 1 - 1 / (1 + exp(X %*% beta0))
y <- rbinom(n, 1, P)
table(y)

library("survey")

## test use log(n_0/n_1) or log(n_1/n_0)
## use log(n_1/n_0)
n <- 1e4
beta0  <- c(rep(1, 7))
d <- length(beta0)
corr  <- 0.5
sigmax  <- matrix(corr, d-1, d-1) + diag(1-corr, d-1)

set.seed(123)
X <- rmvnorm(n, rep(0, d-1), sigmax)
X <- cbind(1, X)
P <- 1 - 1 / (1 + exp(X %*% beta0))
y <- rbinom(n, 1, P)
table(y)

pilot_ssp <- proptional_ssp(N, k, y)
beta_pilot <- rep(0, d)
for (i in 1:1000){
  pilot_indx <- swr_indx(N, 500, pilot_ssp)
  pinv_pilot <- 1/pilot_ssp[pilot_indx]
  beta_pilot <- beta_pilot + logistic_coef_estimate(X, y, 1, pilot_indx)
}
beta_pilot <- beta_pilot/1000
log(sum(y)/(N - sum(y)))
beta_full <- logistic_coef_estimate(X, y, 1, 1:N)




## rare logistic

library(mvtnorm)
n <- 1e4
beta0  <- c(rep(0.5, 7))
beta0[1] <- -8
d <- length(beta0)
corr  <- 0.5
sigmax  <- matrix(corr, d-1, d-1) + diag(1-corr, d-1)

set.seed(123)
X <- rmvnorm(n, rep(0, d-1), sigmax)
X <- cbind(1, X)
P <- 1 - 1 / (1 + exp(X %*% beta0))
y <- rbinom(n, 1, P)
table(y)

criteria <- "optA"
r0 <- 200
r <- 500

glm(y[second_indx]~X[second_indx,]-1,
    family = "binomial", offset = -log(r * ossp[second_indx]))


logistic_coef_estimate(X, y, 1/(ossp[second_indx]), second_indx)




