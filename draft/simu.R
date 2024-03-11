
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

beta_full <- logistic_coef_estimate(X, y, 1, 1:length(y))

library("survey")

r0 <- 300
r <- 500
logistic_optimal_subsampling(X, y, r0, r,
                             criteria = "optL",
                             method = "Poisson",
                             unweighted.estimator = F,
                             b = 2)

mse_optA_swr <- var_optA_swr <- 0
for(i in 1:1000){
  set.seed(i)
  beta_optA_swr <- logistic_optimal_subsampling(X, y, r0, r,
                                                criteria = "optA",
                                                method = "SWR",
                                                unweighted.estimator = F,
                                                b = 2)
  mse_optA_swr <- mse_optA_swr + sum((beta_optA_swr$beta - beta_full)^2)
  var_optA_swr <- var_optA_swr + sum(beta_optA_swr$var_beta)
  cat(i)
}
mse_optA_swr # F: 69.89392 # T: 64.57807
var_optA_swr

mse_optL_swr <- 0
for(i in 1:1000){
  set.seed(i)
  beta_optL_swr <- logistic_optimal_subsampling(X, y, r0, r,
                                                criteria = "optL",
                                                method = "SWR",
                                                unweighted.estimator = T,
                                                b = 2)
  mse_optL_swr <- mse_optL_swr + sum((beta_optL_swr$beta - beta_full)^2)
  cat(i)
}
mse_optL_swr # 81.15607 # 72.02926

mse_LCC_swr <- 0
for(i in 1:1000){
  set.seed(i)
  beta_LCC_swr <- logistic_optimal_subsampling(X, y, r0, r,
                                                criteria = "LCC",
                                                method = "SWR",
                                                unweighted.estimator = F,
                                                b = 2)
  mse_LCC_swr <- mse_LCC_swr + sum((beta_LCC_swr - beta_full)^2)
  cat(i)
}
mse_LCC_swr # 79.56938 # 74.99367

mse_optA_poisson <- 0
for(i in 1:1000){
  set.seed(i)
  beta_optA_poi <- logistic_optimal_subsampling(X, y, r0, r,
                                                criteria = "optA",
                                                method = "Poisson",
                                                unweighted.estimator = T,
                                                b = 2)
  mse_optA_poisson <- mse_optA_poisson + sum((beta_optA_poi - beta_full)^2)
  cat(i)
}
mse_optA_poisson # F: 55.57619 # T: 47.57418


mse_optL_poisson <- 0
for(i in 1:1000){
  set.seed(i)
  beta_optL_poi <- logistic_optimal_subsampling(X, y, r0, r,
                                                criteria = "optL",
                                                method = "Poisson",
                                                unweighted.estimator = T,
                                                b = 2)
  mse_optL_poisson <- mse_optL_poisson + sum((beta_optL_poi - beta_full)^2)
  cat(i)
}
mse_optL_poisson # 57.1524 # 49.71734


mse_LCC_poisson <- 0
for(i in 1:1000){
  set.seed(i)
  beta_LCC_poi <- logistic_optimal_subsampling(X, y, r0, r,
                                                criteria = "LCC",
                                                method = "Poisson",
                                                unweighted.estimator = F,
                                                b = 2)
  mse_LCC_poisson <- mse_LCC_poisson + sum((beta_LCC_poi - beta_full)^2)
  cat(i)
}
mse_LCC_poisson # 76.21859 # 81.89863


mse_uniform <- 0
for(i in 1:1000){
  beta_uniform <- logistic_coef_estimate(X, y, 1,
                                         sample(1:n, r0+r, replace = T))
  mse_uniform <- mse_uniform + sum((beta_uniform - beta_full)^2)
  cat(i)
}
mse_uniform # 113.5921

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

beta_full <- logistic_coef_estimate(X, y, 1, 1:length(y))

mse_optA <- rep(0, 1000)
for (i in 1:1000) {
  set.seed(i)
  tryCatch({
    beta_optA <- rareLogistic(X, y, r0, r, criteria = "optA")
    mse_optA[i] <- sum((beta_optA - beta_full)^2)
  },
  warning = function(w) {cat("warning")})
  cat(i)
}
mean(mse_optA)
# 0.07050239

mse_optL <- rep(0, 1000)
for (i in 1:1000) {
  set.seed(i)
  tryCatch({
    beta_optL <- rareLogistic(X, y, r0, r, criteria = "optL")
    mse_optL[i] <- sum((beta_optL - beta_full)^2)
  },
  warning = function(w) {cat("warning")})
  cat(i)
}
mean(mse_optL)
# 1.139091


## softmax regression

d <- 3
k <- 2
N <- 1e4

beta0_true <- rep(0, d)
beta1_true <- rep(1, d)
beta2_true <- rep(2, d)

## mzNorm
library(MASS)
set.seed(1)
mu <- rep(0, d)
diagonal <- as.vector(rep(1,d))
sigma <- matrix(0.5, nrow = d, ncol = d)
diag(sigma) <- diagonal
X <- mvrnorm(N, mu, sigma)

prob <- cbind(exp(X %*% beta0_true),
              exp(X %*% beta1_true),
              exp(X %*% beta2_true))
set.seed(1)
mChoices = t(apply(prob, 1, rmultinom, n = 1, size = 1))
y = apply(mChoices, 1, function(x) which(x==1)-1)
table(y)

## check nnet package multinom function
library(nnet)
fit_basic <- multinom(y ~ X - 1)
summary(fit_basic)
NR_multi2(X, y, 1)



## check X %*% solve(MN) and

a <- rowSums((X %*% solve(MN1))^2)
b <- colSums(solve(MN1, t(X))^2)

library(microbenchmark)
microbenchmark(rowSums((X %*% solve(MN1))^2), colSums(solve(MN1, t(X))^2))




