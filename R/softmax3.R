proptional_ssp <- function(N, K, y){
  prob <- 1/((K+1) * as.data.frame(table(y))$Freq)
  return(prob[y+1])
}
swr_indx <- function(N, r, pi){
  sample(1:N, r, T, pi)
}
###############################################################################
random.index <- function (N, n, p=NULL) {
  ifelse(is.null(p),
         index <- sample(N, n, replace = TRUE),
         index <- sample(N, n, replace = TRUE, prob = p))
  return(as.vector(index))
}
softmax_coef_estimate <- function(X, y, weights){
  fit <- nnet::multinom(y ~ X - 1, weights = weights)
  t(coef(fit))
}
###############################################################################
pbeta_multi <- function(X, beta){
  exp(X %*% beta - matrixStats::rowLogSumExps(cbind(0, X %*% beta)))
}
###########################
MN_multi <- function(X, p, K, d, pi){
  X_Kd <- X[, rep(1:d, K)]
  diagp <- p[, rep(1:K, each = d)]
  spi <- sqrt(pi)
  diagpX_pd <- (diagp * X_Kd)/spi
  MN_p1 <- t(X_Kd/spi) %*% diagpX_pd
  MN_p2 <- t(diagpX_pd) %*% (diagpX_pd)
  MN <- MN_p1 * (diag(1, K) %x% matrix(1, nrow = d, ncol = d)) - MN_p2
  return(MN)
}
###############################################################################
psi_multi <- function(X, pi, pbeta, r, k, d, I) {
  num <- min(k,2)
  comb <- combn(1:k, num)
  pi_s2_s <- pi/r
  r_pi <- 1 - pi

  psii <- matrix(NA, nrow = d^2, ncol = k^2)
  psii[, (comb[1,]-1)*k+comb[num,]] <-
    psii[,(comb[num,]-1)*k+comb[1,]] <-
    t(X) %*% ((as.matrix((I-pbeta)[,comb[1,]] *
                           (I-pbeta)[,comb[num,]])[,rep(1:dim(comb)[2], each=d)]*
                 (r_pi/(pi_s2_s^2)))*X[,rep(1:d, dim(comb)[2])])
  psii[, seq(1,(k^2), by = (k+1))] <-
    t(X) %*% ((((I-pbeta)^2)[,rep(1:k, each=d)]*
                 (r_pi/(pi_s2_s^2)))*X[,rep(1:d, k)])
  psii <- t(matrix(psii, nrow = d))
  psi <- matrix(psii, nrow = k*d)
  psi <- t(psi[, c(rep(seq(1,k*d,k), k) + rep(0:(k-1), each = d))])
  psi
}
###############################################################################
rm(list = setdiff(ls(), lsf.str()))
d <- 2 # dim of covariates
K <- 4 # K + 1 classes
N <- 1e4
beta0 <- matrix(-1, d, K) # (K+1)
library(MASS)
set.seed(1)
mu <- rep(0, d)
sigma <- matrix(0.5, nrow = d, ncol = d)
diag(sigma) <- rep(1, d)
X <- mvrnorm(N, mu, sigma)
prob.K <- exp(X %*% beta0)
prob.K <- prob.K / (1 + rowSums(prob.K))
prob.Kp1 <- cbind(1 - rowSums(prob.K), prob.K)
Y <- apply(prob.Kp1, 1, function(row) sample(0:K, size = 1, prob = row))
table(Y)
Y.matrix <- matrix(0, nrow = N, ncol = K)
Y.matrix[cbind(seq_along(Y), Y)] <- 1
#
n.plt <- 200
n.ssp <- 500
index.plt <- random.index(N, n.plt)
p.plt <- rep(1 / N, n.plt) # pilot sampling probability
x.plt <- X[index.plt,]
y.plt <- Y[index.plt]
Y.matrix.plt <- Y.matrix[index.plt, ]
beta.plt <- softmax_coef_estimate(x.plt, y.plt, weights = 1 / p.plt) # weigths=1
P.step1 <- pbeta_multi(X, beta.plt) # N*K, i-th row = Prob(Y_i=1) for i= 1 to K
P.plt <- P.step1[index.plt, ]
MN1 <- MN_multi(x.plt, P.plt, K, d, pi = p.plt)
psi1 <- psi_multi(x.plt, pi = p.plt,
                  P.plt, n.plt, K, d, Y.matrix.plt)
sixi <- (Y.matrix - P.step1)[, rep(seq(1:K), each = d)] *
  X[, rep(seq(d), K)]
pi_num <- sqrt(colSums((solve(MN1, t(sixi))^2)))
ossp <- pi_num/sum(pi_num)
index.ssp <- random.index(N, n.ssp, ossp)
p.ssp <- ossp[index.ssp]
x.ssp <- X[index.ssp,]
y.ssp <- Y[index.ssp]
beta.ssp <- softmax_coef_estimate(x.ssp, y.ssp, weights = 1 / p.ssp)
P.step2  <- pbeta_multi(X, beta.ssp)
P.ssp <- P.step2[index.ssp, ]
MN2 <- MN_multi(x.ssp, P.ssp, K, d, pi = p.ssp)
beta_cmb <- solve(MN1 + MN2, (MN1 %*% c(beta.plt) + MN2 %*% c(beta.ssp)))
beta_cmb <- matrix(beta_cmb, ncol = K)
