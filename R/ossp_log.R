Hest <- function(vec, b, r, N){
  quantile(vec, 1 - r/(b*N))
}

MN <- function(X, pbeta, pinv){
  phi <- pbeta * (1 - pbeta)
  t(X) %*% (X * phi * pinv)
}

Psi <- function(X, y, pbeta, pinv){
  psi <- (y - pbeta)^2
  t(X) %*% (X * psi * pinv^2)
}


pbeta <- function(X, beta){
  1 - 1 / (1 + exp(c(X %*% beta)))
}

ossp_num_optA <- function(X, y, pbeta, MN){
  sqrt((y - pbeta)^2 * colSums(solve(MN, t(X))^2))
}

ossp_num_optL <- function(X, y, pbeta){
  sqrt((y - pbeta)^2 * rowSums(X^2))
}

pbeta_multi <- function(X, beta){
  exp(X %*% beta - matrixStats::rowLogSumExps(cbind(0, X %*% beta)))
}

MN_multi <- function(X, pbeta, k, d, pinv){
  X_kd <- X[, rep(1:d, k)]
  diagp <- pbeta[, rep(1:k, each = d)]
  spinv <- sqrt(pinv)
  diagpX_pd <- (diagp * X_kd) * spinv
  MN_p1 <- t(X_kd * spinv) %*% diagpX_pd

  MN_p2 <- t(diagpX_pd) %*% (diagpX_pd)
  MN <- MN_p1 * (diag(1, k) %x% matrix(1, nrow = d, ncol = d)) - MN_p2
  MN
}


