
## functions for softmax regression

MN_multi <- function(X, p, k, d, pi){
  X_kd <- X[, rep(1:d, k)]
  diagp <- p[, rep(1:k, each = d)]
  spi <- sqrt(pi)
  diagpX_pd <- (diagp * X_kd)/spi
  MN_p1 <- t(X_kd/spi) %*% diagpX_pd

  MN_p2 <- t(diagpX_pd) %*% (diagpX_pd)
  MN <- MN_p1 * (diag(1, k) %x% matrix(1, nrow = d, ncol = d)) - MN_p2
  MN
}



algorithm2_multi <- function(X, y, r0, r,
                             criteria = c("optA", "optL"),
                             method = c("SWR", "Poisson"),
                             constraint = c("baseline", "summation", "mMSPE"),
                             b = 2){

  k <- length(unique(y))-1
  N <- nrow(X)
  d <- ncol(X)
  I <- matrix(0, nrow = N, ncol = k)
  I[cbind(seq_along(y), y)] <- 1
  pilot_ssp <- proptional_ssp(N, k, y)

  if (method == "SWR"){
    pilot_indx <- swr_indx(N, r0, pilot_ssp)
    pinv_pilot <- 1/pilot_ssp[pilot_indx]
    invisible(capture.output(
      beta_pilot <- softmax_coef_estimate(X, y, pinv_pilot, pilot_indx)))
    pbeta_pilot <- pbeta_multi(X, beta_pilot)
    if (criteria == "optA"){
      MN <- MN_multi(X[pilot_indx,], pbeta_pilot[pilot_indx,],
                     k, d, pinv = pinv_pilot)
      sixi <- (I-pbeta_pilot)[,rep(seq(1:k), each = ncol(X))] *
        X[,rep(seq(ncol(X)), k)]
      pi_mMSE_num <- sqrt(colSums((solve(MN, t(sixi))^2)))
      ossp <- pi_mMSE_num/sum(pi_mMSE_num)
    } else {
      pi_mVc_num <- sqrt(rowSums((I-pbeta_pilot)^2)* rowSums(X^2))
      ossp <- pi_mVc_num/sum(pi_mVc_num)
    }
    second_indx <- swr_indx(N, r, ossp)
  } else {
    pilot_indx <- poisson_indx(N, r0, pilot_ssp)
    pinv_pilot <- 1/pilot_ssp[pilot_indx]
    invisible(capture.output(
      beta_pilot <- softmax_coef_estimate(X, y, pinv_pilot, pilot_indx)))
    pbeta_pilot <- pbeta_multi(X, beta_pilot)
    r0s <- length(pilot_indx)

    if (criteria == "optA"){
      MN <- MN_multi(X[pilot_indx,], pbeta_pilot[pilot_indx,],
                     k, d, pinv = pinv_pilot)
      sixi <- (I-pbeta_pilot)[,rep(seq(1:k), each = ncol(X))] *
        X[,rep(seq(ncol(X)), k)]
      pi_num <- sqrt(colSums((solve(MN, t(sixi))^2)))
      ti1 <- pi_num[pilot_indx]
    } else {
      pi_num <- sqrt(rowSums((I-pbeta_pilot)^2)* rowSums(X^2))
      ti1 <- pi_num[pilot_indx]
    }
    H <- quantile(ti1, 1-r/(b*N))
    ti1[ti1 > H] <- H
    ti.num <- (r0s/(r0s - d*k))*sum(ti1 * pinv_pilot/r0s)

    pi_num[pi_num > H] <- H
    ossp <- pi_num/ti.num
    second_indx <- poisson_indx(N, r, ossp)
  }
  pinv_cmb <- 1/c(pilot_ssp[pilot_indx], ossp[second_indx])
  invisible(capture.output(
    beta_cmb <- softmax_coef_estimate(X, y, pinv_cmb,
                                      c(pilot_indx, second_indx))))
  beta_cmb
}




