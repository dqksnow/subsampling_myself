
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



NBIBI <- function(X, k, d, r0, pbeta, pinv){
  NBIBI <- matrix(0, nrow = k*d, ncol = k*d)
  p0 <- 1 - rowSums(pbeta)
  pi_num <- rep(0, N)
  for(obs in 1 : r0){
    phi <- diag(pbeta[obs,]) - pbeta[obs,] %*% t(pbeta[obs,])
    fdP <- rbind(-p0[obs] * pbeta[obs,], phi) %x% t(X[obs, ])
    NBIBI <- NBIBI + t(fdP) %*% fdP * pinv[obs]
  }
  NBIBI
}



# criteria = c("optA", "optL", "mMSPE", "custom")
# method = c("SWR", "Poisson")
# constraint = c("baseline", "summation")

algorithm2_multi <- function(X, y, r0, r,
                             criteria = "optL",
                             method = "Poisson",
                             constraint = "baseline",
                             MN_custom,
                             b = 2){

  k <- length(unique(y))-1
  N <- nrow(X)
  d <- ncol(X)
  I <- matrix(0, nrow = N, ncol = k)
  I[cbind(seq_along(y), y)] <- 1
  pilot_ssp <- proptional_ssp(N, k, y)

  if (method == "SWR"){
    pilot_indx <- swr_indx(N, r0, pilot_ssp)
  } else if (method == "Poisson"){
    pilot_indx <- poisson_indx(N, r0, pilot_ssp)
  }
  pinv_pilot <- 1/pilot_ssp[pilot_indx]
  invisible(capture.output(
    beta_pilot <- softmax_coef_estimate(X, y, pinv_pilot, pilot_indx)))
  pbeta_pilot <- pbeta_multi(X, beta_pilot)
  G <- rbind(rep(-1/(k+1), k), diag(k) - 1/(k+1)) %x% diag(d)
  MN1 <- MN_multi(X[pilot_indx,], pbeta_pilot[pilot_indx,],
                  k, d, pinv = pinv_pilot)
  psi1 <- psi_multi(X[pilot_indx,], pilot_ssp[pilot_indx],
                    pbeta_pilot[pilot_indx,], r0, k, d, I[pilot_indx,])

  if (criteria == "optA"){
    sixi <- (I-pbeta_pilot)[,rep(seq(1:k), each = ncol(X))] *
      X[,rep(seq(ncol(X)), k)]

    if (constraint == "baseline"){
      pi_num <- sqrt(colSums((solve(MN1, t(sixi))^2)))
    } else if (constraint == "summation"){
      pi_num <- sqrt(colSums(((G %*% solve(MN1, t(sixi)))^2)))
    }
  } else if (criteria == "optL") {
    if (constraint == "baseline"){
      pi_num <- sqrt(rowSums((I-pbeta_pilot)^2)* rowSums(X^2))
    } else if (constraint == "summation"){
      pi_num <- sqrt((rowSums((I-pbeta_pilot)^2)+
                            (rowSums(pbeta_pilot)-rowSums(I))^2) *
                           rowSums(X^2))
    }
  } else if (criteria == "mMSPE") {
    NBIBI_pilot <- NBIBI(X[pilot_indx,], k, d, r0,
                         pbeta_pilot[pilot_indx,], pinv_pilot)
    sixi <- (I-pbeta_pilot)[,rep(seq(1:k), each = ncol(X))] *
      X[,rep(seq(ncol(X)), k)]
    Msixi <- solve(MN, t(sixi))
    pi_num <- sqrt(colSums((NBIBI_pilot %*% Msixi) * Msixi))
  } else if (criteria == "custom") {
    sixi <- (I-pbeta_pilot)[,rep(seq(1:k), each = ncol(X))] *
      X[,rep(seq(ncol(X)), k)]
    pi_num <- sqrt(colSums((solve(MN_custom, t(sixi))^2)))
  }

  if (method == "SWR"){
    ossp <- pi_num/sum(pi_num)
    second_indx <- swr_indx(N, r, ossp)
  } else if (method == "Poisson"){
    ti1 <- pi_num[pilot_indx]
    H <- quantile(ti1, 1-r/(b*N))
    ti1[ti1 > H] <- H
    ti.num <- (length(pilot_indx)/(length(pilot_indx) - d*k))*
      sum(ti1 * pinv_pilot/length(pilot_indx))
    pi_num[pi_num > H] <- H
    ossp <- pi_num/ti.num
    second_indx <- poisson_indx(N, r, ossp)
  }
  pinv_s2 <- 1/ossp[second_indx]
  invisible(capture.output(
    beta_s2 <- softmax_coef_estimate(X, y, pinv_s2, second_indx)))
  pbeta_s2  <- pbeta(X[second_indx,], beta_s2)

  MN2 <- MN_multi(X[second_indx,], pbeta_pilot[second_indx,],
                  k, d, pinv = pinv_s2)
  psi2 <- psi_multi(X[second_indx,], pilot_ssp[second_indx],
                    pbeta_pilot[second_indx,], r, k, d, I[second_indx,])

  beta_cmb <- solve(MN1 + MN2, (MN1 %*% c(beta_pilot) + MN2 %*% c(beta_s2)))
  beta_cmb <- matrix(beta_cmb, ncol = k)
  var <- (solve(MN1 + MN2, psi1 + psi2) %*% solve(MN1 + MN2))
  trvar <- diag(var)
  return(list(beta = beta_cmb,
              var_beta = trvar,
              ossp = ossp,
              pilot_indx = pilot_indx,
              second_indx = second_indx))
}




