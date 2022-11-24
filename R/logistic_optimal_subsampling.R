

# criteria = c("optA", "optL", "LCC", "custom")
# method = c("SWR", "Poisson")

logistic_optimal_subsampling <- function(X, y, r0, r,
                                         criteria = "optL",
                                         method = "Poisson",
                                         MN_custom,
                                         unweighted.estimator = F,
                                         b = 2){
  # criteria <- match.arg(criteria)
  # method <- match.arg(method)

  N <- length(y)
  k <- length(unique(y)) - 1
  pilot_indx <- plt_indx(N, y, r0)
  pinv_pilot <- 1/c(rep(1/(2 * (N - sum(y))), r0/2),
                    rep(1/(2 * sum(y)), r0/2))

  if (unweighted.estimator == F){
    beta_pilot <- logistic_coef_estimate(X, y, pinv_pilot, pilot_indx)
    pbeta_pilot  <- pbeta(X, beta_pilot)
    MN1 <- MN(X[pilot_indx, ], pbeta_pilot[pilot_indx], pinv_pilot)
    Psi1 <- Psi(X[pilot_indx, ], y[pilot_indx],
                pbeta_pilot[pilot_indx], pinv_pilot)
    if (criteria == "optA"){
      ossp <- ossp_num_optA(X, y, pbeta_pilot, MN1)
    } else if (criteria == "optL"){
      ossp <- ossp_num_optL(X, y, pbeta_pilot)
    } else if (criteria == "LCC"){
      ossp <- abs(y - pbeta_pilot)
    } else if (criteria == "custom"){
      ossp <- ossp_num_optA(X, y, pbeta_pilot, MN_custom)
    }
    if (method == "SWR"){
      ossp <- c(ossp / sum(ossp))
      second_indx <- swr_indx(N, r, ossp)
    } else if (method == "Poisson"){
      H <- Hest(ossp, b, r, N)
      ossp[ossp > H] <- H
      ossp_pilot <- ossp[pilot_indx]
      NPhi <- sum(ossp_pilot * pinv_pilot/r0)
      ossp <- ossp/NPhi
      second_indx <- poisson_indx(N, r, ossp)
    }
    pinv_s2 <- 1/(ossp[second_indx])
    beta_s2 <- logistic_coef_estimate(X, y, pinv_s2, second_indx)
    pbeta_s2  <- pbeta(X[second_indx,], beta_s2)
    MN_s2 <- MN(X[second_indx, ], pbeta_s2, pinv_s2)
    Psi2 <- Psi(X[second_indx, ], y[second_indx], pbeta_s2, pinv_s2)
  } else if (unweighted.estimator == T){
    beta_pilot_llk <- beta_pilot <-
      logistic_coef_estimate(X, y, 1, pilot_indx)
    beta_pilot[1] <- beta_pilot[1] + log(sum(y)/(N - sum(y)))
    pbeta_pilot  <- pbeta(X, beta_pilot)

    pbeta_pilot_llk <- pbeta(X[pilot_indx, ], beta_pilot_llk)
    MN1 <- MN(X[pilot_indx, ], pbeta_pilot_llk, 1)
    Psi1 <- Psi(X[pilot_indx, ], y[pilot_indx], pbeta_pilot_llk, 1)

    if (criteria == "optA"){
      MN_pilot <- MN(X[pilot_indx, ], pbeta_pilot[pilot_indx], pinv_pilot)
      ossp <- ossp_num_optA(X, y, pbeta_pilot, MN_pilot)
    } else if (criteria == "optL"){
      ossp <- ossp_num_optL(X, y, pbeta_pilot)
    } else if (criteria == "LCC"){
      ossp <- abs(y - pbeta_pilot)
    } else if (criteria == "custom"){
      ossp <- ossp_num_optA(X, y, pbeta_pilot, MN_custom)
    }
    if (method == "SWR"){
      ossp <- c(ossp / sum(ossp))
      second_indx <- swr_indx(N, r, ossp)
    } else if (method == "Poisson"){
      H <- Hest(ossp, b, r, N)
      ossp[ossp > H] <- H
      ossp_pilot <- ossp[pilot_indx]
      NPhi <- sum(ossp_pilot * pinv_pilot/r0)
      ossp <- ossp/NPhi
      second_indx <- poisson_indx(N, r, ossp)
    }
    beta_s2_llk <- logistic_coef_estimate(X, y, 1, second_indx)
    beta_s2 <- beta_s2_llk + beta_pilot
    pbeta_s2_llk  <- pbeta(X[second_indx,], beta_s2_llk)
    MN_s2 <- MN(X[second_indx, ], pbeta_s2_llk, 1)
    Psi2 <- Psi(X[second_indx, ], y[second_indx], pbeta_s2_llk, 1)
  }
  MNsolve <- solve(MN1 + MN_s2)
  beta_cmb <- c(MNsolve %*% (MN1 %*% beta_pilot + MN_s2 %*% beta_s2))
  var_beta <- MNsolve %*% (Psi1+Psi2) %*% MNsolve
  return(list(beta = beta_cmb,
              var_beta = var_beta,
              ossp = ossp,
              pilot_indx = pilot_indx,
              second_indx = second_indx))
}






rareLogistic <- function(X, y, r0, r,
                         criteria = c("optA", "optL")){
  criteria <- match.arg(criteria)

  N <- length(y)
  k <- length(unique(y)) - 1
  pilot_ssp <- proptional_ssp(N, k, y)

  pilot_indx <- poisson_indx(N, r0, pilot_ssp)
  pinv_pilot <- 1/(length(pilot_indx) * pilot_ssp[pilot_indx])
  theta_pilot <- logistic_coef_estimate(X, y, 1, pilot_indx)
  theta_pilot[1] <- theta_pilot[1] + log(sum(y)/(N - sum(y)))

  pbeta_pilot <- exp(c(X[,-1] %*% theta_pilot[-1]))
  ptheta_pilot <- 1 - 1 / (1 + exp(theta_pilot[1]) * pbeta_pilot)

  if(criteria == "optA"){
    # Mf <- t(X) %*% (X * pbeta_pilot)
    Mf_tilde <- t(X[pilot_indx,]) %*%
      (X[pilot_indx,] * pbeta_pilot[pilot_indx] * pinv_pilot)
    # Mf_tilde <- MN(X[pilot_indx,], ptheta_pilot[pilot_indx], pinv_pilot)
    num <- ptheta_pilot * sqrt(rowSums((X %*% solve(Mf_tilde))^2))
    omega <- sum(ptheta_pilot[pilot_indx] *
                   sqrt(rowSums((X[pilot_indx,] %*% solve(Mf_tilde))^2)) *
                   pinv_pilot)
    ossp <- num/omega
  } else {
    num <- ptheta_pilot * sqrt(rowSums(X^2))
    omega <- sum(num[pilot_indx] * pinv_pilot)
    ossp <- num/omega
  }
  ossp[y==1] <- 1/r
  ossp[r * ossp > 1] <- 1/r
  second_indx <- poisson_indx(N, r, ossp)

  glm(y[second_indx]~X[second_indx,]-1,
      family = "binomial",
      offset = -log(r * ossp[second_indx]))$coefficients
}





