#' Compute optimal subsampling probabilities for
#' logistic regression
#'
#' @param X A N \times d matrix
#' @param criteria Either "optA" or "optL"
#' @param y Response vector
#' @param pilot_ssp Subsampling probabilities obtaining pilot sample
#' @param pilot_indx Data index of the pilot sample
#'
#' @return Optimal subsampling probabilities under A- and L-
#' optimality criterion
#' @export
#'
#' @examples
logistic_swr_ossp <- function(X, y, pilot_ssp, pilot_indx,
                          criteria = c("optA", "optL")){

  pinv_pilot <- 1/pilot_ssp[pilot_indx]
  beta_pilot <- logistic_coef_estimate(X, y, pinv_pilot, pilot_indx)
  pbeta  <- 1 - 1 / (1 + exp(c(X %*% beta_pilot)))

  if (criteria == "optA"){
    MN <- MN(X, pilot_indx, pbeta, pinv_pilot)
    ossp <- ossp_num_optA(X, y, pbeta, MN)
  } else {
    ossp <- ossp_num_optL(X, y, pbeta)
  }
  ossp <- c(ossp / sum(PI.opt))
  ossp
}


Hest <- function(vec, b, r, N){
  quantile(vec, 1 - r/(b*N))
}

MN <- function(X, indx, pbeta, pinv){
  pbeta_smp <- pbeta[indx]
  phi_smp <- pbeta_smp * (1 - pbeta_smp)
  X_smp <- X[indx,]
  t(X_smp) %*% (X_smp * phi_smp * pinv)
}

pbeta <- function(X, beta){
  1 - 1 / (1 + exp(c(X %*% beta)))
}

ossp_num_optA <- function(X, y, pbeta, MN){
  sqrt((y - pbeta)^2 * rowSums((X %*% solve(MN))^2))
}

ossp_num_optL <- function(X, y, pbeta){
  sqrt((y - pbeta)^2 * rowSums(X^2))
}

logistic_poisson_ossp <- function(X, y, pilot_ssp, pilot_indx,
                                  criteria = c("optA", "optL"), b, r, N){
  pinv_pilot <- 1/(r * pilot_ssp[pilot_indx])
  beta_pilot <- logistic_coef_estimate(X_pilot, y_pilot, pinv_pilot, pilot_indx)
  pbeta  <- pbeta(X, beta_pilot)

  if (criteria == "optA"){
    MN <- MN(X, pilot_indx, pbeta, pinv_pilot)
    ossp <- ossp_num_optA(X, y, pbeta, MN)
  } else {
    ossp <- sqrt((y - pbeta)^2 * rowSums(X^2))
    ossp <- ossp_num_optL(X, y, pbeta)
  }
  H <- Hest(ossp_pilot, b, r, N)
  ossp[ossp > H] <- H
  ossp_pilot <- ossp[pilot_indx]
  NPhi <- sum(ossp_pilot * pinv_pilot)
  ossp/NPhi
}



