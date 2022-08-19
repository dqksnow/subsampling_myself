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
  criteria = match.arg(criteria)
  pinv_pilot <- 1/pilot_ssp[pilot_indx]
  beta_pilot <- logistic_coef_estimate(X_pilot, y_pilot, pinv_pilot, pilot_indx)
  pbeta  <- 1 - 1 / (1 + exp(c(X %*% beta_pilot)))

  if (criteria == "optA"){
    pbeta_pilot <- pbeta[pilot_indx]
    phi_pilot <- pbeta_pilot * (1 - pbeta_pilot)
    MN_solve_pilot <- solve(t(X_pilot) %*% (X_pilot * phi_pilot * pinv_pilot))
    ossp <- sqrt((y - pbeta)^2 * rowSums((X %*% MN_solve_pilot)^2))
  } else {
    ossp <- sqrt((y - pbeta)^2 * rowSums(X^2))
  }
  ossp <- c(ossp / sum(PI.opt))
  ossp
}








