



logistic_optimal_subsampling <- function(X, y, r0, r,
                                         criteria = c("optA", "optL", "LCC"),
                                         method = c("SWR", "Poisson"),
                                         unweighted.estimator = F,
                                         b){
  criteria <- match.arg(criteria)
  method <- match.arg(method)

  N <- length(y)
  k <- length(unique(y)) - 1
  pilot_ssp <- proptional_ssp(N, k, y)

  if (unweighted.estimator == F){
    if (method == "SWR"){
      pilot_indx <- swr_indx(N, r0, pilot_ssp)
      pinv_pilot <- 1/pilot_ssp[pilot_indx]
      beta_pilot <- logistic_coef_estimate(X, y, pinv_pilot, pilot_indx)
      pbeta_pilot  <- pbeta(X, beta_pilot)
      if (criteria == "optA"){
        MN <- MN(X[pilot_indx, ], pbeta_pilot[pilot_indx], pinv_pilot)
        ossp <- ossp_num_optA(X, y, pbeta_pilot, MN)
      } else if (criteria == "optL"){
        ossp <- ossp_num_optL(X, y, pbeta_pilot)
      } else {
        ossp <- abs(y - pbeta_pilot)
      }
      ossp <- c(ossp / sum(ossp))
      second_indx <- swr_indx(N, r, ossp)
    } else {
      pilot_indx <- poisson_indx(N, r0, pilot_ssp)
      pinv_pilot <- 1/(r * pilot_ssp[pilot_indx])
      beta_pilot <- logistic_coef_estimate(X, y, pinv_pilot, pilot_indx)
      pbeta_pilot  <- pbeta(X, beta_pilot)
      if (criteria == "optA"){
        MN <- MN(X[pilot_indx, ], pbeta_pilot[pilot_indx], pinv_pilot)
        ossp <- ossp_num_optA(X, y, pbeta_pilot, MN)
        H <- Hest(ossp_pilot, b, r, N)
        ossp[ossp > H] <- H
        ossp_pilot <- ossp[pilot_indx]
        NPhi <- sum(ossp_pilot * pinv_pilot)
        ossp <- ossp/NPhi
      } else if (criteria == "optL"){
        ossp <- ossp_num_optL(X, y, pbeta_pilot)
        H <- Hest(ossp_pilot, b, r, N)
        ossp[ossp > H] <- H
        ossp_pilot <- ossp[pilot_indx]
        NPhi <- sum(ossp_pilot * pinv_pilot)
        ossp <- ossp/NPhi
      } else {
        ossp <- abs(y - pbeta_pilot)
        ossp <- c(ossp / sum(ossp))
      }
      second_indx <- poisson_indx(N, r, ossp)
    }
    pinv_cmb <- 1/c(pilot_ssp[pilot_indx], ossp[second_indx])
    beta_cmb <- logistic_coef_estimate(X, y, pinv_cmb,
                                       c(pilot_indx, second_indx))
  } else {
    if (method == "SWR"){
      pilot_indx <- swr_indx(N, r0, pilot_ssp)
      pinv_pilot <- 1/pilot_ssp[pilot_indx]
      beta_pilot <- logistic_coef_estimate(X, y, 1, pilot_indx) +
        log(sum(y)/(N - sum(y)))
      pbeta_pilot  <- pbeta(X, beta_pilot)
      if (criteria == "optA"){
        MN <- MN(X[pilot_indx, ], pbeta_pilot[pilot_indx], pinv_pilot)
        ossp <- ossp_num_optA(X, y, pbeta_pilot, MN)
      } else if (criteria == "optL"){
        ossp <- ossp_num_optL(X, y, pbeta_pilot)
      } else {
        ossp <- abs(y - pbeta_pilot)
      }
      ossp <- c(ossp / sum(ossp))
      MN_pilot <- MN(X[pilot_indx, ], pbeta_pilot[pilot_indx], 1)
      second_indx <- swr_indx(N, r, ossp)
    } else {
      pilot_indx <- poisson_indx(N, r0, pilot_ssp)
      pinv_pilot <- 1/(r * pilot_ssp[pilot_indx])
      beta_pilot <- logistic_coef_estimate(X, y, 1, pilot_indx) +
        log(sum(y)/(N - sum(y)))
      pbeta_pilot  <- pbeta(X, beta_pilot)
      if (criteria == "optA"){
        MN <- MN(X[pilot_indx, ], pbeta_pilot[pilot_indx], pinv_pilot)
        ossp <- ossp_num_optA(X, y, pbeta_pilot, MN)
        H <- Hest(ossp_pilot, b, r, N)
        ossp[ossp > H] <- H
        ossp_pilot <- ossp[pilot_indx]
        NPhi <- sum(ossp_pilot * pinv_pilot)
        ossp <- ossp/NPhi
      } else if (criteria == "optL"){
        ossp <- ossp_num_optL(X, y, pbeta_pilot)
        H <- Hest(ossp_pilot, b, r, N)
        ossp[ossp > H] <- H
        ossp_pilot <- ossp[pilot_indx]
        NPhi <- sum(ossp_pilot * pinv_pilot)
        ossp <- ossp/NPhi
      } else {
        ossp <- abs(y - pbeta_pilot)
        ossp <- c(ossp / sum(ossp))
      }
      second_indx <- poisson_indx(N, r, ossp)
    }
    beta_s2 <- logistic_coef_estimate(X, y, 1, second_indx) + beta_pilot
    pbeta_s2  <- pbeta(X, beta_s2)
    MN_s2 <- MN(X[second_indx, ], pbeta_s2[second_indx], 1)
    beta_cmb <- c(solve(MN_pilot + MN_s2,
                        MN_pilot %*% beta_pilot + MN_s2 %*% pbeta_s2))
  }
  beta_cmb
}






rareLogistic <- function(X, y){

}







