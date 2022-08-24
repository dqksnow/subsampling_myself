


logistic_optimal_subsampling <- function(X, y, r0, r,
                                         criteria = c("optA", "optL"),
                                         method = c("SWR", "Poisson"),
                                         b){
  criteria <- match.arg(criteria)
  method <- match.arg(criteria)

  N <- length(y)
  k <- length(unique(y)) - 1
  pilot_ssp <- proptional_ssp(N, k, y)

  if(method == "SWR"){
    pilot_indx <- swr_indx(N, r0, pilot_ssp)
    ossp <- logistic_swr_ossp(X, y,
                              pilot_ssp, pilot_indx,
                              criteria = criteria)
    second_indx <- swr_indx(N, r, ossp)
  } else {
    pilot_indx <- poisson_indx(N, r0, pilot_ssp)
    ossp <- logistic_swr_ossp(X, y,
                              pilot_ssp, pilot_indx,
                              criteria = criteria, b, r, N)
    second_indx <- poisson_indx(N, r, ossp)
  }

  pinv_cmb <- 1/c(pilot_ssp[pilot_indx], ossp[second_indx])
  beta_cmb <- logistic_coef_estimate(X, y, pinv_cmb, c(pilot_indx, second_indx))
  return(beta_cmb)
}



