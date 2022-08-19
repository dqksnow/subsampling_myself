


logistic_optimal_subsampling <- function(X, y, r0, r, 
                                         criteria = c("optA", "optL")){
  criteria = match.arg(criteria)
  
  N <- length(y)
  k <- length(unique(y)) - 1
  
  pilot_ssp <- proptional_ssp(N, k, y)
  pilot_indx <- swr_indx(N, r0, pilot_ssp)
  ossp <- logistic_swr_ossp(X, y, pilot_ssp, pilot_indx, criteria = criteria)
  
  second_indx <- swr_indx(N, r, ossp)
  pinv_cmb <- 1/c(pilot_ssp[pilot_indx], ossp[second_indx])
  beta_cmb <- logistic_coef_estimate(X_cmb, y_cmb, pinv_cmb, cmb_indx)
  return(beta_cmb)
}



