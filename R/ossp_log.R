Hest <- function(vec, b, r, N){
  quantile(vec, 1 - r/(b*N))
}

MN <- function(X, pbeta, pinv){
  phi <- pbeta * (1 - pbeta)
  t(X) %*% (X * phi * pinv)
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


# logistic_swr_ossp <- function(X, y, pilot_indx,
#                               pinv_pilot, beta_pilot,
#                               criteria = c("optA", "optL")){
#
#   pbeta  <- pbeta(X, beta_pilot)
#   if (criteria == "optA"){
#     MN <- MN(X[pilot_indx, ], pbeta[pilot_indx], pinv_pilot)
#     ossp <- ossp_num_optA(X, y, pbeta, MN)
#   } else if (criteria == "optL"){
#     ossp <- ossp_num_optL(X, y, pbeta)
#   } else {
#     ossp <- abs(y - pbeta)
#   }
#   ossp <- c(ossp / sum(ossp))
#   ossp
# }
#
#
#
# logistic_poisson_ossp <- function(X, y, pilot_indx,
#                                   pinv_pilot, beta_pilot,
#                                   criteria = c("optA", "optL"), b, r, N){
#   pbeta  <- pbeta(X, beta_pilot)
#   if (criteria == "optA"){
#     MN <- MN(X[pilot_indx, ], pbeta[pilot_indx], pinv_pilot)
#     ossp <- ossp_num_optA(X, y, pbeta, MN)
#     H <- Hest(ossp_pilot, b, r, N)
#     ossp[ossp > H] <- H
#     ossp_pilot <- ossp[pilot_indx]
#     NPhi <- sum(ossp_pilot * pinv_pilot)
#     ossp <- ossp/NPhi
#   } else if (criteria == "optL"){
#     ossp <- ossp_num_optL(X, y, pbeta)
#     H <- Hest(ossp_pilot, b, r, N)
#     ossp[ossp > H] <- H
#     ossp_pilot <- ossp[pilot_indx]
#     NPhi <- sum(ossp_pilot * pinv_pilot)
#     ossp <- ossp/NPhi
#   } else {
#     ossp <- abs(y - pbeta)
#     ossp <- c(ossp / sum(ossp))
#   }
#
# }
#
#

