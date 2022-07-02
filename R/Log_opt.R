
## optimal subsampling for logistic regression

## Questions?
## use first stage sample to estimate M_N or use the full data to compute M_N?

Log_optA <- function(X, y, beta.prop, idx.prop, PI.prop){
  
  x.prop <- X[idx.prop,]
  y.prop <- y[idx.prop]
  pinv.prop <- 1/PI.prop[idx.prop]
  
  P.prop  <- 1 - 1 / (1 + exp(c(X %*% beta.prop)))
  p.prop <- P.prop[idx.prop]
  phi.prop <- p.prop * (1 - p.prop)
  
  W.prop <- solve(t(x.prop) %*% (x.prop * phi.prop * pinv.prop))
  PI.opt <- sqrt((y - P.prop)^2 * rowSums((X %*% W.prop)^2))
  PI.opt <- c(PI.opt / sum(PI.opt))
  PI.opt
}



Log_optL <- function(X, y, beta.prop){
  
  P.prop  <- 1 - 1 / (1 + exp(c(X %*% beta.prop)))
  PI.opt <- sqrt((y - P.prop)^2 * rowSums(X^2))
  PI.opt <- c(PI.opt / sum(PI.opt))
  PI.opt
}


Log_Poi_optA <- function(X, y, beta.prop, idx.prop, PI.prop){
  
  x.prop <- X[idx.prop,]
  y.prop <- y[idx.prop]
  pinv.prop <- 1/pmin(r0 * PI.prop[idx.prop], 1)
  
  P.prop  <- 1 - 1 / (1 + exp(c(X %*% beta.prop)))
  p.prop <- P.prop[idx.prop]
  phi.prop <- p.prop * (1 - p.prop)
  
  W.prop <- solve(t(x.prop) %*% (x.prop * phi.prop * pinv.prop))
  PI.opt <- sqrt((y - P.prop)^2 * rowSums((X %*% W.prop)^2))
  PHI.opt <- sum(PI.opt[idx.prop] * pinv.prop)
  PI.opt <- c(PI.opt / PHI.opt)
  PI.opt
}



Log_Poi_optL <- function(X, y, beta.prop, idx.prop, PI.prop){
  
  pinv.prop <- 1/pmin(r0 * PI.prop[idx.prop], 1)
  
  P.prop  <- 1 - 1 / (1 + exp(c(X %*% beta.prop)))
  p.prop <- P.prop[idx.prop]
  
  PI.opt <- sqrt((y - P.prop)^2 * rowSums(X^2))
  PHI.opt <- sum(PI.opt[idx.prop] * pinv.prop)
  PI.opt <- c(PI.opt / PHI.opt)
  PI.opt
}






