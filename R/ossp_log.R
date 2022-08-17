#' Title
#'
#' @param X Covariate matrix
#' @param y Response vector
#' @param beta.prop
#' @param idx.prop
#' @param PI.prop
#'
#' @return
#' @export
#'
#' @examples
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

