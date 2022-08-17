
#' Title
#'
#' @param X Covariate
#' @param y Response
#' @param r Sample size
#' @param intercept Wheter contains intercept for regression
#'
#' @return
#' @export
#'
#' @examples
uniform.Samp <- function(X, y, r, intercept = F){
  N <- length(y)
  d <- dim(X)[2]
  cva <- ifelse(is.null(colnames(X)),
                paste(paste0("V", 1:d), collapse = "+"),
                paste(colnames(X), collapse = "+"))
  rsp <- ifelse(is.null(names(y)), "y", names(y))
  fmla <- ifelse(intercept,
                 paste0(rsp, "~", cva),
                 paste0(rsp, "~", cva, "-1"))

  idx.unif <- sample(1:N, r, T)
  x.unif <- X[idx.unif,]
  y.unif <- y[idx.unif]
  p.unif <- rep(1/N, r)
  pinv.unif <- 1/p.unif

  design <- survey::svydesign(ids=~1,
                              weights=~pinv.unif,
                              data=as.data.frame(cbind(X, y)[idx.unif,]))
  beta.unif <- survey::svyglm(as.formula(fmla),
                              design=design,
                              family=quasibinomial)$coefficients

  return(list(idx.smpl = idx.unif,
              p.smpl = p.unif,
              beta.smpl = beta.unif))
}


#' Title
#'
#' @param X Covariate
#' @param y Response
#' @param r Sample size
#' @param intercept Wheter contains intercept for regression
#'
#' @return
#' @export
#'
#' @examples
prop.Samp <- function(X, y, r, intercept = F){
  N <- length(y)
  d <- dim(X)[2]
  cva <- ifelse(is.null(colnames(X)),
                paste(paste0("V", 1:d), collapse = "+"),
                paste(colnames(X), collapse = "+"))
  rsp <- ifelse(is.null(names(y)), "y", names(y))
  fmla <- ifelse(intercept,
                 paste0(rsp, "~", cva),
                 paste0(rsp, "~", cva, "-1"))

  n1 <- sum(y)
  n0 <- N - n1
  PI.prop <- rep(1/(2*n1), N)
  PI.prop[y==0] <- 1/(2*n0)
  idx.prop <- sample(1:N, r, T, PI.prop)
  x.prop <- X[idx.prop,]
  y.prop <- y[idx.prop]
  p.prop <- rep(1/N, r)
  pinv.prop <- 1/p.prop

  design <- survey::svydesign(ids=~1,
                              weights=~pinv.prop,
                              data=as.data.frame(cbind(X, y)[idx.prop,]))
  beta.prop <- survey::svyglm(as.formula(fmla),
                              design=design,
                              family=quasibinomial)$coefficients

  return(list(idx.smpl = idx.prop,
              p.smpl = p.prop,
              beta.smpl = beta.prop))
}

