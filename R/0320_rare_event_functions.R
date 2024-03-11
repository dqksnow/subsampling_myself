getMSLE <- function(x,y,pi,pilot) {
  n <- nrow(x)
  d <- ncol(x)
  beta <- pilot # initial value?
  # beta <- rep(0,d)
  S <- matrix(0,n,d)
  H <- matrix(0,d,d)
  loop <- 1
  Loop <- 100
  msg <- "NA"
  while (loop <= Loop) {
    p = 1 / (1 + exp(-as.vector(x %*% beta)) * pi)
    S = (y - p) * x
    phi = p * (1 - p)
    H = t(x) %*% (phi * x)
    ss = colSums(S)
    shs = tryCatch(
      solve(H,ss),
      error = function(e) {
        msg <<- "H is singular"; cat(msg,"\n")
        beta <<- rep(NA,d)
        break
      }
    )
    beta_new = beta + shs * 0.5
    tlr = sum(shs ^ 2)
    beta = beta_new
    if (tlr < .00001) {
      msg <<- "Successful convergence"
      break
    }
    if (loop == Loop) {
      msg <<- "Maximum iteration reached"; cat(msg,"\n")
      beta <<- rep(NA,d)
      break
    }
    loop <- loop + 1
  }

  return(list(beta=beta,msg=msg,
              loop=loop,H=H,SSt=t(S)%*%S))
}
################################################################################
logistic.coef.estimate <-
  function(X, Y, weights = rep(1, length(Y)), offset = rep(0, length(Y)),
                                     method = FALSE) {
  if (method == 'slik' && !missing(weights) && !missing(offset)){
    warning("The 'slik' method is only applicable with the unweighted likelihood function, so 'weights' are set to 1.")
    weights <- rep(1, length(Y))
  }
  if (method == 'popt' && !missing(weights) && !missing(offset)){
    warning("The 'popt' method is only applicable with no offset, so 'offset' are set to 0.")
    offset <- rep(0, length(Y))
  }
  # if (missing(start))
  #   start <- rep(0, ncol(X))
  if (length(weights) == 1)
    weights <- rep(weights, length(Y))
  data <- as.data.frame(cbind(Y, X))
  formula <- paste(colnames(data)[1], "~",
                   paste(colnames(data)[-1], collapse = "+"), "-1")
  design <- survey::svydesign(ids =  ~ 1,
                              weights =  ~ weights,
                              data = data)
  beta <- survey::svyglm(
    as.formula(formula),
    design = design,
    offset = offset,
    family = quasibinomial(link = "logit") # binomial(link=logit)
  )$coefficients
  return(list(beta = beta))
}
###############################################################################
# logistic.coef.estimate <-
#   function(X, Y, weights = rep(1, length(Y)), offset = rep(0, length(Y)),
#            method = FALSE) {
#   if (method == 'slik' && !missing(weights) && !missing(offset)){
#     warning("The 'slik' method is only applicable with the unweighted likelihood function, so 'weights' are set to 1.")
#     weights <- rep(1, length(Y))
#   }
#   if (method == 'popt' && !missing(weights) && !missing(offset)){
#     warning("The 'popt' method is only applicable with no offset, so 'offset' are set to 0.")
#     offset <- rep(0, length(Y))
#   }
#   if (length(weights) == 1)
#     weights <- rep(weights, length(Y))
#   beta <- glm(Y~X-1,
#               family = quasibinomial(link = "logit"),   # binomial(link = "logit")
#               weights = weights,
#               offset = offset)$coefficients
#   return(list(beta = beta))
# }
################################
MN <- function(X, pbeta, pinv){
  phi <- pbeta * (1 - pbeta)
  t(X) %*% (X * phi * pinv)
}
pbeta <- function(X, beta){
  1 / (1 + exp(-c(X %*% beta)))
}
poisson.indx <- function(N, r, pi){
  runif(N) <= (r * pi)
}
################################################################################
pilot.estimate <- function(X, Y, n0){
  N <- nrow(X)
  d <- ncol(X)
  N1 <- sum(Y)
  N0 <- sum(Y == 0)
  p.pilot <- rep(1 / (2 * N1), N)  # case control
  p.pilot[Y == 0] <- 1 / (2 * N0)
  pilot.indx <- poisson.indx(N, n0, p.pilot)
  x.plt <- X[pilot.indx,]
  y.plt <- Y[pilot.indx]
  cc.plt <- n0 * p.pilot[pilot.indx]

  beta.pilot <- logistic.coef.estimate(X = x.plt, Y = y.plt)$beta
  pbeta.pilot <- pbeta(x.plt, beta.pilot)
  MN.pilot <- MN(x.plt, pbeta.pilot, 1)
  beta.pilot[1] <- beta.pilot[1] - log(N0 / N1)
  P.plt <- pbeta(X, beta.pilot)
  dm <-
    P.plt * sqrt(1 - P.plt) * sqrt(rowSums((X %*% solve(MN.pilot)) ^ 2))
  pi.P <- dm / (N0 * sum(dm[pilot.indx] / cc.plt) / N)
  dm.lcc <- abs(Y - P.plt)
  pi.Plcc = dm.lcc / sum(dm.lcc[pilot.indx]/cc.plt)
  return(list(
    beta.pilot = beta.pilot,
    MN.pilot = MN.pilot,
    pi.P = pi.P,
    pi.Plcc=pi.Plcc,
    pilot.indx = pilot.indx
  ))
}
###############
# pilot.estimate <- function(X, Y, n0){
#   N <- nrow(X)
#   d <- ncol(X)
#   N1 <- sum(Y)
#   loc0 <- Y == 0
#   N0 <- sum(loc0)
#   cc <- rep(n0 / (2 * N1), N) # case control
#   cc[Y == 0] <- n0 / (2 * N0)
#   idx.plt <- runif(N) <= cc #pilot sample index
#   x.plt <- X[idx.plt , ]
#   y.plt <- Y[idx.plt]
#   cc.plt <- cc[idx.plt]
#   # getEst.results = getEst(x.plt, y.plt) # calculate likelihood with weight=1
#   # beta.pilot <- getEst.results$beta
#   # ddm.plt <- getEst.results$H
#   beta.pilot <- logistic.coef.estimate(X=x.plt, Y=y.plt)$beta
#   p = 1 / (1 + exp(-as.vector(x.plt %*% beta.pilot)))
#   MN.pilot <-  t(x.plt) %*% (p * (1 - p) * x.plt)
#
#   beta.pilot[1] <- beta.pilot[1] - log(N0/N1) # correction
#   P.plt <- 1 / (1 + exp(-X %*% beta.pilot))
#   dm = P.plt * sqrt(1 - P.plt) * sqrt(rowSums((X %*% solve(MN.pilot))^2)) # ?
#   pi.P <- dm / (N0 * sum(dm[idx.plt]/cc.plt) / N)
#
#   dm.lcc <- abs(Y - P.plt)
#   pi.Plcc = dm.lcc / sum(dm.lcc[idx.plt]/cc.plt)
#   return(
#     list(
#       beta.pilot = beta.pilot,
#       MN.pilot = MN.pilot,
#       pi.P = pi.P,
#       pi.Plcc = pi.Plcc,
#       pilot.indx = idx.plt
#     )
#   )
# }

################################################################################
RareLogistic <- function(X, Y, n0, nss, method = 'slik') {
  N <- nrow(X)
  d <- ncol(X)
  N1 <- sum(Y)
  if (method %in% c("slik", "popt", 'lcc')){
    pilot.estimate.results <- pilot.estimate(X, Y, n0)
    beta.pilot <- pilot.estimate.results$beta.pilot
    # MN.pilot <- pilot.estimate.results$MN.pilot
    pi.P <- pilot.estimate.results$pi.P
    pi.Plcc <- pilot.estimate.results$pi.Plcc
    pilot.index <- pilot.estimate.results$pilot.indx

    idx <- runif(N) <= Y + (1 - Y) * nss * pi.P  # subsample index
    x <- X[idx, ]
    y <- Y[idx]
    p <- pmin(nss * pi.P[idx],1)
    py <- y + (1 - y) * p
    if (method == "popt"){
      beta.est <-
        logistic.coef.estimate(X = x, Y = y, weights = 1 / py)$beta
    } else if (method == "slik"){
      beta.est <-
        logistic.coef.estimate(
          X = x,
          Y = y,
          weights = rep(1, length(y)),
          offset = -log(p)
        )$beta
      # beta.est <- getMSLE(x,y,p,beta.pilot)$beta
    } else if (method == 'lcc'){
      PLCC <- (nss + N1) * pi.Plcc
      idx <- runif(N) <= PLCC
      x.lcc <- X[idx , ]
      y.lcc <- Y[idx]
      w.lcc <- pmax(1 , PLCC)[idx]
      beta.est <-
        logistic.coef.estimate(X = x.lcc, Y = y.lcc, weights = w.lcc)$beta + beta.pilot
    }
    return(
      list(
        beta.pilot = beta.pilot,
        beta.est = beta.est,
        pilot.index = which(pilot.index == T),
        subsample.idx = which(idx == T)
      )
    )
  }

  if (method %in% c("Uni", "UniW")){
    loc0 = Y == 0
    N0 <- sum(loc0)
    Betas <- rep(NA, d)
    pi.uni <- rep(1,N)
    pi.uni[loc0] <- nss/N0
    idx.uni <- runif(N) <= pi.uni
    x.uni <- X[idx.uni, ]
    y.uni <- Y[idx.uni]
    if (method == "Uni"){
      Betas <- logistic.coef.estimate(X = x.uni, Y = y.uni)$beta
      n.star <- sum(idx.uni[loc0])
      Betas[1] <- Betas[1] + log(nss / N0)
    } else if (method == "UniW"){
      w.uni <- 1 / pi.uni[idx.uni]
      Betas <-
        logistic.coef.estimate(X = x.uni, Y = y.uni, weights = w.uni)$beta
    }
    return(list(beta.est = Betas))
  }
}
