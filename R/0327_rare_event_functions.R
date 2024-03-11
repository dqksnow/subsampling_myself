################################################################################
logistic.coef.estimate <-
  function(X, Y, weights = rep(1, length(Y)), offset = rep(0, length(Y)), #order
           method = FALSE) {
# stop()
    if (method == 'LogOddsCorrection' && !missing(weights) && !missing(offset)){
      warning("The 'LogOddsCorrection' method is only applicable with the unweighted likelihood function, so 'weights' are set to 1.")
      weights <- rep(1, length(Y))
    }
    if (method == 'Weighted' && !missing(weights) && !missing(offset)){
      warning("The 'Weighted' method is only applicable with no offset, so 'offset' are set to 0.")
      offset <- rep(0, length(Y))
    }
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
      family = quasibinomial(link = "logit")
    )$coefficients
    return(list(beta = beta))
  }
###############################################################################
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
pilot.estimate <- function(X, Y, n0, criterion = 'optA'){
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
  if (criterion == 'optA'){
    dm <- P.plt * sqrt(1 - P.plt) * sqrt(rowSums((X %*% solve(MN.pilot)) ^ 2))
  } else if (criterion == 'optL'){
    dm <- P.plt * sqrt(1 - P.plt) * sqrt(rowSums(X ^ 2))
  }
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
################################################################################
RareLogistic <- function(X, Y, n0, nss, method = 'LogOddsCorrection', criterion = 'optL') { #criterion
  N <- nrow(X)
  d <- ncol(X)
  N1 <- sum(Y)
  if (method %in% c("LogOddsCorrection", "Weighted", 'LCC')){
    if (criterion == 'optA'){
      pilot.estimate.results <- pilot.estimate(X, Y, n0, criterion = 'optA')
    } else if (criterion == 'optL'){
      pilot.estimate.results <- pilot.estimate(X, Y, n0, criterion = 'optL')
    }
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
    if (method == "Weighted"){
      beta.est <-
        logistic.coef.estimate(X = x, Y = y, weights = 1 / py)$beta
    } else if (method == "LogOddsCorrection"){
      beta.est <-
        logistic.coef.estimate(
          X = x,
          Y = y, #start value
          weights = rep(1, length(y)),
          offset = -log(p)
        )$beta
      # beta.est <- getMSLE(x,y,p,beta.pilot)$beta
    } else if (method == 'LCC'){
      PLCC <- (nss + N1) * pi.Plcc
      idx <- runif(N) <= PLCC # paper
      x.lcc <- X[idx , ]
      y.lcc <- Y[idx]
      w.lcc <- pmax(1 , PLCC)[idx]
      beta.est <-
        logistic.coef.estimate(X = x.lcc, Y = y.lcc, weights = w.lcc)$beta + beta.pilot
    }
    return(
      list(
        beta.polot = beta.pilot, # combine
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
