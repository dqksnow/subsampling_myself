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
poisson.index <- function(N, r, pi){
  runif(N) <= (r * pi)
}
################################################################################
pilot.estimate <- function(X, Y, n0){
  N <- nrow(X)
  d <- ncol(X)
  N1 <- sum(Y)
  N0 <- N - N1
  p.pilot <- rep(1 / (2 * N1), N)  # case control
  p.pilot[Y == 0] <- 1 / (2 * N0)
  index.pilot <- poisson.index(N, n0, p.pilot)
  x.pilot <- X[index.pilot,]
  y.pilot <- Y[index.pilot]
  cc.pilot <- n0 * p.pilot[index.pilot]

  beta.pilot <- logistic.coef.estimate(X = x.pilot, Y = y.pilot)$beta
  pbeta.pilot <- pbeta(x.pilot, beta.pilot)
  MN.pilot <- MN(x.pilot, pbeta.pilot, 1)
  beta.pilot[1] <- beta.pilot[1] - log(N0 / N1)
  P.pilot <- pbeta(X, beta.pilot)
  return(list(
    beta.pilot = beta.pilot,
    MN.pilot = MN.pilot,
    P.pilot = P.pilot,
    index.pilot = index.pilot,
    cc.pilot = cc.pilot
  ))
}
subsampling <- function(X, Y, nss, criterion,
                        P.pilot, MN.pilot, index.pilot, cc.pilot
                        ){
  N <- nrow(X)
  N1 <- sum(Y)
  N0 <- N - N1

  if (criterion == 'optA'){
    dm <- P.pilot * sqrt(1 - P.pilot) * sqrt(rowSums((X %*% solve(MN.pilot)) ^ 2))
    pi <- dm / (N0 * sum(dm[index.pilot] / cc.pilot) / N)
  } else if (criterion == 'optL'){
    dm <- P.pilot * sqrt(1 - P.pilot) * sqrt(rowSums(X ^ 2))
    pi <- dm / (N0 * sum(dm[index.pilot] / cc.pilot) / N)
  } else if (criterion == 'LCC'){
    dm <- abs(Y - P.pilot)
    pi <- dm / sum(dm[index.pilot] / cc.pilot)
  }

  if(criterion == 'LCC'){
    PLCC <- (nss + N1) * pi
    index.subsample <- runif(N) <= PLCC # paper
    w.subsample <- pmax(1 , PLCC)[index.subsample]
    return (list(index.subsample = index.subsample, w.subsample = w.subsample))
  } else if(criterion %in% c('optL', 'optA')){
    index.subsample <- runif(N) <= Y + (1 - Y) * nss * pi  # subsample index
    y.subsample <- Y[index.subsample]
    p.subsample <- pmin(nss * pi[index.subsample], 1)
    w.subsample <- y.subsample + (1 - y.subsample) * p.subsample
    return (
      list(
        index.subsample = index.subsample,
        w.subsample = w.subsample,
        p.subsample = p.subsample
      )
    )
  }

}

################################################################################
RareLogistic <- function(X, Y, n0, nss, method = 'LogOddsCorrection', criterion = 'optL') { #criterion
  N <- nrow(X)
  d <- ncol(X)
  N1 <- sum(Y)
  N0 <- N - N1
  loc0 = Y == 0

  # Uni: uniform sampling + correction
  # UniW: uniform sampling + weighted likelihood
  if (method %in% c("Uni", "UniW")){
    beta.est <- rep(NA, d)
    pi.uni <- rep(1,N)
    pi.uni[loc0] <- nss/N0
    index.uni <- runif(N) <= pi.uni
    x.uni <- X[index.uni, ]
    y.uni <- Y[index.uni]
    if (method == "Uni"){
      beta.est <- logistic.coef.estimate(X = x.uni, Y = y.uni)$beta
      # n.star <- sum(index.uni[loc0])
      beta.est[1] <- beta.est[1] + log(nss / N0)
    } else if (method == "UniW"){
      beta.est <-
        logistic.coef.estimate(X = x.uni, Y = y.uni, weights = 1 / pi.uni[index.uni])$beta
    }
    return(list(index.subsample = which(index.uni == TRUE),
                beta.est = beta.est))
  } else if (method %in% c("Weighted", "LogOddsCorrection")){
    # pilot step
    pilot.estimate.results <- pilot.estimate(X, Y, n0)
    beta.pilot <- pilot.estimate.results$beta.pilot
    MN.pilot <- pilot.estimate.results$MN.pilot
    P.pilot <- pilot.estimate.results$P.pilot
    index.pilot <- pilot.estimate.results$index.pilot
    cc.pilot <- pilot.estimate.results$cc.pilot

    # subsampling step
    subsampling.results <- subsampling(X, Y, nss, criterion,
                                       P.pilot, MN.pilot, index.pilot, cc.pilot)
    index.subsample <- subsampling.results$index.subsample
    w.subsample <- subsampling.results$w.subsample
    x.subsample <- X[index.subsample, ]
    y.subsample <- Y[index.subsample]

    # subsample estimate step
    if (criterion == 'LCC'){
      beta.est <-
        logistic.coef.estimate(X = x.subsample, Y = y.subsample, weights = w.subsample)$beta + beta.pilot
    } else if (method == "Weighted"){
      beta.est <-
        logistic.coef.estimate(X = x.subsample, Y = y.subsample, weights = 1 / w.subsample)$beta
    } else if (method == "LogOddsCorrection"){
      p.subsample <- subsampling.results$p.subsample
      beta.est <-
        logistic.coef.estimate(
          X = x.subsample,
          Y = y.subsample, #start balue
          weights = rep(1, length(y.subsample)),
          offset = -log(p.subsample)
        )$beta
      # beta.est <- getMSLE(x,y,p,beta.pilot)$beta
    }
    return(
      list(
        beta.polot = beta.pilot, # combine
        beta.est = beta.est,
        index.pilot = which(index.pilot == TRUE),
        index.subsample = which(index.subsample == TRUE)
      )
    )
    }


}

