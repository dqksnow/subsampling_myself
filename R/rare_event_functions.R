getMSLE <- function(x, y, offset, start) {
  n <- nrow(x)
  d <- ncol(x)
  beta <- start
  S <- matrix(0, n, d)
  H <- matrix(0, d, d)
  loop <- 1
  Loop <- 100
  step.length <- 0.5
  msg <- "NA"
  while (loop <= Loop) {
    p <- 1 / (1 + exp(-as.vector(x %*% beta) - offset))
    # p <- 1 / (1 + exp(-as.vector(x %*% beta)) * offset)
    S <- (y - p) * x
    phi <- p * (1 - p)
    H <- t(x) %*% (phi * x)
    ss <- colSums(S)
    shs <- tryCatch(
      solve(H, ss),
      error = function(e) {
        msg <- "H is singular"; cat(msg,"\n")
        beta <- rep(NA, d)
        break
      }
    )
    beta.new <- beta + shs * step.length
    tlr <- sum(shs ^ 2)
    beta <- beta.new
    if (tlr < .00001) {
      msg <- "Successful convergence"
      break
    }
    if (loop == Loop) {
      msg <- "Maximum iteration reached"; cat(msg, "\n")
      beta <- rep(NA, d)
      break
    }
    loop <- loop + 1
  }
  return(list(beta = beta))
}
###############################################################################
logistic.coef.estimate <-
  function(X,
           Y,
           offset = rep(0, length(Y)),
           start = rep(0, ncol(X)),
           weights = rep(1, length(Y))) {
    if (length(weights) == 1)
      weights <- rep(weights, length(Y))
    if (offset != rep(0, length(Y)) & weights != rep(1, length(Y))) {
      stop("offset will only be needed in the unweighted likelihood function.")
    }
    data <- as.data.frame(cbind(Y, X))
    formula <- paste(colnames(data)[1], "~",
                     paste(colnames(data)[-1], collapse = "+"), "-1")
    design <- survey::svydesign(ids =  ~ 1,
                                weights =  ~ weights,
                                data = data)
    fit <- try(
      beta <- survey::svyglm( as.formula(formula),
                              design = design,
                              offset = offset,
                              start = start,
                              family = quasibinomial(link = "logit"))$coefficients,
      silent=TRUE
    )
    if ("try-error" %in% class(fit)){
      message("Warning: an error occurred while calling 'glm.fit': ", geterrmessage(),
              "This is probably due to the iteratively reweighted least squares algorithm called by 'glm.fit' not getting converged. Tring another function 'getMSLE' to replace 'glm.fit'.")
      beta <- getMSLE(X, Y, offset = offset, start = start)$beta
    }
    return(list(beta = beta))
  }
###############################################################################
MN <- function(X, pbeta, pinv){
  phi <- pbeta * (1 - pbeta)
  t(X) %*% (X * phi * pinv)
}
Psi <- function(X, Y, pbeta, pinv){
  psi <- (Y - pbeta)^2
  t(X) %*% (X * psi * pinv^2)
}
pbeta <- function(X, beta){
  1 / (1 + exp(-c(X %*% beta)))
}
poisson.index <- function(N, r, pi){
  runif(N) <= (r * pi)
}
###############################################################################
pilot.estimate <- function(X, Y, n0){
  N <- nrow(X)
  d <- ncol(X)
  N1 <- sum(Y)
  N0 <- N - N1
  p.plt <- rep(1 / (2 * N1), N)  # case control
  p.plt[Y == 0] <- 1 / (2 * N0)
  index.plt <- poisson.index(N, n0, p.plt)
  x.plt <- X[index.plt,]
  y.plt <- Y[index.plt]
  cc.plt <- n0 * p.plt[index.plt]

  beta.plt <- logistic.coef.estimate(X = x.plt, Y = y.plt)$beta
  pbeta.plt <- pbeta(x.plt, beta.plt)
  MN.plt <- MN(x.plt, pbeta.plt, 1)
  beta.plt[1] <- beta.plt[1] - log(N0 / N1)
  P.plt <- pbeta(X, beta.plt)
  Psi.plt <- Psi(x.plt, y.plt, pbeta.plt, 1)
  return(list(
    beta.plt = beta.plt,
    MN.plt = MN.plt,
    Psi.plt = Psi.plt,
    P.plt = P.plt,
    index.plt = index.plt,
    cc.plt = cc.plt
  ))
}
###############################################################################
subsampling <- function(X, Y, nss, criterion,
                        P.plt, MN.plt, index.plt, cc.plt
){
  N <- nrow(X)
  N1 <- sum(Y)
  N0 <- N - N1

  if (criterion == 'optA'){
    dm <- P.plt * sqrt(1 - P.plt) * sqrt(rowSums((X %*% solve(MN.plt)) ^ 2))
    pi <- dm / (N0 * sum(dm[index.plt] / cc.plt) / N)
  } else if (criterion == 'optL'){
    dm <- P.plt * sqrt(1 - P.plt) * sqrt(rowSums(X ^ 2))
    pi <- dm / (N0 * sum(dm[index.plt] / cc.plt) / N)
  } else if (criterion == 'LCC'){
    dm <- abs(Y - P.plt)
    pi <- dm / sum(dm[index.plt] / cc.plt)
  }

  if(criterion == 'LCC'){
    PLCC <- (nss + N1) * pi
    index.ssp <- runif(N) <= PLCC # paper
    w.ssp <- pmax(PLCC[index.ssp], 1)
    return (list(index.ssp = index.ssp, w.ssp = w.ssp))
  } else if(criterion %in% c('optL', 'optA')){
    P0 <- nss * pi
    index.ssp <- runif(N) <= Y + (1 - Y) * P0
    y.ssp <- Y[index.ssp]
    p.ssp <- pmin(P0[index.ssp], 1)
    w.ssp <- y.ssp + (1 - y.ssp) * p.ssp
    return (
      list(
        index.ssp = index.ssp,
        w.ssp = w.ssp,
        p.ssp = p.ssp
      )
    )
  }
}
###############################################################################
subsample.estimate <- function(x.ssp, y.ssp, w.ssp, p.ssp, beta.plt, method, criterion){
  if (criterion == 'LCC'){
    beta.ssp <-
      logistic.coef.estimate(X = x.ssp, Y = y.ssp, weights = w.ssp)$beta + beta.plt
  } else if (method == "Weighted"){
    beta.ssp <- logistic.coef.estimate(X = x.ssp, Y = y.ssp, start = beta.plt, weights = 1 / w.ssp)$beta
  } else if (method == 'LogOddsCorrection'){
    # beta.ssp <-  getMSLE(x.ssp, y.ssp, p.ssp, beta.plt)$beta
    beta.ssp <-
      logistic.coef.estimate(
        X = x.ssp,
        Y = y.ssp,
        start = beta.plt,
        weights = rep(1, length(y.ssp)),
        offset = -log(p.ssp)
      )$beta
  }
  P.ssp  <- pbeta(x.ssp, beta.ssp)
  MN.ssp <- MN(x.ssp, P.ssp, 1 / w.ssp) #1/w.ssp
  Psi.ssp <- Psi(x.ssp, y.ssp, P.ssp, 1 / w.ssp)
  return (
    list(
      beta.ssp = beta.ssp,
      MN.ssp = MN.ssp,
      Psi.ssp = Psi.ssp
    )
  )
}
###############################################################################
RareLogistic <- function(X, Y, n0, nss, method = 'LogOddsCorrection', criterion = 'optL') {
  N <- nrow(X)
  d <- ncol(X)
  N1 <- sum(Y)
  N0 <- N - N1
  loc0 = Y == 0
  if (method %in% c("Weighted", "LogOddsCorrection")){
    # pilot step
    plt.estimate.results <- pilot.estimate(X = X, Y = Y, n0 = n0)
    beta.plt <- plt.estimate.results$beta.plt
    MN.plt <- plt.estimate.results$MN.plt
    Psi.plt <- plt.estimate.results$Psi.plt
    P.plt <- plt.estimate.results$P.plt
    index.plt <- plt.estimate.results$index.plt
    cc.plt <- plt.estimate.results$cc.plt

    # subsampling step
    ssp.results <- subsampling(X, Y, nss, criterion = criterion,
                               P.plt, MN.plt, index.plt, cc.plt)
    index.ssp <- ssp.results$index.ssp
    w.ssp <- ssp.results$w.ssp
    p.ssp <- ssp.results$p.ssp
    x.ssp <- X[index.ssp, ]
    y.ssp <- Y[index.ssp]

    # subsample estimation step
    ssp.estimate.reults <-
      subsample.estimate(x.ssp, y.ssp, w.ssp, p.ssp,
                         beta.plt, method = method, criterion = criterion)
    beta.ssp <- ssp.estimate.reults$beta.ssp
    MN.ssp <- ssp.estimate.reults$MN.ssp
    Psi.ssp <- ssp.estimate.reults$Psi.ssp

    # combine step
    MNsolve <- solve(MN.plt + MN.ssp)
    beta.cmb <- c(MNsolve %*% (MN.plt %*% beta.plt + MN.ssp %*% beta.ssp))
    var.beta <- MNsolve %*% (Psi.plt + Psi.ssp) %*% MNsolve

    return(
      list(
        beta.ssp = beta.ssp,
        beta.cmb = beta.cmb,
        var.beta = var.beta,
        index.plt = which(index.plt == TRUE),
        index.ssp = which(index.ssp == TRUE)
      )
    )
  } else if (method %in% c("Uni", "UniW")){
    # Uni:  poisson sampling + unweighted likelihood estimation + pilot correction
    # UniW: poisson sampling + weighted likelihood estimation
    beta.ssp <- rep(NA, d)
    pi.uni <- rep(1,N)
    pi.uni[loc0] <- nss/N0
    index.uni <- runif(N) <= pi.uni
    x.uni <- X[index.uni, ]
    y.uni <- Y[index.uni]
    if (method == "Uni"){
      beta.ssp <- logistic.coef.estimate(X = x.uni, Y = y.uni)$beta
      # n.star <- sum(index.uni[loc0])
      beta.ssp[1] <- beta.ssp[1] + log(nss / N0)
    } else if (method == "UniW"){
      beta.ssp <-
        logistic.coef.estimate(X = x.uni, Y = y.uni, weights = 1 / pi.uni[index.uni])$beta
    }
    return(list(index.ssp = which(index.uni == TRUE),
                beta.ssp = beta.ssp))
  }
}

