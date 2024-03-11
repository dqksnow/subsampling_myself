getMSLE <- function(X, Y, offset = NA, start) {
  n <- nrow(X)
  d <- ncol(X)
  beta <- start
  S <- matrix(0, n, d)
  H <- matrix(0, d, d)
  loop <- 1
  Loop <- 100
  step.length <- 0.5
  msg <- "NA"
  while (loop <= Loop) {
    ifelse(is.na(offset),
           p <- 1 / (1 + exp(-as.vector(X %*% beta))),
           p <- 1 / (1 + exp(-as.vector(X %*% beta) - offset)))
    S <- (Y - p) * X
    phi <- p * (1 - p)
    MN <- t(X) %*% (phi * X)
    ss <- colSums(S)
    shs <- tryCatch(
      solve(MN, ss),
      error = function(e) {
        msg <- "MN is singular"; cat(msg,"\n")
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
  return(list(beta = beta, MN = MN, Psi = t(S) %*% S)) # not used
}
###############################################################################
logistic.coef.estimate <- function(X,
                                   Y,
                                   offset = rep(0, length(Y)),
                                   start = rep(0, ncol(X)),
                                   weights = rep(1, length(Y))) {
  # if (length(weights) == 1)
  #   weights <- rep(weights, length(Y))
  # if (offset != rep(0, length(Y)) && weights != rep(1, length(Y))) {
  #   stop("offset will only be needed in the unweighted likelihood function.")
  # }
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
uni.index <- function(N, Y, n.plt){ # same
  N1 <- sum(Y)
  N0 <- N - N1
  if (N0 < n.plt/2 | N1 < n.plt/2){
    warning(paste("n.plt/2 exceeds the number of Y=1 or Y=0 in the full data.",
                  "All rare events will be drawn into the pilot sample."))
  }
  n.plt.0 <- min(N0, n.plt/2)
  n.plt.1 <- min(N1, n.plt/2)
  index.plt <- c(sample(which(Y == 0), n.plt.0, replace = FALSE),
                 sample(which(Y == 1), n.plt.1, replace = FALSE))
  return (list(index.plt = index.plt, n.plt.0 = n.plt.0, n.plt.1 = n.plt.1))
}
poisson.index <- function(N, pi){
  index <- runif(N) <= pi
  return(which(index == TRUE))
}
MN <- function(X, pbeta, w){
  phi <- pbeta * (1 - pbeta)
  t(X) %*% (X * (phi * w))
}
Psi <- function(X, Y, pbeta, w){
  psi <- (Y - pbeta)^2
  t(X) %*% (X * (psi * w^2))
}
pbeta <- function(X, beta, offset = NA){
  ifelse(is.na(offset),
         p <- 1 / (1 + exp(-as.vector(X %*% beta))),
         p <- 1 / (1 + exp(-as.vector(X %*% beta) - offset)))
  return (p)
}
###############################################################################
pilot.estimate <- function(X, Y, n.plt){
  N <- nrow(X)
  d <- ncol(X)
  N1 <- sum(Y)
  N0 <- N - N1

  uni.index.results <- uni.index(N, Y, n.plt)
  index.plt <- uni.index.results$index.plt
  n.plt.0 <- uni.index.results$n.plt.0
  n.plt.1 <- uni.index.results$n.plt.1
  x.plt <- X[index.plt,]
  y.plt <- Y[index.plt]

  beta.plt <- logistic.coef.estimate(X = x.plt, Y = y.plt)$beta
  pbeta.plt <- pbeta(x.plt, beta.plt)
  MN.plt <- MN(x.plt, pbeta.plt, 1/n.plt) # use to be 1
  Psi.plt <- Psi(x.plt, y.plt, pbeta.plt, 1)
  beta.plt[1] <- beta.plt[1] - log(N0 / N1)
  P.plt <- pbeta(X, beta.plt)

  return(list(
    beta.plt = beta.plt,
    MN.plt = MN.plt,
    Psi.plt = Psi.plt,
    P.plt = P.plt,
    index.plt = index.plt
  ))
}
###############################################################################
subsampling <- function(X, Y, n.plt, n.ssp, criterion,
                        P.plt, MN.plt, index.plt, cc.plt){
  N <- nrow(X)
  N1 <- sum(Y)
  N0 <- N - N1

  cc.plt <- n.plt * rep(1 / (2 * N1), length(index.plt))
  cc.plt[which(Y[index.plt] == 0)] <- n.plt * 1 / (2 * N0)

  if (criterion == 'optA'){
    nm <- P.plt * sqrt(1 - P.plt) * sqrt(colSums(solve(MN.plt, t(X))^2)) #
    dm <- (N0 * sum(nm[index.plt] / cc.plt) / N) # Q： why * N0 / N ?
    # print(sum(nm))
    # print(dm)
    pi <- nm / dm
    # print(paste('sum(pi): ', sum(pi))) # close to 1
  } else if (criterion == 'optL'){
    nm <- P.plt * sqrt(1 - P.plt) * sqrt(rowSums(X ^ 2)) #
    dm <- (N0 * sum(nm[index.plt] / cc.plt) / N)
    # print(sum(nm))
    # print(dm)
    pi <- nm / dm
    # print(paste('sum(pi): ', sum(pi)))
  } else if (criterion == 'LCC'){
    nm <- abs(Y - P.plt)
    dm <- sum(nm[index.plt] / cc.plt)
    pi <- nm / dm
    # print(paste('sum(pi): ', sum(pi)))
  }
  if(criterion == 'LCC'){
    P.lcc <- (n.ssp + N1) * pi
    index.ssp <- poisson.index(N, P.lcc)
    p.ssp <- pmin(P.lcc[index.ssp], 1)
    w.ssp <- 1 / p.ssp
    # offset:
    nm.1 <- abs(1 - P.plt[index.ssp])
    nm.0 <- abs(P.plt[index.ssp])
    pi.1 <- pmin( (n.ssp + N1) * nm.1 / dm, 1)
    pi.0 <- pmin( (n.ssp + N1) * nm.0 / dm, 1)
    offset <- log(pi.1 / pi.0)
  } else if(criterion %in% c('optL', 'optA')){
    p.ssp <- n.ssp * pi
    # plot(p.ssp)
    # print(sum(p.ssp))
    # print(sum(pmin(p.ssp, 1)))
    index.ssp <- poisson.index(N, Y + (1 - Y) * p.ssp)
    # weights and offset:
    p.ssp <- pmin(p.ssp[index.ssp], 1)
    w.ssp <- 1 / (Y[index.ssp] + (1 - Y[index.ssp]) * p.ssp)
    offset <- -log(p.ssp)
    # print(offset)
  }
  return (
    list(
      index.ssp = index.ssp,
      w.ssp = w.ssp,
      offset = offset,
      p.ssp = n.ssp * pi
    )
  )
}
###############################################################################
subsample.estimate <- function(x.ssp, y.ssp, n.ssp, w.ssp, offset, beta.plt, estimate.method){
  if (estimate.method == "Weighted"){
    beta.ssp <- logistic.coef.estimate(X = x.ssp, Y = y.ssp, weights = w.ssp)$beta #should start = beta.plt?
    P.ssp  <- pbeta(x.ssp, beta.ssp)
    MN.ssp <- MN(x.ssp, P.ssp, w.ssp/n.ssp)
    Psi.ssp <- Psi(x.ssp, y.ssp, P.ssp, w.ssp)
  } else if (estimate.method == 'LogOddsCorrection'){
    beta.ssp <-
      logistic.coef.estimate(
        X = x.ssp,
        Y = y.ssp,
        start = beta.plt,
        offset = offset
      )$beta
    P.ssp  <- pbeta(x.ssp, beta.ssp, offset)
    MN.ssp <- MN(x.ssp, P.ssp, 1/n.ssp)
    Psi.ssp <- Psi(x.ssp, y.ssp, P.ssp, 1)
  }
  return (
    list(
      beta.ssp = beta.ssp,
      MN.ssp = MN.ssp,
      Psi.ssp = Psi.ssp
    )
  )
}
###############################################################################
RareLogistic <- function(X, Y, n.plt, n.ssp, estimate.method = 'LogOddsCorrection', criterion = 'optL') {
  N <- nrow(X)
  d <- ncol(X)
  N1 <- sum(Y)
  N0 <- N - N1

  # if (!all(Y %in% c(0, 1))) {
  #   message("Y contains values other than 0 and 1. Please check Y.")
  # }

  if (estimate.method %in% c("Weighted", "LogOddsCorrection")){
    # pilot step
    plt.estimate.results <- pilot.estimate(X = X, Y = Y, n.plt = n.plt)
    beta.plt <- plt.estimate.results$beta.plt
    MN.plt <- plt.estimate.results$MN.plt
    Psi.plt <- plt.estimate.results$Psi.plt
    P.plt <- plt.estimate.results$P.plt
    index.plt <- plt.estimate.results$index.plt
    cc.plt <- plt.estimate.results$cc.plt
    # subsampling step
    ssp.results <- subsampling(X, Y, n.plt, n.ssp, criterion = criterion,
                               P.plt, MN.plt, index.plt, cc.plt)
    index.ssp <- ssp.results$index.ssp
    w.ssp <- ssp.results$w.ssp
    offset <- ssp.results$offset
    p.ssp <- ssp.results$p.ssp
    # subsample estimation step
    ssp.estimate.results <-
      subsample.estimate(X[index.ssp, ], Y[index.ssp], n.ssp, w.ssp, offset,
                         beta.plt, estimate.method = estimate.method)
    beta.ssp <- ssp.estimate.results$beta.ssp
    MN.ssp <- ssp.estimate.results$MN.ssp
    Psi.ssp <- ssp.estimate.results$Psi.ssp

    # combine step
    # print(MN.plt)
    # print(MN.ssp)
    MN.plt <- MN.plt * n.plt
    MN.ssp <- MN.ssp * n.ssp
    # print(MN.plt)
    # print(MN.ssp)
    MNsolve <- solve(MN.plt + MN.ssp)
    beta.cmb <- c(MNsolve %*% (MN.plt %*% beta.plt + MN.ssp %*% beta.ssp))
    var.cmb <- MNsolve %*% (Psi.plt + Psi.ssp) %*% MNsolve

    return(list(beta.plt = beta.plt,
                beta.ssp = beta.ssp,
                beta.cmb = beta.cmb,
                var.cmb = var.cmb,
                index.plt = index.plt,
                index.ssp = index.ssp,
                p.ssp = p.ssp
    )
    )
  } else if (estimate.method %in% c("Uni", "UniW")){
    # Uni:  poisson sampling + unweighted likelihood estimation + pilot correction
    # UniW: poisson sampling + weighted likelihood estimation
    loc0 = Y == 0
    beta.ssp <- rep(NA, d)
    pi.uni <- rep(1,N)
    pi.uni[loc0] <- n.ssp / N0
    index.uni <- runif(N) <= pi.uni
    x.uni <- X[index.uni, ]
    y.uni <- Y[index.uni]
    if (estimate.method == "Uni"){
      beta.ssp <- logistic.coef.estimate(X = x.uni, Y = y.uni)$beta
      # n.star <- sum(index.uni[loc0])
      beta.ssp[1] <- beta.ssp[1] + log(n.ssp / N0)
    } else if (estimate.method == "UniW"){
      beta.ssp <-
        logistic.coef.estimate(X = x.uni, Y = y.uni, weights = 1 / pi.uni[index.uni])$beta
    }
    return(list(index.ssp = which(index.uni == TRUE),
                beta.ssp = beta.ssp))
  }
}

