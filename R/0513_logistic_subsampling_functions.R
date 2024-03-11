getMSLE <- function(X, Y, offset, beta.start) {
  n <- nrow(X)
  d <- ncol(X)
  beta <- beta.start
  S <- matrix(0, n, d)
  H <- matrix(0, d, d)
  loop <- 1
  Loop <- 100
  step.length <- 0.5
  msg <- "NA"
  while (loop <= Loop) {
    ifelse(is.null(offset),
           p <- 1 / (1 + exp(-as.vector(X %*% beta))),
           p <- 1 / (1 + exp(-as.vector(X %*% beta) - offset)))

    # p <- 1 / (1 + exp(-as.vector(X %*% beta) - offset))
    S <- (Y - p) * X
    phi <- p * (1 - p)
    H <- t(X) %*% (phi * X)
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
logistic.coef.estimate <- function(X,
                                   Y,
                                   offset = NULL, #  rep(0, length(Y))
                                   beta.start = rep(0, ncol(X)),
                                   weights = 1) {
  if (!is.null(offset) && weights != 1) {
    stop("'offset' will only be needed in the unweighted method 'LogOddsCorrection'.")
  }
  data <- as.data.frame(cbind(Y, X))
  formula <- paste(colnames(data)[1], "~",
                   paste(colnames(data)[-1], collapse = "+"), "-1")
  design <- survey::svydesign(ids =  ~ 1,
                              weights = ~ weights,
                              data = data)
  fit <- try(
    beta <- do.call(survey::svyglm, list(formula = as.formula(formula),
                                         design = design,
                                         start = beta.start,
                                         family = quasibinomial(link = "logit"),
                                         if(!is.null(offset)) offset = offset
    ))$coefficients,
    silent = TRUE
  )
  if ("try-error" %in% class(fit)){
    message("Warning: an error occurred while calling 'svyglm': ",
            geterrmessage(),
            paste("This is probably due to the iteratively reweighted least",
                  "squares algorithm called by 'svyglm' not getting converged.",
                  "Tring another function 'getMSLE' to replace 'svyglm'."))
    beta <- getMSLE(X, Y, if(!is.null(offset)) offset = offset, beta.start = beta.start)$beta
  }
  return(list(beta = beta))
}
###############################################################################
uni.index <- function(N, Y, n.plt){
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
MN <- function(X, pbeta, pinv){
  phi <- pbeta * (1 - pbeta)
  # t(X) %*% diag(phi * pinv) %*% X
  t(X) %*% (X * phi * pinv)
}
Psi <- function(X, Y, pbeta, pinv){
  psi <- (Y - pbeta)^2
  # t(X) %*% diag(psi * pinv^2) %*% X
  t(X) %*% (X * psi * pinv^2)
}
pbeta <- function(X, beta){
  1 / (1 + exp(-c(X %*% beta)))
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
  pbeta.plt <- pbeta(x.plt, beta.plt) # used to calculate derivative
  MN.plt <- MN(x.plt, pbeta.plt, 1)
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
subsampling <- function(X, Y, n.ssp, alpha, b, criterion, sampling.method,
                        MN.plt, P.plt, index.plt){

  N <- nrow(X)
  N1 <- sum(Y)
  N0 <- N - N1
  n.plt <- length(index.plt) # length(index.plt) might be smaller than n.plt you set.

  # criterion = optA, optL, LCC
  # sampling.method = WithReplacement, Poisson
  # estimate.method = Weighted, LogOddsCorrection

  if (criterion == "optA"){
    norm <- sqrt(colSums(solve(MN.plt, t(X))^2))
    nm <- abs(Y - P.plt) * norm # numerator
    nm.1 <- abs(1 - P.plt) * norm # Q: use X[index.plt, ] in norm?
    nm.0 <- abs(P.plt) * norm
  } else if (criterion == "optL"){
    norm <- sqrt(rowSums(X^2))
    nm <- abs(Y - P.plt) * norm
    nm.1 <- abs(1 - P.plt) * norm
    nm.0 <- abs(P.plt) * norm
  } else if (criterion == "LCC"){
    nm <- abs(Y - P.plt)
    nm.1 <-  abs(1 - P.plt)
    nm.0 <-  abs(P.plt)
  }

  if (sampling.method == "WithReplacement"){
    dm <- sum(nm) # denominator
    # dm <- N * sum(nm[index.plt]) / n.plt
    p.ssp <- (1 - alpha) * nm / dm + alpha / N
    index.ssp <- sample(1:N, n.ssp, replace = T, prob = p.ssp)
    w.ssp <- 1 / p.ssp[index.ssp]
    pi.1 <- nm.1 / dm
    pi.0 <- nm.0 / dm
    offset <- log(pi.1[index.ssp] / pi.0[index.ssp])
  } else if (sampling.method == "Poisson"){
    NPhi <- N * sum(nm[index.plt]) / n.plt # the estimator of sum(nm)
    p.ssp <- (1 - alpha) * nm / NPhi + alpha / N
    P0 <- n.ssp * p.ssp
    index.ssp <- poisson.index(N, P0)
    w.ssp <- 1 / pmin(P0[index.ssp], 1) # Q: might be wrong

    print(nm.1[index.ssp])

    pi.1 <- pmin(n.ssp * nm.1 / NPhi, 1)
    pi.0 <- pmin(n.ssp * nm.0 / NPhi, 1)
    offset <- log(pi.1[index.ssp] / pi.0[index.ssp])
    print(offset)
  }

  return (list(index.ssp = index.ssp,
               offset = offset,
               w.ssp = w.ssp))
}
###############################################################################
subsample.estimate <-  function(x.ssp, y.ssp,
                                w.ssp, offset,
                                beta.plt, estimate.method) {
  if (estimate.method == "Weighted"){
    beta.ssp <- logistic.coef.estimate(x.ssp, y.ssp, weights = w.ssp)$beta # get similar results with start
    beta.ssp <- unlist(beta.ssp)
    P.ssp  <- pbeta(x.ssp, beta.ssp)
    MN.ssp <- MN(x.ssp, P.ssp, w.ssp)
    Psi.ssp <- Psi(x.ssp, y.ssp, P.ssp, w.ssp)
  } else if (estimate.method == 'LogOddsCorrection'){
    beta.ssp <- logistic.coef.estimate(X = x.ssp,
                                       Y = y.ssp,
                                       beta.start = beta.plt, # get similar results without start
                                       offset = offset)$beta
    beta.ssp <- unlist(beta.ssp)
    P.ssp  <- pbeta(x.ssp, beta.ssp)
    MN.ssp <- MN(x.ssp, P.ssp, 1)
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
logistic.optimal.subsampling <- function(X, Y, n.plt, n.ssp,
                                         criterion = "optL",
                                         sampling.method = 'Poisson',
                                         estimate.method = "LogOddsCorrection",
                                         # MN_custom,
                                         alpha = 0,
                                         b = 2){
  N <- nrow(X)
  d <- ncol(X)
  N1 <- sum(Y)
  N0 <- N - N1

  if (!all(Y %in% c(0, 1))) {
    message("Y contains values other than 0 and 1. Please check Y.")
  }

  # criterion = optA, optL, LCC
  # sampling.method = WithReplacement, Poisson
  # estimate.method = Weighted, LogOddsCorrection

  if (estimate.method %in% c("Weighted", "LogOddsCorrection")){

    # pilot step
    plt.estimate.results <- pilot.estimate(X = X, Y = Y, n.plt = n.plt)
    beta.plt <- plt.estimate.results$beta.plt
    MN.plt <- plt.estimate.results$MN.plt
    Psi.plt <- plt.estimate.results$Psi.plt
    P.plt <- plt.estimate.results$P.plt
    index.plt <- plt.estimate.results$index.plt

    # subsampling step
    ssp.results <- subsampling(X = X,
                               Y = Y,
                               n.ssp = n.ssp,
                               alpha = alpha,
                               b = b,
                               criterion = criterion,
                               sampling.method = sampling.method,
                               MN.plt = MN.plt,
                               P.plt = P.plt,
                               index.plt = index.plt)
    index.ssp <- ssp.results$index.ssp
    w.ssp <- ssp.results$w.ssp
    offset <- ssp.results$offset

    # subsample estimation step
    ssp.estimate.reults <- subsample.estimate(X[index.ssp, ],
                                              Y[index.ssp],
                                              w.ssp = w.ssp,
                                              offset = offset,
                                              beta.plt = beta.plt,
                                              estimate.method = estimate.method)
    beta.ssp <- ssp.estimate.reults$beta.ssp
    MN.ssp <- ssp.estimate.reults$MN.ssp
    Psi.ssp <- ssp.estimate.reults$Psi.ssp

    # combine step
    # print(MN.plt)
    # print(MN.ssp)
    MNsolve <- solve(MN.plt + MN.ssp)
    beta.cmb <- c(MNsolve %*% (MN.plt %*% beta.plt + MN.ssp %*% beta.ssp))
    var.beta <- MNsolve %*% (Psi.plt + Psi.ssp) %*% MNsolve

    return(list(beta.plt = beta.plt,
                beta.ssp = beta.ssp,
                beta.cmb = beta.cmb,
                var.beta = var.beta,
                index.plt = index.plt,
                index.ssp = index.ssp))
  } else if (estimate.method %in% c("Uni", "UniW")){
    # Uni:  poisson sampling + unweighted likelihood estimation + pilot correction
    # UniW: poisson sampling + weighted likelihood estimation
    plt.estimate.results <- pilot.estimate(X = X,
                                           Y = Y,
                                           n.plt = n.ssp)
    beta.ssp <- plt.estimate.results$beta.plt
    index.ssp <- plt.estimate.results$index.plt
    return(list(index.ssp = index.ssp,
                beta.ssp = beta.ssp))
  }
}
