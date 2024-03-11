# criterion = optA, optL, LCC
# sampling.method = WithReplacement, Poisson
# estimate.method = Weighted, LogOddsCorrection
getMSLE <- function(X, Y, offset = NULL, start) {
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
    ifelse(is.null(offset),
           p <- 1 / (1 + exp(-as.vector(X %*% beta))),
           p <- 1 / (1 + exp(-as.vector(X %*% beta) - offset)))
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
                                   offset = NULL,
                                   start = rep(0, ncol(X)),
                                   weights = 1) {
  if (length(weights) == 1)
    weights <- rep(weights, length(Y))
  beta <- rep(0, length(Y))
  data <- as.data.frame(cbind(Y, X))
  formula <- paste(colnames(data)[1], "~",
                   paste(colnames(data)[-1], collapse = "+"), "-1")
  design <- survey::svydesign(ids =  ~ 1,
                              weights =  ~ weights,
                              data = data)
  fit <- try(
    ifelse(is.null(offset),
           beta <- survey::svyglm(as.formula(formula),
                                  design = design,
                                  start = start,
                                  family = quasibinomial(link = "logit"))$coefficients,
           beta <- survey::svyglm(as.formula(formula),
                                  design = design,
                                  offset = offset,
                                  start = start,
                                  family = quasibinomial(link = "logit"))$coefficients),
    silent = TRUE)
  if ("try-error" %in% class(fit)) {
    message("Warning: an error occurred while calling 'svyglm': ", geterrmessage(),
            "This is probably due to the iteratively reweighted least squares algorithm called by 'glm.fit' not getting converged. Tring another function 'getMSLE' to replace 'glm.fit'.")
    ifelse(is.null(offset),
           beta <- getMSLE(X, Y, start = start),
           beta <- getMSLE(X, Y, offset = offset, start = start))
  }
  return(list(beta = beta))
}
###############################################################################
halfhalf.index <- function(N, Y, n.plt) {
  N1 <- sum(Y)
  N0 <- N - N1
  if (N0 < n.plt / 2 | N1 < n.plt / 2) { # if it happens, sampling from the rare events would be unnecessary
    warning(paste("n.plt/2 exceeds the number of Y=1 or Y=0 in the full data.",
                  "In this case, all rare events will be drawn into the pilot sample."))
  }
  n.plt.0 <- min(N0, n.plt / 2)
  n.plt.1 <- min(N1, n.plt / 2)
  if (n.plt.0 == n.plt.1) {
    index.plt <- c(sample(which(Y == 0), n.plt.0, replace = FALSE),
                   sample(which(Y == 1), n.plt.1, replace = FALSE))
  } else if (n.plt.0 < n.plt.1) {
    index.plt <- c(which(Y == 0), sample(which(Y == 1), n.plt.1, replace = FALSE))
  } else if (n.plt.0 > n.plt.1) {
    index.plt <- c(sample(which(Y == 0), n.plt.0, replace = FALSE), which(Y == 1))
  }
  return(list(index.plt = index.plt, n.plt.0 = n.plt.0, n.plt.1 = n.plt.1))
}
random.index <- function (N, n, p) {
  index <- sample(N, n, replace = TRUE, prob = p)
  return(index)
}
poisson.index <- function (N, pi) {
  index <- runif(N) <= pi
  return(which(index == TRUE))
}
MN <- function (X, pbeta, w) {
  phi <- pbeta * (1 - pbeta)
  MN <- t(X) %*% (X * (phi * w))
  return(MN)
}
Psi <- function (X, Y, pbeta, w) {
  psi <- (Y - pbeta)^2
  Psi <- t(X) %*% (X * (psi * w))
  return(Psi)
}
pbeta <- function (X, beta, offset = NULL) {
  ifelse(is.null(offset),
         p <- 1 / (1 + exp(-c(X %*% beta))),
         p <- 1 / (1 + exp(-c(X %*% beta) - offset)))
  return(p)
}
calculate.offset <- function (criterion, sampling.method, P.plt, dm, NPhi = NULL, n.ssp = NULL, norm = NULL) {
  if (criterion %in% c("optA", "optL")) {
    nm.1 <- abs(1 - P.plt) * norm
    nm.0 <- abs(P.plt) * norm
  } else if (criterion == "LCC") {
    nm.1 <-  abs(1 - P.plt)
    nm.0 <-  abs(P.plt)
  }
  if (sampling.method == 'WithReplacement') {
    pi.1 <- nm.1 / dm
    pi.0 <- nm.0 / dm
  } else if (sampling.method == 'Poisson') {
    pi.1 <- pmin(n.ssp * nm.1 / NPhi, 1)
    pi.0 <- pmin(n.ssp * nm.0 / NPhi, 1)
  }
  offset <- log(pi.1 / pi.0)
  return(offset)
}
###############################################################################
pilot.estimate <- function(X, Y, n.plt){
  N <- nrow(X)
  d <- ncol(X)
  N1 <- sum(Y)
  N0 <- N - N1

  # half half pilot estimator
  halfhalf.index.results <- halfhalf.index(N, Y, n.plt)
  index.plt <- halfhalf.index.results$index.plt
  n.plt.0 <- halfhalf.index.results$n.plt.0
  n.plt.1 <- halfhalf.index.results$n.plt.1
  x.plt <- X[index.plt,]
  y.plt <- Y[index.plt]
  p.plt <- c(rep(1 / (2 * N0), n.plt.0), rep(1 / (2 * N1), n.plt.1))

  ## case control pilot estimator
  ## if sampling uniformly, p.plt <- rep(1 / N, N), and correction is unnecessary
  # p.plt <- rep(1 / (2 * N1), N)
  # p.plt[which(Y == 0)] <- 1 / (2 * N0)
  # random.index.results <- random.index(N, n.plt, p.plt)
  # index.plt <- random.index.results$index
  # x.plt <- X[index.plt,]
  # y.plt <- Y[index.plt]
  # p.plt <- p.plt[index.plt]

  # unweighted log likelihood estimation
  beta.plt <- logistic.coef.estimate(X = x.plt, Y = y.plt)$beta # not corrected yet
  pbeta.plt <- pbeta(x.plt, beta.plt)
  ddL.plt <- MN(x.plt, pbeta.plt, 1 / n.plt) # 2nd derivative of the log likelihood function.
  Psi.plt <- Psi(x.plt, y.plt, pbeta.plt, 1 / n.plt^2)
  beta.plt[1] <- beta.plt[1] - log(N0 / N1)
  P.plt <- pbeta(X, beta.plt)
  MN.plt <- MN(x.plt, P.plt[index.plt], 1 / n.plt) # the pilot estimator of MN, used for calculating 'optA' subsampling probability

  return(
    list(
      p.plt = p.plt,
      beta.plt = beta.plt,
      ddL.plt = ddL.plt,
      Psi.plt = Psi.plt,
      MN.plt = MN.plt,
      P.plt = P.plt,
      index.plt = index.plt
    )
  )
}
###############################################################################
subsampling <- function(X, Y, n.ssp, alpha, b, criterion, sampling.method,
                        p.plt, MN.plt, P.plt, index.plt) {

  N <- nrow(X)
  N1 <- sum(Y)
  N0 <- N - N1
  n.plt <- length(index.plt) # length(index.plt) might be smaller than the n.plt user sets, so here we reset n.plt.

  if (criterion == "optA"){
    norm <- sqrt(colSums(solve(MN.plt, t(X))^2)) # the norm term of the criterion
    nm <- abs(Y - P.plt) * norm # numerator
  } else if (criterion == "optL"){
    norm <- sqrt(rowSums(X^2))
    nm <- abs(Y - P.plt) * norm
  } else if (criterion == "LCC"){
    nm <- abs(Y - P.plt)
  }

  if (sampling.method == "WithReplacement"){
    dm <- sum(nm) # denominator
    p.ssp <- (1 - alpha) * nm / dm + alpha / N
    index.ssp <- random.index(N, n.ssp, p.ssp)
    w.ssp <- 1 / p.ssp[index.ssp]
    offset <- calculate.offset(criterion = criterion,
                               sampling.method = sampling.method,
                               P.plt = P.plt[index.ssp],
                               dm = dm,
                               norm = norm[index.ssp])
  } else if (sampling.method == "Poisson"){
    # H <- quantile(nm, 1 - n.ssp / (b * N)) # if consider threshold
    # nm[nm > H] <- H
    NPhi <- sum(nm[index.plt] / p.plt) / n.plt
    p.ssp <- n.ssp * ((1 - alpha) * nm / NPhi + alpha / N)
    index.ssp <- poisson.index(N, p.ssp)
    w.ssp <- 1 / pmin(p.ssp[index.ssp], 1)
    offset <- calculate.offset(criterion = criterion,
                               sampling.method = sampling.method,
                               P.plt = P.plt[index.ssp],
                               dm = dm,
                               NPhi = NPhi,
                               n.ssp = n.ssp,
                               norm = norm[index.ssp])
  }

  return(list(index.ssp = index.ssp,
              offset = offset,
              w.ssp = w.ssp))
}
###############################################################################
subsample.estimate <- function(x.ssp, y.ssp, n.ssp,
                               w.ssp, offset,
                               beta.plt, sampling.method, estimate.method) {
  if (estimate.method == "Weighted") {
    beta.ssp <- logistic.coef.estimate(x.ssp, y.ssp, weights = w.ssp)$beta
    P.ssp <- pbeta(x.ssp, beta.ssp)
    if (sampling.method == 'Poisson') {
      ddL.ssp <- MN(x.ssp, P.ssp, w.ssp / N) # 2nd derivative of the log likelihood function.
      Psi.ssp <- Psi(x.ssp, y.ssp, P.ssp, w.ssp ^ 2 / N ^ 2)
      Lambda.ssp <- 0 # holding
    } else if (sampling.method == "WithReplacement") {
      ddL.ssp <- MN(x.ssp, P.ssp, w.ssp / (N * n.ssp))
      Psi.ssp <- Psi(x.ssp, y.ssp, P.ssp, w.ssp ^ 2 / (N ^ 2 * n.ssp ^ 2))
      c <- n.ssp / N
      Lambda.ssp <- c * Psi(x.ssp, y.ssp, P.ssp, w.ssp / (N * n.ssp ^ 2))
      # print(Lambda.ssp / Psi.ssp)
    }
  } else if (estimate.method == 'LogOddsCorrection') {
    beta.ssp <- logistic.coef.estimate(X = x.ssp,
                                       Y = y.ssp,
                                       start = beta.plt,
                                       offset = offset)$beta
    P.ssp  <- pbeta(x.ssp, beta.ssp, offset)
    ddL.ssp <- MN(x.ssp, P.ssp, 1 / n.ssp)
    Psi.ssp <- Psi(x.ssp, y.ssp, P.ssp, 1 / n.ssp^2)
    Lambda.ssp <- 0 # holding
  }
  return(list(beta.ssp = beta.ssp,
              ddL.ssp = ddL.ssp,
              Psi.ssp = Psi.ssp,
              Lambda.ssp = Lambda.ssp)
  )
}
###############################################################################
combining <- function(ddL.plt, ddL.ssp, Psi.plt, Psi.ssp, Lambda.ssp, n.plt, n.ssp, beta.plt, beta.ssp) {
  ddL.plt <- n.plt * ddL.plt
  ddL.ssp <- n.ssp * ddL.ssp
  Psi.plt <- n.plt ^ 2 * Psi.plt
  Psi.ssp <- n.ssp ^ 2 * Psi.ssp
  Lambda.ssp <- n.ssp ^ 2 * Lambda.ssp
  MNsolve <- solve(ddL.plt + ddL.ssp)
  beta.cmb <- c(MNsolve %*% (ddL.plt %*% beta.plt + ddL.ssp %*% beta.ssp))
  var.cmb <- MNsolve %*% (Psi.plt + Psi.ssp) %*% MNsolve
  var.cmb.true <- MNsolve %*% (Psi.plt + Psi.ssp + Lambda.ssp) %*% MNsolve
  return(list(beta.cmb = beta.cmb,
              var.cmb = var.cmb,
              var.cmb.true = var.cmb.true)
  )
}
###############################################################################
logistic.optimal.subsampling <- function(formula, data, n.plt, n.ssp,
                                         criterion = "optL",
                                         sampling.method = 'Poisson',
                                         estimate.method = "LogOddsCorrection",
                                         # MN_custom, # holding
                                         alpha = 0,
                                         b = 2) {

  model.call <- match.call()
  mf <- model.frame(formula, data)
  # Extract the response and predictor variables from the model frame
  Y <- model.response(mf, "any") # if we only consider logistic model ,it should not be 'any'.
  X <- model.matrix(formula, mf)

  N <- nrow(X)
  d <- ncol(X)
  N1 <- sum(Y)
  N0 <- N - N1

  if (estimate.method %in% c("Weighted", "LogOddsCorrection")) {

    # pilot step
    plt.estimate.results <- pilot.estimate(X = X, Y = Y, n.plt = n.plt)
    p.plt <- plt.estimate.results$p.plt
    beta.plt <- plt.estimate.results$beta.plt
    ddL.plt <- plt.estimate.results$ddL.plt
    Psi.plt <- plt.estimate.results$Psi.plt
    MN.plt <- plt.estimate.results$MN.plt
    P.plt <- plt.estimate.results$P.plt
    index.plt <- plt.estimate.results$index.plt
    var.plt <- solve(ddL.plt) %*% Psi.plt %*% solve(ddL.plt)

    # subsampling step
    ssp.results <- subsampling(X = X,
                               Y = Y,
                               n.ssp = n.ssp,
                               alpha = alpha,
                               b = b,
                               criterion = criterion,
                               sampling.method = sampling.method,
                               p.plt = p.plt,
                               MN.plt = MN.plt,
                               P.plt = P.plt,
                               index.plt = index.plt)
    index.ssp <- ssp.results$index.ssp
    w.ssp <- ssp.results$w.ssp
    offset <- ssp.results$offset

    # subsample estimating step
    ssp.estimate.reults <- subsample.estimate(X[index.ssp, ],
                                              Y[index.ssp],
                                              n.ssp = n.ssp,
                                              w.ssp = w.ssp,
                                              offset = offset,
                                              beta.plt = beta.plt,
                                              sampling.method = sampling.method,
                                              estimate.method = estimate.method)
    beta.ssp <- ssp.estimate.reults$beta.ssp
    ddL.ssp <- ssp.estimate.reults$ddL.ssp
    Psi.ssp <- ssp.estimate.reults$Psi.ssp
    Lambda.ssp <- ssp.estimate.reults$Lambda.ssp
    var.ssp <- solve(ddL.ssp) %*% Psi.ssp %*% solve(ddL.ssp)
    var.ssp.true <- solve(ddL.ssp) %*% (Psi.ssp + Lambda.ssp) %*% solve(ddL.ssp)

    # combining step
    combining.results <- combining(ddL.plt = ddL.plt,
                                   ddL.ssp = ddL.ssp,
                                   Psi.plt = Psi.plt,
                                   Psi.ssp = Psi.ssp,
                                   Lambda.ssp = Lambda.ssp,
                                   n.plt = n.plt,
                                   n.ssp = n.ssp,
                                   beta.plt = beta.plt,
                                   beta.ssp = beta.ssp)
    beta.cmb <- combining.results$beta.cmb
    var.cmb <- combining.results$var.cmb
    var.cmb.true <- combining.results$var.cmb.true

    return(list(beta.plt = beta.plt,
                beta.ssp = beta.ssp,
                beta.cmb = beta.cmb,
                var.plt = var.plt,
                var.ssp = var.ssp,
                var.ssp.true = var.ssp.true,
                var.cmb = var.cmb,
                var.cmb.true = var.cmb.true,
                index.plt = index.plt,
                index.ssp = index.ssp,
                ddL.plt = ddL.plt,
                ddL.ssp = ddL.ssp))
  } else if (estimate.method == "Uni"){
    index.ssp <- random.index(N, n.plt + n.ssp, rep(1 / N, N))
    beta.ssp <- logistic.coef.estimate(X = X[index.ssp,], Y = Y[index.ssp])$beta
    return(list(index.ssp = index.ssp,
                beta.ssp = beta.ssp))
  }
}
