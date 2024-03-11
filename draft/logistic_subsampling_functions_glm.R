# criterion = optA, optL, LCC
# sampling.method = WithReplacement, Poisson
# estimate.method = Weighted, LogOddsCorrection
getMSLE <- function(X, Y, offset = NULL, start, weights = 1) {
  n <- nrow(X)
  d <- ncol(X)
  beta <- start
  wX <- weights * X
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
    S <- (Y - p) * wX
    phi <- p * (1 - p)
    H <- t(X) %*% (phi * wX)
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

  H.solve <- solve(H)
  cov <- H.solve %*% (t(S) %*% S) %*% H.solve

  return(list(coefficients = beta,  # Use the same variable names as 'glm' ouput
              cov.unscaled = cov,
              fitted.values = p
              ))
}
###############################################################################
logistic.coef.estimate <- function(X,
                                   Y,
                                   offset = NULL,
                                   start = rep(0, ncol(X)),
                                   weights = 1) {

  data <- as.data.frame(cbind(Y, X))
  formula <- paste(colnames(data)[1], "~",
                   paste(colnames(data)[-1], collapse = "+"), "-1") # use '-1' to avoid adding intercept column again.
  design <- survey::svydesign(ids =  ~ 1,
                              weights =  ~ weights,
                              data = data)
  beta <- rep(0, length(Y))
  fit <- try(
    ifelse(is.null(offset),
           results <- survey::svyglm(as.formula(formula),
                                     design = design,
                                     start = start,
                                     family = quasibinomial(link = "logit")),
           results <- survey::svyglm(as.formula(formula),
                                     design = design,
                                     start = start,
                                     offset = offset,
                                     family = quasibinomial(link = "logit"))),
           silent = TRUE
   )
  if ("try-error" %in% class(fit)) {
    message("Warning: an error occurred while calling 'svyglm': ",
            geterrmessage(),
            "Tring another function 'getMSLE' to replace 'glm.fit'.")
    ifelse(is.null(offset),
           results <- getMSLE(X, Y, start = start, weights = weights),
           results <- getMSLE(X, Y, offset = offset, start = start, weights = weights))
  }

  beta <- results$coefficients
  cov <- results$cov.unscaled
  pbeta <- as.vector(results$fitted.values)

  return(list(
    beta = beta,
    cov = cov,
    pbeta = pbeta
  ))
}
###############################################################################
halfhalf.index <- function(N, Y, n.plt) {
  N1 <- sum(Y)
  N0 <- N - N1
  if (N0 < n.plt / 2 | N1 < n.plt / 2) {
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
MN <- function (X, pbeta, w = 1) {
  phi <- pbeta * (1 - pbeta)
  MN <- t(X) %*% (X * (phi * w))
  return(MN)
}
Psi <- function (X, Y, pbeta, w = 1) {
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
calculate.offset <- function (criterion, sampling.method, P.plt, dm, NPhi = NULL, n.ssp = NULL, norm = NULL,
                              alpha = alpha, N = N) {
  if (criterion %in% c("optA", "optL")) {
    nm.1 <- abs(1 - P.plt) * norm
    nm.0 <- abs(P.plt) * norm
  } else if (criterion == "LCC") {
    nm.1 <-  abs(1 - P.plt)
    nm.0 <-  abs(P.plt)
  }
  if (sampling.method == 'WithReplacement') {
    warning('estimate.method="LogOddsCorrection" is not well defined when sampling.method="WithReplacement"')
    pi.1 <- (1 - alpha) * nm.1 / dm +  alpha / N
    pi.0 <- (1 - alpha) * nm.0 / dm +  alpha / N
  } else if (sampling.method == 'Poisson') {
    pi.1 <- pmin(n.ssp * ((1 - alpha) * nm.1 / NPhi + alpha / N), 1)
    pi.0 <- pmin(n.ssp * ((1 - alpha) * nm.0 / NPhi + alpha / N), 1)
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

  ## half half pilot estimator
  halfhalf.index.results <- halfhalf.index(N, Y, n.plt)
  index.plt <- halfhalf.index.results$index.plt
  n.plt.0 <- halfhalf.index.results$n.plt.0
  n.plt.1 <- halfhalf.index.results$n.plt.1
  x.plt <- X[index.plt,]
  y.plt <- Y[index.plt]
  p.plt <- c(rep(1 / (2 * N0), n.plt.0), rep(1 / (2 * N1), n.plt.1))

  # unweighted log likelihood estimation
  results.plt <- logistic.coef.estimate(X = x.plt, Y = y.plt)
  beta.plt <- results.plt$beta # not corrected yet
  pbeta.plt <- results.plt$pbeta
  var.plt <- results.plt$cov

  # pbeta.plt <- pbeta(x.plt, beta.plt)
  ddL.plt <- MN(x.plt, pbeta.plt, 1 / n.plt) # 2nd derivative of the log likelihood function.
  Psi.plt <- Psi(x.plt, y.plt, pbeta.plt, 1 / n.plt^2)
  # var.plt <- solve(ddL.plt) %*% Psi.plt %*% solve(ddL.plt)
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
      index.plt = index.plt,
      var.plt = var.plt
    )
  )
}
###############################################################################
subsampling <- function(X, Y, n.ssp, alpha, b, criterion, estimate.method, sampling.method,
                        p.plt, MN.plt, P.plt, index.plt) {

  N <- nrow(X)
  N1 <- sum(Y)
  N0 <- N - N1
  n.plt <- length(index.plt) # length(index.plt) might be smaller than the n.plt user sets, so here we reset n.plt.
  w.ssp <- offset <- NA
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
    # boxplot(nm / dm)
    p.ssp <- (1 - alpha) * nm / dm + alpha / N
    # boxplot(p.ssp)
    index.ssp <- random.index(N, n.ssp, p.ssp)

    if (estimate.method == 'LogOddsCorrection') {
      offset <- calculate.offset(criterion = criterion,
                                 sampling.method = sampling.method,
                                 P.plt = P.plt[index.ssp],
                                 dm = dm,
                                 norm = norm[index.ssp],
                                 alpha = alpha,
                                 N = N)
    } else if (estimate.method == 'Weighted') {
      w.ssp <- 1 / p.ssp[index.ssp]
    }

  } else if (sampling.method == "Poisson"){
    # H <- quantile(nm, 1 - n.ssp / (b * N)) # if consider threshold
    # nm[nm > H] <- H
    NPhi <- sum(nm[index.plt] / p.plt) / n.plt
    # boxplot(n.ssp * nm / NPhi)
    p.ssp <- n.ssp * ((1 - alpha) * nm / NPhi + alpha / N)
    # boxplot(p.ssp)
    index.ssp <- poisson.index(N, p.ssp)

    if (estimate.method == 'LogOddsCorrection') {
      offset <- calculate.offset(criterion = criterion,
                                 sampling.method = sampling.method,
                                 P.plt = P.plt[index.ssp],
                                 dm = dm,
                                 NPhi = NPhi,
                                 n.ssp = n.ssp,
                                 norm = norm[index.ssp],
                                 alpha = alpha,
                                 N = N)
    } else if (estimate.method == 'Weighted') {
      w.ssp <- 1 / pmin(p.ssp[index.ssp], 1)
    }
  }

  return(list(index.ssp = index.ssp,
              offset = offset,
              w.ssp = w.ssp))
}
###############################################################################
subsample.estimate <- function(x.ssp, y.ssp, n.ssp,
                               w.ssp, offset,
                               beta.plt, sampling.method, estimate.method) {

  d <- ncol(x.ssp)

  if (estimate.method == "Weighted") {

    results.ssp <- logistic.coef.estimate(x.ssp, y.ssp, weights = w.ssp)
    beta.ssp <- results.ssp$beta
    P.ssp <- results.ssp$pbeta
    var.ssp <- results.ssp$cov

    # beta.ssp <- logistic.coef.estimate(x.ssp, y.ssp, weights = w.ssp)$beta
    # P.ssp <- pbeta(x.ssp, beta.ssp)

    if (sampling.method == 'Poisson') {
      ddL.ssp <- MN(x.ssp, P.ssp, w.ssp / N) # 2nd derivative of the log likelihood function.
      Psi.ssp <- Psi(x.ssp, y.ssp, P.ssp, w.ssp^2 / N^2) # the estimator of the variance of 1st derivative
      Lambda.ssp <- matrix(0, d, d) # holding
    } else if (sampling.method == "WithReplacement") {
      ddL.ssp <- MN(x.ssp, P.ssp, w.ssp / (N * n.ssp))
      Psi.ssp <- Psi(x.ssp, y.ssp, P.ssp, w.ssp ^ 2 / (N ^ 2 * n.ssp ^ 2))
      c <- n.ssp / N
      Lambda.ssp <- c * Psi(x.ssp, y.ssp, P.ssp, w.ssp / (N * n.ssp ^ 2))
      # print(Lambda.ssp / Psi.ssp)

    }
  } else if (estimate.method == 'LogOddsCorrection') {
    results.ssp <- logistic.coef.estimate(X = x.ssp,
                                          Y = y.ssp,
                                          start = beta.plt,
                                          offset = offset)
    beta.ssp <- results.ssp$beta
    P.ssp <- results.ssp$pbeta
    var.ssp <- results.ssp$cov
    # beta.ssp <- logistic.coef.estimate(X = x.ssp,
    #                                    Y = y.ssp,
    #                                    start = beta.plt,
    #                                    offset = offset)$beta
    # P.ssp  <- pbeta(x.ssp, beta.ssp, offset)
    ddL.ssp <- MN(x.ssp, P.ssp, 1 / n.ssp)
    Psi.ssp <- Psi(x.ssp, y.ssp, P.ssp, 1 / n.ssp ^ 2)
    Lambda.ssp <- matrix(0, d, d) # holding
  }

  # var.ssp <- solve(ddL.ssp) %*% Psi.ssp %*% solve(ddL.ssp)
  # var.ssp.true <- solve(ddL.ssp) %*% (Psi.ssp + Lambda.ssp) %*% solve(ddL.ssp)
  var.ssp.true <- var.ssp + solve(ddL.ssp) %*% (Lambda.ssp) %*% solve(ddL.ssp)
  # print('var.ssp + solve(ddL.ssp) %*% (Lambda.ssp) %*% solve(ddL.ssp)')
  # print(var.ssp + solve(ddL.ssp) %*% (Lambda.ssp) %*% solve(ddL.ssp))
  return(list(beta.ssp = beta.ssp,
              ddL.ssp = ddL.ssp,
              Psi.ssp = Psi.ssp,
              var.ssp = var.ssp,
              var.ssp.true = var.ssp.true,
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
  # print(MNsolve %*% ddL.plt %*% beta.plt)
  # print(MNsolve %*% ddL.ssp %*% beta.ssp)
  var.cmb <- MNsolve %*% (Psi.plt + Psi.ssp) %*% MNsolve
  var.cmb.true <- MNsolve %*% (Psi.plt + Psi.ssp + Lambda.ssp) %*% MNsolve
  return(list(beta.cmb = beta.cmb,
              var.cmb = var.cmb,
              var.cmb.true = var.cmb.true)
  )
}
###############################################################################
format.p.values <- function(p.values, threshold = 0.0001) {
  formatted <- sapply(p.values, function(p.value) {
    if (p.value < threshold) {
      return(sprintf("<%.4f", threshold))
    } else {
      return(sprintf("%.4f", p.value))
    }
  })
  return(formatted)
}
###############################################################################
subsampling.summary <- function(object) {
  print(object$beta)
  coef <- object$beta
  se <- sqrt(diag(object$var))
  pred <- object$prediction.results
  cat("Model Summary\n\n")
  cat("\nCall:\n")
  cat("\n")
  print(object$model.call)
  cat("\n")
  cat("Coefficients:\n")
  cat("\n")
  coef.table <- data.frame(
    Estimate = round(coef, digits = 4),
    `Std. Error` = round(se, digits = 4),
    `z value` = round(coef / se, digits = 4),
    `Pr(>|z|)` = format.p.values(2 * (1 - pnorm(abs(coef / se))), threshold = 0.0001),
    check.names = FALSE
  )
  rownames(coef.table) <- names(coef)
  print(coef.table)
  cat("\n")
  cat("Predictions:\n")
  pred.table <- data.frame(Value = format(unlist(pred), digits = 4))
  print(pred.table)

  # Add more summary information as needed
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
  colnames(X)[1] <- "intercept"
  N <- nrow(X)
  d <- ncol(X)
  N1 <- sum(Y)
  N0 <- N - N1
  # print(d)
  # print(X[1:5,])
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
    var.plt <- plt.estimate.results$var.plt

    # subsampling step
    ssp.results <- subsampling(X = X,
                               Y = Y,
                               n.ssp = n.ssp,
                               alpha = alpha,
                               b = b,
                               criterion = criterion,
                               estimate.method = estimate.method,
                               sampling.method = sampling.method,
                               p.plt = p.plt,
                               MN.plt = MN.plt,
                               P.plt = P.plt,
                               index.plt = index.plt)
    index.ssp <- ssp.results$index.ssp
    w.ssp <- ssp.results$w.ssp
    offset <- ssp.results$offset

    # subsample estimating step
    ssp.estimate.results <- subsample.estimate(X[index.ssp, ],
                                               Y[index.ssp],
                                               n.ssp = n.ssp,
                                               w.ssp = w.ssp,
                                               offset = offset,
                                               beta.plt = beta.plt,
                                               sampling.method = sampling.method,
                                               estimate.method = estimate.method)
    beta.ssp <- ssp.estimate.results$beta.ssp
    ddL.ssp <- ssp.estimate.results$ddL.ssp
    Psi.ssp <- ssp.estimate.results$Psi.ssp
    Lambda.ssp <- ssp.estimate.results$Lambda.ssp
    var.ssp <- ssp.estimate.results$var.ssp
    var.ssp.true <- ssp.estimate.results$var.ssp.true

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

    # setClass("ModelSummary",
    #          slots = c(
    #            coefficients = "numeric",
    #            se = "numeric",
    #            call = "call"
    #          )
    # )
    # result <- model.summary(beta.cmb, sqrt(diag(var.cmb.true)), model.call)
    # return(result)

    # prediction.results <- in.sample.prediction(X, Y, beta.cmb)

    return(list(model.call = model.call,
                beta.plt = beta.plt,
                beta.ssp = beta.ssp,
                beta.cmb = beta.cmb,
                var.plt = var.plt,
                var.ssp = var.ssp,
                var.ssp.true = var.ssp.true,
                var.cmb = var.cmb,
                var.cmb.true = var.cmb.true,
                index.plt = index.plt,
                index.ssp = index.ssp))
  } else if (estimate.method == "Uni"){
    n.uni <- n.plt + n.ssp
    index.uni <- random.index(N, n.uni, rep(1 / N, N))
    X.uni <- X[index.uni, ]
    Y.uni = Y[index.uni]

    results.uni <- logistic.coef.estimate(X = X.uni, Y = Y.uni) # not corrected yet
    beta.uni <- results.uni$beta
    pbeta.uni <- results.uni$pbeta
    var.uni <- results.uni$cov

    # pbeta.uni <- pbeta(X[index.uni, ], beta.uni)
    # ddL.uni <- MN(X.uni, pbeta.uni, 1 / n.uni) # 2nd derivative of the log likelihood function.
    # Psi.uni <- Psi(X.uni, Y.uni, pbeta.uni, 1 / n.uni^2)
    # var.uni <- solve(ddL.uni) %*% Psi.uni %*% solve(ddL.uni)

    # print(var.uni)
    # print(var.uni.glm)
    return(list(model.call = model.call,
                index = index.uni,
                beta = beta.uni,
                var = var.uni))
  }
}
