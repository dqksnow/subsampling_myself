getMSLE <- function(X, Y, offset = NULL, start, weights = 1) {
  n <- nrow(X)
  d <- ncol(X)
  beta <- start
  S <- matrix(0, n, d)
  H <- matrix(0, d, d)
  wX <- weights * X
  loop <- 1
  Loop <- 100
  step.length <- 0.5
  msg <- "NA"
  convergence_flag <- FALSE # determine whether the loop should proceed or terminate.
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
        convergence_flag <- TRUE
        return(rep(NA, d))
      }
    )
    if (!convergence_flag) {
      beta.new <- beta + shs * step.length
      tlr <- sum(shs ^ 2)
      beta <- beta.new
      if (tlr < .00001) {
        msg <- "Successful convergence"
        convergence_flag <- TRUE
      }
      if (loop == Loop) {
        msg <- "Maximum iteration reached"; cat(msg, "\n")
        convergence_flag <- TRUE
      }
    }
    if (convergence_flag) {
      break
    }
    loop <- loop + 1
  }

  H.solve <- solve(H)
  cov <- H.solve %*% (t(S) %*% S) %*% H.solve

  return(list(coefficients = beta,  # Use the same variable names as 'glm' ouput
              cov.unscaled = cov,
              fitted.values = p))
}
###############################################################################
glm.coef.estimate <- function(X,
                              Y,
                              offset = NULL,
                              start = rep(0, ncol(X)),
                              weights = 1) {
  data <- as.data.frame(cbind(Y, X))
  formula <- as.formula(paste(colnames(data)[1], "~",
                   paste(colnames(data)[-1], collapse = "+"), "-1"))
  # use '-1' to avoid adding intercept column again.
  design <- survey::svydesign(ids =  ~ 1,
                              weights =  ~ weights,
                              data = data)
  fit <- ifelse(is.null(offset),
                results <- survey::svyglm(formula,
                                          design = design,
                                          start = start,
                                          family = quasibinomial(link="logit")),
                results <- survey::svyglm(formula,
                                          design = design,
                                          start = start,
                                          offset = offset,
                                          family = quasibinomial(link="logit"))
  )
  beta <- results$coefficients
  cov <- results$cov.unscaled
  pbeta <- as.vector(results$fitted.values)
  return(list(beta = beta,
              cov = cov,
              pbeta = pbeta))
}
###############################################################################
halfhalf.index <- function(N, Y, n.plt) {
  N1 <- sum(Y)
  N0 <- N - N1
  if (N0 < n.plt / 2 | N1 < n.plt / 2) {
    warning(paste("n.plt/2 exceeds the number of Y=1 or Y=0 in the full data.",
                  "All rare events will be drawn into the pilot sample."))
  }
  n.plt.0 <- min(N0, n.plt / 2)
  n.plt.1 <- min(N1, n.plt / 2)
  if (n.plt.0 == n.plt.1) {
    index.plt <- c(sample(which(Y == 0), n.plt.0, replace = FALSE),
                   sample(which(Y == 1), n.plt.1, replace = FALSE))
  } else if (n.plt.0 < n.plt.1) {
    index.plt <- c(which(Y == 0),
                   sample(which(Y == 1),
                          n.plt.1,
                          replace = FALSE))
  } else if (n.plt.0 > n.plt.1) {
    index.plt <- c(sample(which(Y == 0),
                          n.plt.0,
                          replace = FALSE),
                   which(Y == 1))
  }
  return(list(index.plt = index.plt, n.plt.0 = n.plt.0, n.plt.1 = n.plt.1))
}
poisson.index <- function(N, pi){
  return(which(runif(N) <= pi))
}
ddL <- function(X, P, w=1){
  phi <- P * (1 - P)
  t(X) %*% (X * (phi * w))
}
dL.sq <- function(X, Y, P, w=1){
  dL.sq <- (Y - P)^2
  t(X) %*% (X * (dL.sq * w))
}
pbeta <- function(X, beta, offset = NA){
  ifelse(is.na(offset),
         p <- 1 / (1 + exp(-as.vector(X %*% beta))),
         p <- 1 / (1 + exp(-as.vector(X %*% beta) - offset)))
  return (p)
}
calculate.nm <- function(X, Y, ddL.plt.correction, P.plt, criterion){
  if (criterion == "OptA"){
    temp <- X %*% t(solve(ddL.plt.correction))
    nm <- sqrt(rowSums(temp^2)) # norm
    nm <- P.plt * sqrt(1 - P.plt) * nm # numerator
  } else if (criterion == "OptL"){
    nm <- sqrt(rowSums(X^2))
    nm <- P.plt * sqrt(1 - P.plt) * nm
  } else if (criterion == "LCC"){
    nm <- abs(Y - P.plt)
  }
  return(nm)
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
  coef <- object$beta
  se <- sqrt(diag(object$var))
  # pred <- object$prediction.results
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
  # cat("\n")
  # cat("Predictions:\n")
  # pred.table <- data.frame(Value = format(unlist(pred), digits = 4))
  # print(pred.table)

  # Add more summary information as needed
}
###############################################################################
# in.sample.prediction <- function(X, Y, beta){
#   N <- length(Y)
#   pbeta.estimate <- pbeta(X, beta)
#   Y.pred <- rep(0, N)
#   Y.pred[which(pbeta.estimate > 0.5)] <- 1
#   plot(pbeta.estimate)
#   # Calculate True Positives, True Negatives, False Positives, False Negatives
#   tp <- sum(Y == 1 & Y.pred == 1)
#   tn <- sum(Y == 0 & Y.pred == 0)
#   fp <- sum(Y == 0 & Y.pred == 1)
#   fn <- sum(Y == 1 & Y.pred == 0)
#   print(sum(Y.pred == 1))
#   # accuracy measures the model's overall performance of prediction.
#   accuracy <- (tp + tn) / N
#   # precision measures the model's performance when it predicts positive outcomes(Y=1).
#   precision <- tp / (tp + fp)
#   # recall measures that of all the actual positive (Y=1) instances, how many did the model correctly predict as positive.
#   recall <- tp / (tp + fn) # tp + fn = N1
#   # f1.score is the harmonic mean of precision and recall.
#   f1.score <- 2 * (precision * recall) / (precision + recall)
#   roc_obj <- roc(Y, pbeta.estimate)
#   print(auc(roc_obj))
#   plot(roc_obj, main = "ROC Curve", print.auc = TRUE)
#   return(list(
#     Accuracy = accuracy,
#     Precision = precision,
#     Recall = recall,
#     F1.Score = f1.score
#   ))
# }
###############################################################################
pilot.estimate <- function(X, Y, n.plt){
  N <- nrow(X)
  d <- ncol(X)
  N1 <- sum(Y)
  N0 <- N - N1
  # half half sampling
  halfhalf.index.results <- halfhalf.index(N, Y, n.plt)
  index.plt <- halfhalf.index.results$index.plt
  n.plt.0 <- halfhalf.index.results$n.plt.0
  n.plt.1 <- halfhalf.index.results$n.plt.1
  x.plt <- X[index.plt, ]
  y.plt <- Y[index.plt]
  p.plt <- c(rep(1 / (2 * N0), n.plt.0), rep(1 / (2 * N1), n.plt.1))
  # lower p: sampling probability
  # unweighted likelihood estimation
  results.plt <- glm.coef.estimate(X = x.plt, Y = y.plt)
  beta.plt <- results.plt$beta
  pbeta.plt <- results.plt$pbeta

  print(beta.plt)
  # print(pbeta.plt)

  ddL.plt <- ddL(x.plt, pbeta.plt, 1 / n.plt)
  dL.sq.plt <- dL.sq(x.plt, y.plt, pbeta.plt, 1 / n.plt ^ 2)
  # correction
  beta.plt[1] <- beta.plt[1] - log(N0 / N1)
  P.plt <- pbeta(X, beta.plt)
  # capital P: Prob(Y=1)
  ddL.plt.correction <- ddL(x.plt, P.plt[index.plt], 1 / n.plt) # the estimator of MN, used for calculating 'OptA' subsampling probability

  ## weighted likelihood
  # weights = 1 / p.plt
  # beta.plt <- glm.coef.estimate(X = x.plt, Y = y.plt, weights = weights)$beta
  # pbeta.plt <- pbeta(x.plt, beta.plt)
  # ddL.plt <- MN(x.plt, pbeta.plt, weights / (N * n.plt)) # the 2nd derivative of the log likelihood function.
  # dL.sq.plt <- dL.sq(x.plt, y.plt, pbeta.plt, weights^2 / (N ^ 2 * n.plt ^ 2))
  # P.plt <- pbeta(X, beta.plt)
  # ddL.plt.correction <- MN(x.plt, P.plt[index.plt], weights / (N * n.plt)) # the pilot estimator of MN, used for calculating 'OptA' subsampling probability


  return(list(p.plt = p.plt,
              beta.plt = beta.plt,
              ddL.plt = ddL.plt,
              dL.sq.plt = dL.sq.plt,
              ddL.plt.correction = ddL.plt.correction,
              P.plt = P.plt,
              index.plt = index.plt
  ))
}
###############################################################################
subsampling <- function(X, Y, n.ssp, alpha, b, criterion, estimate.method,
                        p.plt, ddL.plt.correction, P.plt, index.plt) {
  N <- nrow(X)
  N1 <- sum(Y)
  N0 <- N - N1
  n.plt <- length(index.plt)
  # length(index.plt) might be smaller than n.plt, so here we reset n.plt.
  w.ssp <- offset <- NA
  nm <- calculate.nm(X, Y, ddL.plt.correction, P.plt, criterion)

  # Currently only Poisson sampling method has been implemented.
  if(criterion %in% c('OptL', 'OptA')){
    H <- quantile(nm[index.plt], 1 - n.ssp / (b * N)) # if consider threshold
    nm[nm > H] <- H
    NPhi <- (N0 / N) * sum(nm[index.plt] / p.plt) / n.plt
    p.ssp <- n.ssp * ((1 - alpha) * nm / NPhi + alpha / N)
    index.ssp <- poisson.index(N, Y + (1 - Y) * p.ssp)
    p.ssp <- pmin(p.ssp[index.ssp], 1)
    # calculate offset or weights:
    if (estimate.method == 'LogOddsCorrection') {
      offset <- -log(p.ssp)
    } else if (estimate.method == 'Weighted') {
      w.ssp <- 1 / (Y[index.ssp] + (1 - Y[index.ssp]) * p.ssp)
    }
  } else if (criterion == 'LCC'){
    dm <- sum(nm[index.plt] / p.plt) / n.plt
    p.ssp <- (n.ssp + N1) * nm / dm
    index.ssp <- poisson.index(N, p.ssp)
    p.ssp <- pmin(p.ssp[index.ssp], 1)
    # calculate offset or weights:
    if (estimate.method == 'LogOddsCorrection') {
      nm.1 <- abs(1 - P.plt[index.ssp])
      nm.0 <- abs(P.plt[index.ssp])
      pi.1 <- pmin((n.ssp + N1) * nm.1 / dm, 1)
      pi.0 <- pmin((n.ssp + N1) * nm.0 / dm, 1)
      offset <- log(pi.1 / pi.0)
    } else if (estimate.method == 'Weighted') {
      w.ssp <- 1 / p.ssp
    }
  }

  # if (criterion == "OptA"){
  #   norm <- sqrt(colSums(solve(ddL.plt.correction, t(X))^2)) # the norm term of the criterion
  #   nm <- P.plt * sqrt(1 - P.plt) * norm # numerator
  # } else if (criterion == "OptL"){
  #   norm <- sqrt(rowSums(X^2))
  #   nm <- P.plt * sqrt(1 - P.plt) * norm
  # } else if (criterion == "LCC"){
  #   nm <- abs(Y - P.plt)
  # }

  # # only Poisson sampling method is considered
  # # H <- quantile(nm, 1 - n.ssp / (b * N)) # if consider threshold
  # # nm[nm > H] <- H
  # NPhi <- sum(nm[index.plt] / p.plt) / n.plt
  # # NPhi <- sum(nm)
  # # NPhi <- (N0 / N) * sum(nm[index.plt] / p.plt) / n.plt
  # p.ssp <- n.ssp * ((1 - alpha) * nm / NPhi + alpha / N)
  # print(paste('sum(nm / NPhi): ', sum(nm / NPhi))) # 1
  # print(paste('sum(Y + (1 - Y) * pmin(p.ssp, 1))): ', sum(Y + (1 - Y) * pmin(p.ssp, 1)))) # N1 + p.ssp
  # index.ssp <- poisson.index(N, Y + (1 - Y) * p.ssp)
  # print(paste('length(index.ssp)-N1: ', length(index.ssp)-N1))
  # p.ssp <- pmin(p.ssp[index.ssp], 1)
  # w.ssp <- 1 / (Y[index.ssp] + (1 - Y[index.ssp]) * p.ssp)
  # print(paste('sum(w.ssp): ', sum(w.ssp)))
  # print(paste('sum(n.ssp * w.ssp / N): ', sum(n.ssp * w.ssp / N)))
  #
  # if (estimate.method == 'LogOddsCorrection'){
  #   offset <- -log(p.ssp)
  # } else (
  #   offset <- NA
  # )


  return (
    list(
      index.ssp = index.ssp,
      w.ssp = w.ssp,
      offset = offset
    )
  )
}
###############################################################################
subsample.estimate <- function(x.ssp,
                               y.ssp,
                               n.ssp,
                               N,
                               w.ssp,
                               offset,
                               beta.plt,
                               estimate.method) {
  if (estimate.method == "Weighted"){
    results.ssp <- glm.coef.estimate(x.ssp, y.ssp, weights = w.ssp)
    beta.ssp <- results.ssp$beta
    P.ssp <- results.ssp$pbeta
    # P.ssp <- pbeta(x.ssp, beta.ssp)
    var.ssp <- results.ssp$cov
    ddL.ssp <- ddL(x.ssp, P.ssp, w = w.ssp * n.ssp / N)
    dL.sq.ssp <- dL.sq(x.ssp, y.ssp, P.ssp, w = w.ssp ^ 2 * n.ssp ^ 2 / N ^2)
  } else if (estimate.method == 'LogOddsCorrection'){
    results.ssp <- glm.coef.estimate(X = x.ssp,
                                     Y = y.ssp,
                                     start = beta.plt,
                                     offset = offset)
    beta.ssp <- results.ssp$beta
    P.ssp <- results.ssp$pbeta # equal to pbeta(x.ssp, beta.ssp, offset)
    var.ssp <- results.ssp$cov
    ddL.ssp <- ddL(x.ssp, P.ssp, w = 1 / n.ssp)
    dL.sq.ssp <- dL.sq(x.ssp, y.ssp, P.ssp, w = 1 / n.ssp ^ 2)
  }

  # var.ssp <- solve(ddL.ssp) %*% (dL.sq.ssp) %*% solve(ddL.ssp)

  return (
    list(
      beta.ssp = beta.ssp,
      ddL.ssp = ddL.ssp,
      dL.sq.ssp = dL.sq.ssp,
      var.ssp = var.ssp
    )
  )
}
###############################################################################
combining <- function(ddL.plt,
                      ddL.ssp,
                      dL.sq.plt,
                      dL.sq.ssp,
                      n.plt,
                      n.ssp,
                      beta.plt,
                      beta.ssp) {
  # print(ddL.plt)
  # print(ddL.ssp)
  ddL.plt <- n.plt * ddL.plt
  ddL.ssp <- n.ssp * ddL.ssp
  # print(ddL.plt)
  # print(ddL.ssp)
  dL.sq.plt <- n.plt ^ 2 * dL.sq.plt
  dL.sq.ssp <- n.ssp ^ 2 * dL.sq.ssp
  MNsolve <- solve(ddL.plt + ddL.ssp)
  # print(MNsolve %*% ddL.plt %*% beta.plt)
  # print(MNsolve %*% ddL.ssp %*% beta.ssp)
  beta.cmb <- c(MNsolve %*% (ddL.plt %*% beta.plt + ddL.ssp %*% beta.ssp))
  var.cmb <- MNsolve %*% (dL.sq.plt + dL.sq.ssp) %*% MNsolve
  return(list(beta.cmb = beta.cmb,
              var.cmb = var.cmb)
  )
}
###############################################################################
rare.logistic.subsampling <- function(formula,
                                      data,
                                      n.plt,
                                      n.ssp,
                                      criterion = c('OptL', 'OptA', 'LCC'),
                                      estimate.method = c('LogOddsCorrection',
                                                          'Weighted', 'Uni'),
                                      alpha = 0.1,
                                      b = 2) {

  model.call <- match.call()
  mf <- model.frame(formula, data)
  Y <- model.response(mf, "numeric") # if we only consider logistic model ,it should not be 'any'.
  X <- model.matrix(formula, mf)
  colnames(X)[1] <- "intercept"
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
    p.plt <- plt.estimate.results$p.plt
    beta.plt <- plt.estimate.results$beta.plt
    ddL.plt <- plt.estimate.results$ddL.plt
    dL.sq.plt <- plt.estimate.results$dL.sq.plt
    ddL.plt.correction <- plt.estimate.results$ddL.plt.correction
    P.plt <- plt.estimate.results$P.plt
    index.plt <- plt.estimate.results$index.plt
    # var.plt <- solve(ddL.plt) %*% dL.sq.plt %*% solve(ddL.plt)

    # subsampling step
    ssp.results <- subsampling(X = X,
                               Y = Y,
                               n.ssp = n.ssp,
                               alpha = alpha,
                               b = b,
                               criterion = criterion,
                               estimate.method = estimate.method,
                               p.plt = p.plt,
                               ddL.plt.correction = ddL.plt.correction,
                               P.plt = P.plt,
                               index.plt = index.plt)
    index.ssp <- ssp.results$index.ssp
    w.ssp <- ssp.results$w.ssp
    offset <- ssp.results$offset

    # subsample estimation step
    ssp.estimate.results <- subsample.estimate(X[index.ssp, ],
                                               Y[index.ssp],
                                               n.ssp = n.ssp,
                                               N = N,
                                               w.ssp = w.ssp,
                                               offset = offset,
                                               beta.plt = beta.plt,
                                               estimate.method = estimate.method)
    beta.ssp <- ssp.estimate.results$beta.ssp
    ddL.ssp <- ssp.estimate.results$ddL.ssp
    dL.sq.ssp <- ssp.estimate.results$dL.sq.ssp
    var.ssp <- ssp.estimate.results$var.ssp

    # combine step
    combining.results <- combining(ddL.plt = ddL.plt,
                                   ddL.ssp = ddL.ssp,
                                   dL.sq.plt = dL.sq.plt,
                                   dL.sq.ssp = dL.sq.ssp,
                                   n.plt = n.plt,
                                   n.ssp = n.ssp,
                                   beta.plt = beta.plt,
                                   beta.ssp = beta.ssp)
    beta.cmb <- combining.results$beta.cmb
    var.cmb <- combining.results$var.cmb
    # prediction.results <- in.sample.prediction(X, Y, beta.cmb)

    return(list(model.call = model.call,
                beta.plt = beta.plt,
                beta.ssp = beta.ssp,
                beta.cmb = beta.cmb,
                var.ssp = var.ssp,
                var.cmb = var.cmb,
                index.plt = index.plt,
                index.ssp = index.ssp
    )
    )
  } else if (estimate.method == "Uni"){
    # Uni:  poisson sampling + unweighted likelihood estimation + correction
    n.uni <- N1 + n.ssp
    pi.uni <- rep(1, N)
    pi.uni[Y == 0] <- n.ssp / N0
    index.uni <- poisson.index(N, pi.uni)
    X.uni <- X[index.uni, ]
    Y.uni <- Y[index.uni]
    results.uni <- glm.coef.estimate(X = X.uni, Y = Y.uni) # not corrected yet
    beta.uni <- results.uni$beta
    # pbeta.uni <- results.uni$pbeta
    var.uni <- results.uni$cov
    beta.uni[1] <- beta.uni[1] + log(n.ssp / N0)

  # else if (estimate.method %in% c("Uni", "UniW")){
  #   # Uni:  negative uniform sampling. poisson sampling + unweighted likelihood estimation + correction
  #   # UniW: poisson sampling + weighted likelihood estimation
  #   n.uni <- N1 + n.ssp
  #   pi.uni <- rep(1, N)
  #   pi.uni[Y == 0] <- n.ssp / N0
  #   index.uni <- poisson.index(N, pi.uni)
  #   X.uni <- X[index.uni, ]
  #   Y.uni <- Y[index.uni]
  #   # offset <- -log(pi.uni[index.uni])
  #   if (estimate.method == "Uni"){
  #     beta.uni <- glm.coef.estimate(X = X.uni, Y = Y.uni)$beta
  #     pbeta <- pbeta(X, beta.uni)
  #     ddL.uni <- MN(X.uni, pbeta[index.uni], 1 / n.uni) # 2nd derivative of the log likelihood function.
  #     dL.sq.uni <- dL.sq(X.uni, Y.uni, pbeta[index.uni], 1 / n.uni^2)
  #     beta.uni[1] <- beta.uni[1] + log(n.ssp / N0) # or log(n.ssp / N0)?
  #   } else if (estimate.method == "UniW"){
  #     w.uni <- 1 / pi.uni[index.uni]
  #     beta.uni <- glm.coef.estimate(X = X.uni, Y = Y.uni, weights = w.uni)$beta
  #     pbeta <- pbeta(X, beta.uni)
  #     ddL.uni <- MN(X.uni, pbeta[index.uni], w.uni / N) # order is wrong
  #     dL.sq.uni <- dL.sq(X.uni, Y.uni, pbeta[index.uni], w.uni^2 / N^2)
  #   }
  #   var.uni <- solve(ddL.uni) %*% dL.sq.uni %*% solve(ddL.uni)
  #   # prediction.results <- in.sample.prediction(X, Y, beta.uni)

    return(list(model.call = model.call,
                index = index.uni,
                beta = beta.uni,
                var = var.uni
                ))
  }
}
