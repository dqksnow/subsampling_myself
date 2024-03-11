# criterion = OptA, OptL, LCC
# sampling.method = WithReplacement, Poisson
# estimate.method = Weighted, LogOddsCorrection
###############################################################################
glm.coef.estimate <- function(X,
                              Y,
                              offset = NULL,
                              start = rep(0, ncol(X)),
                              weights = 1,
                              family) {
  # if (length(weights) == 1) # remove this?
  #   weights <- rep(weights, length(Y))
  family <- switch(family$family.name,
                  "binomial" = binomial(),
                  "poisson" = poisson(),
                  "negative.binomial" = MASS::negative.binomial(2),
                  "gamma" = Gamma(link = "inverse"))
  data <- as.data.frame(cbind(Y, X))
  formula <- as.formula(paste(colnames(data)[1], "~",
                   paste(colnames(data)[-1], collapse = "+"), "-1")) # use '-1' to avoid adding intercept column again.
  # response_var <- colnames(data)[1]
  # predictor_vars <- paste(colnames(data)[-1], collapse = " + ")
  # formula_string <- paste(response_var, "~", predictor_vars, "- 1")
  # formula <- as.formula(formula_string)
  design <- survey::svydesign(ids = ~ 1,
                              weights =  ~ weights,
                              data = data)
  # print(start)
  fit <- ifelse(is.null(offset),
                results <- survey::svyglm(formula = formula,
                                           design = design,
                                           start = start,
                                           family = family),
                results <- survey::svyglm(formula = formula,
                                           design = design,
                                           start = start,
                                           offset = offset,
                                           family = family))
  beta <- results$coefficients
  return(list(beta = beta))
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
  return(list(index.plt = as.vector(index.plt), n.plt.0 = n.plt.0, n.plt.1 = n.plt.1))
}
random.index <- function (N, n, p=NULL) {
  ifelse(is.null(p), index <- sample(N, n, replace = TRUE), index <- sample(N, n, replace = TRUE, prob = p))
  # index <- sample(N, n, replace = TRUE, prob = p)
  return(as.vector(index))
}
poisson.index <- function (N, pi) {
  index <- runif(N) <= pi
  return(as.vector(which(index == TRUE)))
}
poisson.index.test <- function (N, pi) {
  return(which(runif(N) <= pi))
}
# MN <- function (X, pbeta, w) {
#   phi <- pbeta * (1 - pbeta)
#   MN <- t(X) %*% (X * (phi * w))
#   return(MN)
# }
# Psi <- function (X, Y, pbeta, w) {
#   psi <- (Y - pbeta)^2
#   Psi <- t(X) %*% (X * (psi * w))
#   return(Psi)
# }

##############family version
ddL <- function (X, beta, weights = 1, offset = NULL, family) {
  if (family$canonical == TRUE){
    dd.psi <- family$dd.psi(as.vector(X %*% beta), offset)
    MN <- t(X) %*% (X * (dd.psi * weights))
  } else {
    linear.predictor <- as.vector(X %*% beta)
    temp <- family$dd.psi(linear.predictor, offset) * family$d.u(linear.predictor)^2
    MN <- t(X) %*% (X * (temp * weights))
  }
  return(MN)
}
dL.sq <- function (X, beta, Y, weights = 1, offset = NULL, family) {
  linear.predictor <- as.vector(X %*% beta)
  if (family$canonical == TRUE){
    temp <- (Y - family$d.psi(linear.predictor, offset))^2
  } else {
    temp <- (Y - family$d.psi(linear.predictor, offset))^2 * family$d.u(linear.predictor)^2
  }
  Psi <- t(X) %*% (X * (temp * weights))
  return(Psi)
}
##############end
pbeta <- function (X, beta, offset = NULL) {
  ifelse(is.null(offset),
         p <- 1 / (1 + exp(-c(X %*% beta))),
         p <- 1 / (1 + exp(-c(X %*% beta) - offset)))
  return(p)
}
calculate.offset <- function (X,
                              N,
                              dm,
                              d.psi,
                              alpha = alpha,
                              ddL.plt.correction,
                              NPhi = NULL,
                              n.ssp = NULL,
                              criterion,
                              sampling.method) {
  if (criterion == "OptA") {
    norm <- sqrt(rowSums((X %*% t(solve(ddL.plt.correction)))^2))
    nm.1 <- abs(1 - d.psi) * norm
    nm.0 <- abs(d.psi) * norm
  } else if (criterion == "OptL") {
    norm <- sqrt(rowSums(X^2))
    nm.1 <- abs(1 - d.psi) * norm
    nm.0 <- abs(d.psi) * norm
  } else if (criterion == "LCC") {
    nm.1 <- abs(1 - d.psi)
    nm.0 <- abs(d.psi)
  }
  if (sampling.method == 'WithReplacement') {
    stop("Currently only the 'LogOddsCorrection' estimate method with
         'WithReplacement' sampling method has been implemented.")
  } else if (sampling.method == 'Poisson') {
    pi.1 <- pmin(n.ssp * ((1 - alpha) * nm.1 / NPhi + alpha / N), 1)
    pi.0 <- pmin(n.ssp * ((1 - alpha) * nm.0 / NPhi + alpha / N), 1)
  }
  offset <- log(pi.1 / pi.0)
  return(offset)
}
###############################################################################
pilot.estimate <- function(X, Y, n.plt,  family = family){
  N <- nrow(X)
  d <- ncol(X)
  N1 <- sum(Y)
  N0 <- N - N1
  if (family$family.name == 'binomial'){
    ## half half pilot estimator
    halfhalf.index.results <- halfhalf.index(N, Y, n.plt)
    index.plt <- halfhalf.index.results$index.plt
    n.plt.0 <- halfhalf.index.results$n.plt.0
    n.plt.1 <- halfhalf.index.results$n.plt.1
    x.plt <- X[index.plt,]
    y.plt <- Y[index.plt]
    p.plt <- c(rep(1 / (2 * N0), n.plt.0), rep(1 / (2 * N1), n.plt.1))

    # unweighted log likelihood estimation
    beta.plt <- glm.coef.estimate(X = x.plt, Y = y.plt, family = family)$beta # not corrected yet
    # pbeta.plt <- pbeta(x.plt, beta.plt)
    # ddL.plt <- MN(x.plt, pbeta.plt, 1 / n.plt) # 2nd derivative of the log likelihood function.
    # Psi.plt <- Psi(x.plt, y.plt, pbeta.plt, 1 / n.plt^2)
    #family version
    ddL.plt <- ddL(x.plt, beta.plt, weights = 1 / n.plt, family = family)
    dL.sq.plt <- dL.sq(x.plt, beta.plt, y.plt, weights = 1 / n.plt^2, family = family)

    # after correction
    beta.plt[1] <- beta.plt[1] - log(N0 / N1)
    d.psi <- family$d.psi(X %*% beta.plt)
    ddL.plt.correction <- ddL(x.plt, beta.plt, w = 1 / n.plt, family = family)
  } else {
    index.plt <- random.index(N, n.plt)
    x.plt <- X[index.plt,]
    y.plt <- Y[index.plt]
    p.plt <- rep(1 / N, n.plt)
    # unweighted log likelihood estimation
    beta.plt <- glm.coef.estimate(X = x.plt, Y = y.plt, family = family)$beta
    ddL.plt <- ddL.plt.correction <- ddL(x.plt, beta.plt, weights = 1 / n.plt, family = family)
    dL.sq.plt <- dL.sq(x.plt, beta.plt, y.plt, weights = 1 / n.plt^2, family = family)
    d.psi <- family$d.psi(X %*% beta.plt)

  }

  var.plt <- solve(ddL.plt) %*% dL.sq.plt %*% solve(ddL.plt) #useless

  # # case control pilot estimator
  # # if sampling uniformly, p.plt <- rep(1 / N, N), and correction is unnecessary
  # p.plt <- rep(1 / (2 * N1), N)
  # p.plt[which(Y == 0)] <- 1 / (2 * N0)
  # random.index.results <- random.index(N, n.plt, p.plt)
  # index.plt <- random.index.results$index
  # x.plt <- X[index.plt,]
  # y.plt <- Y[index.plt]
  # p.plt <- p.plt[index.plt]

  # # # uniform pilot estimator
  # p.plt <- rep(1 / N, N)
  # index.plt <- random.index(N, n.plt, p.plt)
  # x.plt <- X[index.plt,]
  # y.plt <- Y[index.plt]
  # p.plt <- p.plt[index.plt]
  # # unweighted log likelihood estimation
  # beta.plt <- glm.coef.estimate(X = x.plt, Y = y.plt)$beta # not corrected yet
  # pbeta.plt <- pbeta(x.plt, beta.plt)
  # ddL.plt <- MN(x.plt, pbeta.plt, 1 / n.plt) # 2nd derivative of the log likelihood function.
  # Psi.plt <- Psi(x.plt, y.plt, pbeta.plt, 1 / n.plt^2)
  # P.plt <- pbeta(X, beta.plt)
  # MN.plt <- MN(x.plt, P.plt[index.plt], 1 / n.plt) # the pilot estimator of MN, used for calculating 'OptA' subsampling probability
  # var.plt <- solve(ddL.plt) %*% Psi.plt %*% solve(ddL.plt)


  return(
    list(
      p.plt = p.plt,
      beta.plt = beta.plt,
      ddL.plt = ddL.plt,
      dL.sq.plt = dL.sq.plt,
      ddL.plt.correction = ddL.plt.correction,
      d.psi = d.psi,
      index.plt = index.plt,
      var.plt = var.plt
    )
  )
}

calculate.nm <- function(X, Y, ddL.plt.correction, d.psi, criterion){
  if (criterion == "OptA"){
    nm <- sqrt(rowSums((X %*% t(solve(ddL.plt.correction)))^2))
    nm <- abs(Y - d.psi) * nm # numerator
  } else if (criterion == "OptL"){
    nm <- sqrt(rowSums(X^2))
    nm <- abs(Y - d.psi) * nm
  } else if (criterion == "LCC"){
    nm <- abs(Y - d.psi)
  }
  return(nm)
}

###############################################################################
subsampling <- function(X, Y, n.ssp, alpha, b, criterion, estimate.method, sampling.method,
                        p.plt, ddL.plt.correction, d.psi, index.plt) {

  N <- nrow(X)
  N1 <- sum(Y)
  N0 <- N - N1
  n.plt <- length(index.plt) # length(index.plt) might be smaller than the n.plt user sets, so here we reset n.plt.
  w.ssp <- offset <- NA
  # microbenchmark(
  #   calculate.norm(X, Y, ddL.plt.correction, d.psi, criterion),
  #   fun2(x)
  # )

  nm <- calculate.nm(X, Y, ddL.plt.correction, d.psi, criterion) # numerator

  # if (criterion == "OptA"){
  #   norm <- sqrt(colSums(solve(ddL.plt.correction, t(X))^2)) # the norm term of the criterion
  #   nm <- abs(Y - d.psi) * norm # numerator
  # } else if (criterion == "OptL"){
  #   norm <- sqrt(rowSums(X^2))
  #   nm <- abs(Y - d.psi) * norm
  # } else if (criterion == "LCC"){
  #   nm <- abs(Y - d.psi)
  # }

  if (sampling.method == "WithReplacement"){
    dm <- sum(nm) # denominator
    # boxplot(nm / dm)
    p.ssp <- (1 - alpha) * nm / dm + alpha / N
    # boxplot(p.ssp)
    index.ssp <- random.index(N, n.ssp, p.ssp)

    if (estimate.method == 'LogOddsCorrection') {
      stop("Currently only the 'LogOddsCorrection' estimate method with
         'WithReplacement' sampling method has been implemented.")
      offset <- calculate.offset(criterion = criterion,
                                 sampling.method = sampling.method,
                                 d.psi = d.psi[index.ssp],
                                 dm = dm,
                                 X = X[index.ssp,],
                                 ddL.plt.correction = ddL.plt.correction
                                 )
    } else if (estimate.method == 'Weighted'){
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
    # mbm <- microbenchmark::microbenchmark(
    #   poisson.index(N, p.ssp),
    #   poisson.index.test(N, p.ssp),
    #   unit = 'us'
    # )
    # print(mbm)
    if (estimate.method == 'LogOddsCorrection') {
      offset <- calculate.offset(X = X[index.ssp,],
                                 N = N,
                                 dm = dm,
                                 d.psi = d.psi[index.ssp],
                                 alpha = alpha,
                                 ddL.plt.correction = ddL.plt.correction,
                                 criterion = criterion,
                                 sampling.method = sampling.method,
                                 NPhi = NPhi,
                                 n.ssp = n.ssp)
    } else if (estimate.method == 'Weighted'){
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
                                beta.plt, sampling.method, estimate.method, family) {
  if (estimate.method == "Weighted") {
    beta.ssp <- glm.coef.estimate(x.ssp, y.ssp, weights = w.ssp, family = family)$beta
    # print(beta.ssp)
    # P.ssp <- pbeta(x.ssp, beta.ssp)
    # d.psi <- family$d.psi(linear.predictor.ssp)
    if (sampling.method == 'Poisson') {
      ddL.ssp <- ddL(x.ssp, beta.ssp, y.ssp, weights = w.ssp / N, family = family)
      # ddL.ssp <- MN(x.ssp, P.ssp, w.ssp / N) # 2nd derivative of the log likelihood function.
      dL.sq.ssp <- dL.sq(x.ssp, beta.ssp, y.ssp, weights = w.ssp^2 / N^2, family = family) # the estimator of the variance of 1st derivative
      Lambda.ssp <- 0 # holding
      # print(ddL.ssp)
    } else if (sampling.method == "WithReplacement") {
      ddL.ssp <- ddL(x.ssp, beta.ssp, y.ssp, weights = w.ssp / (N * n.ssp), family = family)
      dL.sq.ssp <- dL.sq(x.ssp, beta.ssp, y.ssp, weights = w.ssp ^ 2 / (N ^ 2 * n.ssp ^ 2), family = family)
      # ddL.ssp <- MN(x.ssp, P.ssp, w.ssp / (N * n.ssp))
      # Psi.ssp <- Psi(x.ssp, y.ssp, P.ssp, w.ssp ^ 2 / (N ^ 2 * n.ssp ^ 2))
      c <- n.ssp / N
      Lambda.ssp <- c * dL.sq(x.ssp, beta.ssp, y.ssp, weights = w.ssp / (N * n.ssp ^ 2), family = family)
      # Lambda.ssp <- c * Psi(x.ssp, y.ssp, P.ssp, w.ssp / (N * n.ssp ^ 2))
      # print(Lambda.ssp / Psi.ssp)
    }
  } else if (estimate.method == 'LogOddsCorrection') {
    # print(paste('beta.plt:', beta.plt))
    beta.ssp <- glm.coef.estimate(X = x.ssp,
                                       Y = y.ssp,
                                       start = beta.plt,
                                       offset = offset,
                                       family = family)$beta
    # P.ssp  <- pbeta(x.ssp, beta.ssp, offset)
    ddL.ssp <- ddL(x.ssp, beta.ssp, weights = 1 / n.ssp, offset= offset, family = family)
    dL.sq.ssp <- dL.sq(x.ssp, beta.ssp, y.ssp, weights = 1 / n.ssp ^ 2, offset = offset, family = family)
    # ddL.ssp <- MN(x.ssp, P.ssp, 1 / n.ssp)
    # Psi.ssp <- Psi(x.ssp, y.ssp, P.ssp, 1 / n.ssp ^ 2)
    Lambda.ssp <- 0 # holding
  }

  var.ssp <- solve(ddL.ssp) %*% dL.sq.ssp %*% solve(ddL.ssp)
  var.ssp.true <- solve(ddL.ssp) %*% (dL.sq.ssp + Lambda.ssp) %*% solve(ddL.ssp)

  return(list(beta.ssp = beta.ssp,
               ddL.ssp = ddL.ssp,
              dL.sq.ssp = dL.sq.ssp,
               var.ssp = var.ssp,
               var.ssp.true = var.ssp.true,
               Lambda.ssp = Lambda.ssp)
  )
}
###############################################################################
combining <- function(ddL.plt, ddL.ssp, dL.sq.plt, dL.sq.ssp, Lambda.ssp, n.plt, n.ssp, beta.plt, beta.ssp) {
  ddL.plt <- n.plt * ddL.plt
  ddL.ssp <- n.ssp * ddL.ssp

  # print(ddL.ssp)
  Psi.plt <- n.plt ^ 2 * dL.sq.plt
  Psi.ssp <- n.ssp ^ 2 * dL.sq.ssp

  # print(ddL.ssp)
  # print(Psi.ssp)

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
in.sample.prediction <- function(X, Y, beta){
  N <- length(Y)
  pbeta.estimate <- pbeta(X, beta)
  Y.pred <- rep(0, N)
  Y.pred[which(pbeta.estimate > 0.5)] <- 1
  plot(pbeta.estimate)
  # Calculate True Positives, True Negatives, False Positives, False Negatives
  tp <- sum(Y == 1 & Y.pred == 1)
  tn <- sum(Y == 0 & Y.pred == 0)
  fp <- sum(Y == 0 & Y.pred == 1)
  fn <- sum(Y == 1 & Y.pred == 0)
  print(sum(Y.pred == 1))
  # accuracy measures the model's overall performance of prediction.
  accuracy <- (tp + tn) / N
  # precision measures the model's performance when it predicts positive outcomes(Y=1).
  precision <- tp / (tp + fp)
  # recall measures that of all the actual positive (Y=1) instances, how many did the model correctly predict as positive.
  recall <- tp / (tp + fn) # tp + fn = N1
  # f1.score is the harmonic mean of precision and recall.
  f1.score <- 2 * (precision * recall) / (precision + recall)
  # roc_obj <- roc(Y, pbeta.estimate)
  # print(auc(roc_obj))
  # plot(roc_obj, main = "ROC Curve", print.auc = TRUE)

  return(list(
    Accuracy = accuracy,
    Precision = precision,
    Recall = recall,
    F1.Score = f1.score
  ))
}

###############################################################################
glm.optimal.subsampling <- function(formula, data, n.plt, n.ssp,
                                         family = c('binomial','poisson','negative.binomial','gamma'),
                                         criterion = c('OptL', 'OptA', 'LCC'),
                                         sampling.method = c('Poisson', 'WithReplacement'),
                                         estimate.method = c('LogOddsCorrection', 'Weighted'),
                                         # MN_custom, # holding
                                         alpha = 0,
                                         b = 2) {
  model.call <- match.call()
  # family <- match.arg(family)

  if(missing(family)) {
    stop("Specify a valid 'family' from
         c('binomial','poisson','negative.binomial','gamma')")
  }
  if(is.function(family)) family <- family$family
  family <- switch(family,
                   "binomial" = binomial.expand(),
                   "poisson" = poisson.expand(),
                   "negative.binomial" = negative.binomial.expand(),
                   "gamma" = gamma.expand())
  # criterion <- match.arg(criterion)
  # sampling.method <- match.arg(sampling.method)
  # estimate.method <- match.arg(estimate.method)

  mf <- model.frame(formula, data)
  # Extract the response and predictor variables from the model frame
  Y <- model.response(mf, "any")
  # if (is.character(Y) && length(unique(Y)) == 2) {
  #   levels <- unique(Y)
  #   Y <- as.integer(Y == levels[2])  # Assuming levels[2] is the 'success' category
  # }
  X <- model.matrix(formula, mf)
  colnames(X)[1] <- "intercept"
  # print(X)
  N <- nrow(X)
  d <- ncol(X)
  N1 <- sum(Y)
  N0 <- N - N1



  if (estimate.method %in% c("Weighted", "LogOddsCorrection")) {

    # pilot step
    plt.estimate.results <- pilot.estimate(X = X, Y = Y, n.plt = n.plt, family = family)
    p.plt <- plt.estimate.results$p.plt
    beta.plt <- plt.estimate.results$beta.plt
    ddL.plt <- plt.estimate.results$ddL.plt
    dL.sq.plt <- plt.estimate.results$dL.sq.plt
    ddL.plt.correction <- plt.estimate.results$ddL.plt.correction
    d.psi <- plt.estimate.results$d.psi
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
                               ddL.plt.correction = ddL.plt.correction,
                               d.psi = d.psi,
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
                                              estimate.method = estimate.method,
                                              family = family)
    beta.ssp <- ssp.estimate.results$beta.ssp
    ddL.ssp <- ssp.estimate.results$ddL.ssp
    dL.sq.ssp <- ssp.estimate.results$dL.sq.ssp
    Lambda.ssp <- ssp.estimate.results$Lambda.ssp
    var.ssp <- ssp.estimate.results$var.ssp
    var.ssp.true <- ssp.estimate.results$var.ssp.true

    # combining step
    combining.results <- combining(ddL.plt = ddL.plt,
                                   ddL.ssp = ddL.ssp,
                                   dL.sq.plt = dL.sq.plt,
                                   dL.sq.ssp = dL.sq.ssp,
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
    index.uni <- random.index(N, n.uni)
    # index.uni <- sample(N, n.uni, replace = F)
    # index.uni <- poisson.index(N, n.uni/N)
    X.uni <- X[index.uni, ]
    Y.uni = Y[index.uni]
    beta.uni <- glm.coef.estimate(X = X.uni, Y = Y.uni, family = family)$beta

    # pbeta.uni <- pbeta(X[index.uni, ], beta.uni)
    # ddL.uni <- MN(X.uni, pbeta.uni, 1 / n.uni) # 2nd derivative of the log likelihood function.
    # Psi.uni <- Psi(X.uni, Y.uni, pbeta.uni, 1 / n.uni^2)
    ddL.uni <- ddL(X.uni, beta.uni, Y.uni, weights = 1 / n.uni, family = family)
    dL.sq.plt <- dL.sq(X.uni, beta.uni, Y.uni, weights = 1 / n.uni^2, family = family)
    var.uni <- solve(ddL.uni) %*% dL.sq.plt %*% solve(ddL.uni)
    # var.uni <- solve(ddL.uni) %*% (dL.sq.plt * (1+n.uni/N)) %*% solve(ddL.uni)
    # prediction.results <- in.sample.prediction(X, Y, beta.uni)

    return(list(model.call = model.call,
                index = index.uni,
                beta = beta.uni,
                var = var.uni))
  }
}
