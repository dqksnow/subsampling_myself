# criterion = optA, optL, LCC
# sampling.method = WithReplacement, Poisson
# estimate.method = Weighted, LogOddsCorrection
getMSLE <- function(X, Y, offset, start) {
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
logistic.coef.estimate <-  function( X,
                                     Y,
                                     offset = rep(0, length(Y)),
                                     start = rep(0, ncol(X)),
                                     weights = rep(1, length(Y))) {
  if (length(weights) == 1)
    weights <- rep(weights, length(Y))
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
uni.index <- function(N, Y, n.plt){
  N1 <- sum(Y)
  N0 <- N - N1
  if (N0 < n.plt / 2 | N1 < n.plt / 2){ # if it happens, sampling would be unnecessary
    warning(paste("n.plt/2 exceeds the number of Y=1 or Y=0 in the full data.",
                  "In this case, all rare events will be drawn into the pilot sample."))
  }
  n.plt.0 <- min(N0, n.plt / 2)
  n.plt.1 <- min(N1, n.plt / 2)
  index.plt <- c(sample(which(Y == 0), n.plt.0, replace = FALSE),
                 sample(which(Y == 1), n.plt.1, replace = FALSE))
  return (list(index.plt = index.plt, n.plt.0 = n.plt.0, n.plt.1 = n.plt.1))
}
random.index <- function(N, n, p){
  index <- sample(N, n, replace = TRUE, prob = p)
  return (list(index = index))
}
poisson.index <- function(N, pi){
  index <- runif(N) <= pi
  return(which(index == TRUE))
}
MN <- function(X, pbeta, w){
  phi <- pbeta * (1 - pbeta)
  MN <- t(X) %*% (X * (phi * w))
  return(MN)
}
Psi <- function(X, Y, pbeta, w){
  psi <- (Y - pbeta)^2
  Psi <- t(X) %*% (X * (psi * w^2))
  return(Psi)
}
pbeta <- function(X, beta, offset = NA){
  ifelse(is.na(offset),
         p <- 1 / (1 + exp(-c(X %*% beta))),
         p <- 1 / (1 + exp(-c(X %*% beta) - offset)))
  return (p)
}
calculate.offset.nm <- function(criterion, P.plt, norm = NULL){
  if (criterion %in% c("optA", "optL")){
    nm.1 <- abs(1 - P.plt) * norm
    nm.0 <- abs(P.plt) * norm
  } else if (criterion == "LCC"){
    nm.1 <-  abs(1 - P.plt)
    nm.0 <-  abs(P.plt)
  }
  return(list(nm.1 = nm.1, nm.0 = nm.0))
}
###############################################################################
pilot.estimate <- function(X, Y, n.plt){
  N <- nrow(X)
  d <- ncol(X)
  N1 <- sum(Y)
  N0 <- N - N1

  ## half half pilot estimator
  # uni.index.results <- uni.index(N, Y, n.plt)
  # index.plt <- uni.index.results$index.plt
  # n.plt.0 <- uni.index.results$n.plt.0
  # n.plt.1 <- uni.index.results$n.plt.1
  # x.plt <- X[index.plt,]
  # y.plt <- Y[index.plt]


  ## case control pilot estimator
  # if sampling uniformly, p.plt <- rep(1 / N, N), and correction is unnecessary
  p.plt <- rep(1 / (2 * N1), N)
  p.plt[which(Y == 0)] <- 1 / (2 * N0)
  random.index.results <- random.index(N, n.plt, p.plt)
  index.plt <- random.index.results$index
  x.plt <- X[index.plt,]
  y.plt <- Y[index.plt]

  # unweighted log likelihood function
  beta.plt <- logistic.coef.estimate(X = x.plt, Y = y.plt)$beta
  pbeta.plt <- pbeta(x.plt, beta.plt)
  ddL.plt <- MN(x.plt, pbeta.plt, 1 / n.plt) # 2nd derivative of the log likelihood function.
  Psi.plt <- Psi(x.plt, y.plt, pbeta.plt, 1 / n.plt)
  beta.plt[1] <- beta.plt[1] - log(N0 / N1)
  P.plt <- pbeta(X, beta.plt)
  MN.plt <- MN(x.plt, P.plt[index.plt], 1 / n.plt) # the estimator of MN, used for calculating 'optA' subsampling probability

  return(list(
    p.plt = p.plt[index.plt],
    beta.plt = beta.plt,
    ddL.plt = ddL.plt,
    Psi.plt = Psi.plt,
    MN.plt = MN.plt,
    P.plt = P.plt,
    index.plt = index.plt
  ))
}
###############################################################################
subsampling <- function(X, Y, n.ssp, alpha, b, criterion, sampling.method,
                        p.plt, MN.plt, P.plt, index.plt){

  N <- nrow(X)
  N1 <- sum(Y)
  N0 <- N - N1
  n.plt <- length(index.plt) # length(index.plt) might be smaller than the n.plt user sets.

  if (criterion == "optA"){
    norm <- sqrt(colSums(solve(MN.plt, t(X))^2))
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
    index.ssp <- sample(1:N, n.ssp, replace = T, prob = p.ssp)
    w.ssp <- 1 / p.ssp[index.ssp]
    offset.nm <- calculate.offset.nm(criterion = criterion, P.plt = P.plt[index.ssp], norm = norm[index.ssp])
    pi.1 <- offset.nm$nm.1 / dm
    pi.0 <- offset.nm$nm.0 / dm
    offset.cmb <- c(rep(log(N0/N1), n.plt), log(pi.1 / pi.0))
  } else if (sampling.method == "Poisson"){
    NPhi <- sum(nm[index.plt] / p.plt) / n.plt
    p.ssp <- n.ssp * ((1 - alpha) * nm / NPhi + alpha / N)
    index.ssp <- poisson.index(N, p.ssp)
    w.ssp <- 1 / pmin(p.ssp[index.ssp], 1)
    offset.nm <- calculate.offset.nm(criterion = criterion, P.plt = P.plt[index.ssp], norm = norm[index.ssp])
    pi.1 <- pmin((n.ssp) * offset.nm$nm.1 / NPhi, 1)
    pi.0 <- pmin((n.ssp) * offset.nm$nm.0 / NPhi, 1)
    offset.cmb <- c(rep(log(N0/N1), n.plt), log(pi.1 / pi.0))
  }

  return (list(index.ssp = index.ssp,
               offset.cmb = offset.cmb,
               w.ssp = w.ssp))
}
###############################################################################
CombineSubsample.estimate <- function(x.cmb, y.cmb, n.cmb,
                                       w.cmb, offset.cmb,
                                       beta.plt, sampling.method, estimate.method) {
  if (estimate.method == "Weighted"){
    beta.cmb <- logistic.coef.estimate(x.cmb, y.cmb, weights = w.cmb)$beta # get similar results with ‘start’
    beta.cmb <- unlist(beta.cmb)
  } else if (estimate.method == 'LogOddsCorrection'){
    beta.cmb <- logistic.coef.estimate(X = x.cmb,
                                       Y = y.cmb,
                                       start = beta.plt,
                                       offset = offset.cmb)$beta
    beta.cmb <- unlist(beta.cmb)
  }
  return (
    list(
      beta.cmb = beta.cmb
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

  if (estimate.method %in% c("Weighted", "LogOddsCorrection")){

    # pilot step
    plt.estimate.results <- pilot.estimate(X = X, Y = Y, n.plt = n.plt)
    p.plt <- plt.estimate.results$p.plt
    beta.plt <- plt.estimate.results$beta.plt
    ddL.plt <- plt.estimate.results$ddL.plt
    Psi.plt <- plt.estimate.results$Psi.plt
    MN.plt <- plt.estimate.results$MN.plt
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
                               p.plt = p.plt,
                               MN.plt = MN.plt,
                               P.plt = P.plt,
                               index.plt = index.plt)
    index.ssp <- ssp.results$index.ssp
    w.ssp <- ssp.results$w.ssp
    offset.cmb <- ssp.results$offset.cmb
    index.cmb <- c(index.plt,index.ssp)
    if (sampling.method == 'WithReplacement'){
      w.cmb <- c(1/(p.plt), w.ssp)
    } else {w.cmb <- c(1/(n.plt*p.plt), w.ssp)}
    # w.cmb <- ifelse(sampling.method == 'WithReplacement', c(1/(p.plt), w.ssp), c(1/(n.plt*p.plt), w.ssp))
    print(w.cmb)
    # print(offset.cmb)
    cmb.estimate.reults <- CombineSubsample.estimate(X[index.cmb, ],
                                                     Y[index.cmb],
                                                     n.cmb = n.plt + n.ssp,
                                                     w.cmb = w.cmb,
                                                     offset.cmb = offset.cmb,
                                                     beta.plt = beta.plt,
                                                     sampling.method = sampling.method,
                                                     estimate.method = estimate.method)
    beta.cmb <- cmb.estimate.reults$beta.cmb
    print(beta.cmb)
    # sample combine test

    return(list(beta.plt = beta.plt,
                beta.cmb = beta.cmb,
                var = matrix(0,d,d),
                index.plt = index.plt,
                index.ssp = index.ssp))
  } else if (estimate.method == "Uni"){
    random.index.results <- random.index(N, n.plt + n.ssp, rep(1/N, N))
    index.ssp <- random.index.results$index
    beta.ssp <- logistic.coef.estimate(X = X[index.ssp,], Y = Y[index.ssp])$beta
    return(list(index.ssp = index.ssp,
                beta.ssp = beta.ssp))
  }
}
