proptional_ssp <- function(N, K, Y){
  prob <- 1/((K+1) * as.data.frame(table(Y))$Freq)
  return(prob[Y+1])
}
###############################################################################
random.index <- function (N, n, p = NULL) {
  ifelse(is.null(p),
         index <- sample(N, n, replace = TRUE),
         index <- sample(N, n, replace = TRUE, prob = p))
  return(as.vector(index))
}
###############################################################################
poisson.index <- function (N, pi) {
  return(which(runif(N) <= pi))
}
###############################################################################
softmax.coef.estimate <- function(X, y, weights = NULL, offset = NULL){
  if (is.null(offset)) {
    ifelse(is.null(weights),
           fit <- nnet::multinom(y ~ X - 1, trace = FALSE),
           fit <- nnet::multinom(y ~ X - 1, weights = weights, trace = FALSE))
  } else {
    fit <- nnet::multinom(y ~ X - 1 + offset(offset), trace = FALSE)
  }
  return(list(beta = t(coef(fit)), # return a d*K matrix
              # P = fit$fitted.values[, -1], # return a n*K matrix
              P1 = fit$fitted.values,
              var = sum((summary(fit)$standard.errors)^2)) # return a n*(K+1) matrix
         )
}
###############################################################################
pbeta.multi <- function(X, beta){
  P <- exp(X %*% beta - matrixStats::rowLogSumExps(cbind(0, X %*% beta)))
  P1 <- cbind((1 - rowSums(P)), P)
  return(P1)
}
########unefficient version####################################################
# softmax.ddL.old <- function(X, P, p, K, d, scale){
#   X.Kd <- X[, rep(1:d, K)]
#   P.kd <- P[, rep(1:K, each = d)]
#   spi <- sqrt(p)
#   diagpX.pd <- (P.kd * X.Kd) / spi
#   MN.p1 <- t(X.Kd / spi) %*% diagpX.pd
#   MN.p2 <- t(diagpX.pd) %*% diagpX.pd
#   MN <- MN.p1 * (diag(1, K) %x% matrix(1, nrow = d, ncol = d)) - MN.p2
#   MN <- MN / scale
#   return(MN)
# }
###############################################################################
# softmax.ddL.test1 <- function(X, P, p, K, d, scale){
#   ddL <- matrix(0, nrow = K * d, ncol = K * d)
#   p.sqrt <- sqrt(p)
#   for (i in 1:K){
#     XPi.t <- t(X * (P[, i] / p.sqrt))
#     for (j in i:K){
#       XPj <- X * (P[, j] / p.sqrt)
#       ddL[(1+(i-1)*d):(i*d), (1+(j-1)*d):(j*d)] <-
#         ddL[(1+(j-1)*d):(j*d), (1+(i-1)*d):(i*d)] <- -XPi.t %*% XPj
#     }
#     ddL[(1+(i-1)*d):(i*d), (1+(i-1)*d):(i*d)] <- ddL[(1+(i-1)*d):(i*d), (1+(i-1)*d):(i*d)] +
#       XPi.t %*% (X / p.sqrt)
#   }
#   ddL <- ddL / scale
#   return(ddL)
# }
###############################################################################
softmax.ddL <- function(X, P, p, K, d, scale){
  ddL <- matrix(0, nrow = K * d, ncol = K * d)
  for (i in 1:K){
    for (j in i:K){
      ifelse(j==i,
             XPj <- t( X * ((P[, i] - P[, i]*P[, j]) / p) ),
             XPj <- t( X * ((- P[, i]*P[, j]) / p) )
             )
      ddL[(1+(i-1)*d):(i*d), (1+(j-1)*d):(j*d)] <-
        ddL[(1+(j-1)*d):(j*d), (1+(i-1)*d):(i*d)] <- XPj %*% X
    }
  }
  ddL <- ddL / scale
  return(ddL)
}
###############################################################################
softmax.dL.sq <- function(X, Y.matrix, P, p, K, d, scale){
  S <- Y.matrix - P
  dL <- matrix(0, nrow = K * d, ncol = K * d)
  X.t <- t(X)
  for (i in 1:K){
    for (j in i:K){
      XSij <- X * (S[, i] * S[, j] / p^2)
      dL[(1+(i-1)*d):(i*d), (1+(j-1)*d):(j*d)] <-
        dL[(1+(j-1)*d):(j*d), (1+(i-1)*d):(i*d)] <- X.t %*% XSij
    }
  }
  dL <- dL / scale
  return(dL)
}
#######unefficient version#####################################################
# softmax.dL.sq.old <- function(X, Y.matrix, P, p, K, d, scale){
#   X.Kd <- X[, rep(1:d, K)]
#   psi <- Y.matrix[, rep(1:K, each = d)] - P[, rep(1:K, each = d)]
#   spi <- p
#   psiX <- psi * X.Kd / spi
#   PsiX <- t(psiX) %*% (psiX) / scale
#   return(PsiX)
# }
###############################################################################
# Omega <- function(X, P1, p, K, d, scale){
#   Omega <- matrix(0, nrow = K * d, ncol = K * d)
#   P.sq <- rowSums(P1^2) # N * 1
#   P0 <- P1[, -1]
#   X.t <- t(X)
#   for (i in 1:K){
#     for (j in i:K){
#       if (j == i){
#         XP <- X * ((P0[, i] * P0[, j] * P.sq - P0[, i]^2 * P0[, j] -
#                       P0[, i] * P0[, j]^2 + P0[, i]^2) / p)
#       } else {
#         XP <- X * ((P0[, i] * P0[, j] * P.sq - P0[, i]^2 * P0[, j] -
#                       P0[, i] * P0[, j]^2 ) / p)
#       }
#       Omega[(1+(i-1)*d):(i*d), (1+(j-1)*d):(j*d)] <-
#         Omega[(1+(j-1)*d):(j*d), (1+(i-1)*d):(i*d)] <- X.t %*% XP
#     }
#   }
#   Omega <- Omega / scale
#   return(Omega)
# }
###############################################################################
Omega <- function(X, P1, p, K, d, scale){
  Omega <- matrix(0, nrow = K * d, ncol = K * d)
  P.sq <- rowSums(P1^2) # N * 1
  P0 <- P1[, -1]
  X.t <- t(X)
  for (i in 1:K){
    for (j in i:K){
      if (j == i){
        XP <- X * (((P.sq + 1 - 2 * P0[, i]) * P0[, i]^2) / p)
        # XP <- X * ((P0[, i] * P0[, j] * P.sq - P0[, i]^2 * P0[, j] -
        #               P0[, i] * P0[, j]^2 + P0[, i]^2) / p)
      } else {
        XP <- X * ((P.sq - P0[, i] - P0[, j]) * (P0[, i] * P0[, j])  / p)
        # XP <- X * ((P0[, i] * P0[, j] * P.sq - P0[, i]^2 * P0[, j] -
        #               P0[, i] * P0[, j]^2 ) / p)
      }
      Omega[(1+(i-1)*d):(i*d), (1+(j-1)*d):(j*d)] <-
        Omega[(1+(j-1)*d):(j*d), (1+(i-1)*d):(i*d)] <- X.t %*% XP
    }
  }
  Omega <- Omega / scale
  return(Omega)
}
###############################################################################
# Omega.test1 <- function(X, P1, p, K, d, scale){
#   Omega <- matrix(0, nrow = K * d, ncol = K * d)
#   p.sqrt <- sqrt(p)
#   P.sq <- rowSums(P1^2) # N * 1
#   P0 <- P1[, -1]
#   for (i in 1:K){
#     XP0i1 <- X * (P0[, i] / p.sqrt)
#     XP0i2 <- X * ((P0[, i]^2) / p.sqrt)
#     XP0i3 <- X * (P0[, i] * P.sq / p.sqrt)
#     for (j in i:K){
#       XP0j <- X * (P0[, j] / p.sqrt)
#       Omega[(1+(i-1)*d):(i*d), (1+(j-1)*d):(j*d)] <-
#         Omega[(1+(j-1)*d):(j*d), (1+(i-1)*d):(i*d)] <- (t(XP0i3) %*% XP0j) - 2 * (t(XP0j) %*% XP0i2)
#     }
#     Omega[(1+(i-1)*d):(i*d), (1+(i-1)*d):(i*d)] <- Omega[(1+(i-1)*d):(i*d), (1+(i-1)*d):(i*d)] +
#        t(XP0i1) %*% XP0i1
#   }
#   # Omega <- Omega / scale
#   return(Omega)
# }
# ###############################################################################
# Omega.test <- function(X, K, d, n, P1, p){
#   Omega <- matrix(0, nrow = K*d, ncol = K*d)
#   P0 <- P1[, -1]
#   for(i in 1 : n){
#     B3 <- matrix(t(P1[i, ]) %*% P1[i, ], nrow=K, ncol = K)
#     BIBI <- (B3 - 2 * diag(P0[i,])) %*% (P0[i,] %*% t(P0[i,])) +
#             diag(P0[i, ]^2)
#     Omega <- Omega + BIBI %x% (X[i, ] %*% t(X[i, ])) / p[i]
#   }
#   # Omega <- Omega/ (n * N)
#   return(Omega)
# }
###############################################################################
# Omega.test <- function(X, K, d, n, P1, p){
#   Omega <- matrix(0, nrow = K*d, ncol = K*d)
#   P0 <- P1[, 1]
#   PK <- P1[, -1]
#   for(i in 1 : n){
#     BIBI <- (t(P1[i, ] ) %*% P1[i, ]) * (PK[i, ] %*% t(PK[i, ])) + diag(PK[i, ]^2)
#     - 2 * diag(PK[i, ]) %*% (PK[i, ] %*% t(PK[i, ]))
#     Omega <- Omega + BIBI %x% (X[i, ] %*% t(X[i, ])) / p[i]
#   }
#   # Omega <- Omega/ (n * N)
#   return(Omega)
# }
###############################################################################
softmax.plt.estimate <- function(X, Y, Y.matrix, n.plt, N, K, d, criterion){
  # set.seed(3)
  ### uniform sample with rep
  index.plt <- random.index(N, n.plt)
  p.plt <- rep(1 / N, n.plt) # plt sampling probability
  # ### uniform poisson
  # index.plt <- poisson.index(N, n.plt/N)
  # n.plt <- length(index.plt)
  # p.plt <- rep(n.plt / N, length(index.plt)) # plt sampling probability
  ###
  x.plt <- X[index.plt,]
  y.plt <- Y[index.plt]
  results <- softmax.coef.estimate(x.plt, y.plt)
  beta.plt <- results$beta
  P1.plt <-  results$P1 # n.plt*(K+1) matrix

  ddL.plt <- softmax.ddL(X = x.plt, P = P1.plt[, -1], p = p.plt, K, d, scale = N * n.plt) # n.plt
  dL.sq.plt <- softmax.dL.sq(X = x.plt, Y.matrix = Y.matrix[index.plt, ],
                             P = P1.plt[, -1], p = p.plt, K = K, d = d, scale = (N^2 * n.plt^2)) # (n.plt^2)
  c <- n.plt / N
  Lambda.plt <- c * softmax.dL.sq(X = x.plt, Y.matrix = Y.matrix[index.plt, ],
                                  P = P1.plt[, -1], p = p.plt,
                                  K = K, d = d, scale = N^2*n.plt^2)
  # cov.plt <- solve(ddL.plt) %*% (dL.sq.plt * (1 + c)) %*% solve(ddL.plt)
  # Lambda.plt <- 0
  cov.plt <- solve(ddL.plt) %*% (dL.sq.plt + Lambda.plt) %*% solve(ddL.plt)

  if (criterion == "MSPE") {
    Omega.plt <- Omega(x.plt, P1 = P1.plt, p = 1, K, d, scale = n.plt) # n.plt
  } else {Omega.plt <- NA}
  P1.plt.N <- pbeta.multi(X, beta.plt) # N*(K+1)


  ### case control
  # p.plt <- proptional_ssp(N, K, Y)
  # index.plt <- random.index(N, n.plt, p=p.plt)
  # x.plt <- X[index.plt,]
  # y.plt <- Y[index.plt]
  # results <- softmax.coef.estimate(x.plt, y.plt, weights=p.plt[index.plt])
  # beta.plt <- results$beta
  # P1.plt <-  results$P1 # n.plt*(K+1) matrix
  # P.plt <- P1.plt[, -1] # n.plt*K matrix
  # ddL.plt <- softmax.ddL(X = x.plt, P = P1.plt[, -1], p = p.plt[index.plt], K, d, scale = N*n.plt)
  # dL.sq.plt <- softmax.dL.sq(X = x.plt, Y.matrix = Y.matrix[index.plt, ],
  #                            P = P1.plt[, -1], p = p.plt[index.plt], K = K, d = d, scale = (N^2*n.plt^2))

  # cov.plt <- solve(ddL.plt) %*% dL.sq.plt %*% solve(ddL.plt)
  # if (criterion == "MSPE") {
  #   # Omega.plt <- Omega(x.plt, K, d, n.plt, P = P.plt, p = p.plt)
  #   # Omega.test <- Omega.test(x.plt, K, d, n.plt, P1 = P1.plt, p = p.plt)
  #   # Omega.test1 <- Omega.test1(x.plt, P1 = P1.plt, p = p.plt, K, d, scale = N*n.plt)
  #
  #   Omega.plt <- Omega(x.plt, P1 = P1.plt, p = p.plt, K, d, scale = N*n.plt)
  #   # Omega.plt.test <- Omega.test(x.plt, P1 = P1.plt, p = p.plt, K, d, scale = N*n.plt)
  #   # print(Omega.plt == Omega.plt.test)
  #   # mbm <- microbenchmark::microbenchmark(
  #   #   Omega(x.plt, P1 = P1.plt, p = p.plt, K, d, scale = N*n.plt),
  #   #   Omega.test(x.plt, P1 = P1.plt, p = p.plt, K, d, scale = N*n.plt),
  #   #   unit = 'us'
  #   # )
  #   # print(mbm)
  # } else {Omega.plt <- NA}
  # P1.plt.N <- pbeta.multi(X, beta.plt) # N*(K+1)
  return(
    list(
      p.plt = p.plt,
      beta.plt = beta.plt,
      P1.plt = P1.plt.N, # N*(K+1)
      ddL.plt = ddL.plt,
      dL.sq.plt = dL.sq.plt,
      index.plt = index.plt,
      cov.plt = cov.plt,
      Omega.plt = Omega.plt,
      Lambda.plt = Lambda.plt
    )
  )
}
###############################################################################
# softmax.calculate.nm <- function(X, Y, ddL.plt, Omega.plt, sixi, G, criterion, constraint){
#   if (criterion == "OptA") {
#     if (constraint == "baseline"){
#       nm <- sqrt(colSums((solve(ddL.plt, t(sixi))^2)))
#     } else if (constraint == "summation"){
#       nm <- sqrt(colSums(((G %*% solve(ddL.plt, t(sixi)))^2)))
#     }
#   } else if (criterion == "OptL"){
#     if (constraint == "baseline"){
#       nm <- sqrt(rowSums((sixi)^2))
#     } else if (constraint == "summation"){
#       # print(G %*% solve(t(G) %*% G))
#       nm <- sqrt(colSums(((G %*% solve(t(G) %*% G) %*% t(sixi))^2)))
#     }
#   } else if (criterion == "MSPE") {
#     # Msixi <- solve(MN, t(sixi))
#     # pi.num <- sqrt(colSums((Omega.plt %*% Msixi) * Msixi))
#     nm <- sqrt(colSums(((expm::sqrtm(Omega.plt) %*% solve(ddL.plt, t(sixi)))^2)))
#   }
#   return(nm)
# }
###############################################################################
softmax.calculate.nm <- function(X, Y, ddL.plt, Omega.plt, sixi, G, criterion, constraint){

  if (criterion == "OptA") {
    if (constraint == "baseline"){
      nm <- sqrt(rowSums((sixi %*% solve(ddL.plt))^2))
    } else if (constraint == "summation"){
      temp <- sixi %*% solve(ddL.plt) %*% t(G)
      nm <- sqrt(rowSums(temp^2))
    }
  } else if (criterion == "OptL"){
    if (constraint == "baseline"){
      nm <- sqrt(rowSums((sixi)^2))
    } else if (constraint == "summation"){
      tempG <- t(G %*% solve(t(G) %*% G))
      temp <- sixi %*% tempG
      nm <- sqrt(rowSums(temp^2))
    }
  } else if (criterion == "MSPE") {
    temp <- sixi %*% solve(ddL.plt) %*% expm::sqrtm(Omega.plt)
    nm <- sqrt(rowSums(temp^2))
  }

  return(nm)
}
softmax.calculate.offset <- function(P1, X, G, ddL.plt, criterion, constraint) {

  if (criterion == "OptA") {
  } else if (criterion == "OptL"){
    offset <- log(rowSums(P1^2) + 1 - 2*P1) / 2
    # if (constraint == 'baseline') {
    #   offset <- log(rowSums(P1^2) + 1 - 2*P1) / 2
    # } else if (constraint == 'summation') {
    #   offset <- log(rowSums(P1^2) + 1 - 2*P1) / 2
    #   dim <- dim(P1)
    #   n <- dim[1] # subsample size
    #   K <- dim[2]
    #   d <- ncol(X)
    #   offset <- matrix(NA, nrow = n, ncol = K)
    #   for (k in 1:K) {
    #     P.temp <- -P1
    #     P.temp[, k] <- 1 - P.temp[, k]
    #     sixi <- P.temp[, rep(seq(K), each = d)] * X[, rep(seq(d), K)]
    #     offset[, k] <- log(sqrt(rowSums((sixi %*% solve(ddL.plt))^2)))
    #   }

  } else if (criterion == "MSPE") {

  }
  return(offset)
}
###############################################################################
softmax.subsampling <- function(X, Y, Y.matrix, G, n.ssp, N, K, d, alpha, b,
                                criterion, estimate.method, sampling.method, constraint,
                                p.plt, ddL.plt, P1.plt, Omega.plt, index.plt) {
  if (estimate.method == 'MSCLE'){
    E <- Y.matrix - P1.plt[, -1]
    E <- cbind(-rowSums(E), E)
    nm <- sqrt(rowSums((E)^2) * rowSums(X^2))

    # nm <- sqrt(rowSums((Y.matrix.1 - P1.plt)^2) * rowSums(X^2))
  } else if (estimate.method == 'Weighted'){
    sixi <- (Y.matrix-P1.plt[, -1])[, rep(seq(K), each = d)] * X[,rep(seq(d), K)] # N*(Kd))
    nm <- softmax.calculate.nm(X, Y, ddL.plt, Omega.plt, sixi, G, criterion, constraint)
  }

  # nm <- softmax.calculate.nm(X, Y, ddL.plt, Omega.plt, sixi, G, criterion, constraint)
  if (sampling.method == "WithReplacement"){
    dm <- sum(nm) # denominator
    p.ssp <- (1 - alpha) * nm / dm + alpha / N # N*1
    index.ssp <- random.index(N, n.ssp, p.ssp)
  } else if (sampling.method == "Poisson"){
    H <- quantile(nm[index.plt], 1-n.ssp/(b*N))
    nm[nm > H] <- H
    dm <- (sum(nm[index.plt] /p.plt) / n.plt) * (n.plt/(n.plt - K*d))
    # dm <- (sum(nm[index.plt] /p.plt)) * (n.plt/(n.plt - K*d)) # poi pilot
    p.ssp <- n.ssp * ((1 - alpha) * nm / dm + alpha / N)

    index.ssp <- poisson.index(N, p.ssp)
    # calculate offset
    if (estimate.method == 'MSCLE') {
      offset <- softmax.calculate.offset(P1.plt[index.ssp, ], X[index.ssp, ], G,
                                         ddL.plt, criterion, constraint)
    } else if (estimate.method == 'Weighted') {
      offset <- NA
    }
  }

  return(list(index.ssp = index.ssp, p.ssp = p.ssp, offset = offset))
}
###############################################################################
softmax.subsample.estimate <- function(x.ssp, y.ssp, y.matrix.ssp, n.ssp, index.ssp,
                                       p.ssp, offset,
                                       beta.plt, sampling.method, estimate.method,
                                       N, K, d) {
  if (estimate.method == "Weighted"){
    results <- softmax.coef.estimate(x.ssp, y.ssp, weights = 1 / p.ssp)
    beta.ssp <- results$beta
    P.ssp <- (results$P1)[, -1] # n.ssp*K matrix
    if (sampling.method == "WithReplacement") {
      ddL.ssp <- softmax.ddL(X = x.ssp, P = P.ssp, p = p.ssp,
                             K = K, d = d, scale = N*n.ssp)
      dL.sq.ssp <- softmax.dL.sq(X = x.ssp, Y.matrix = y.matrix.ssp,
                                 P = P.ssp, p = p.ssp,
                                 K = K, d = d, scale = (N^2)*(n.ssp^2))
      c <- n.ssp / N
      Lambda.ssp <- c * softmax.dL.sq(X = x.ssp, Y.matrix = y.matrix.ssp,
                                      P = P.ssp, p = sqrt(p.ssp),
                                      K = K, d = d, scale = N*(n.ssp^2))
    } else if (sampling.method == 'Poisson') {
      ddL.ssp <- softmax.ddL(X = x.ssp, P = P.ssp, p = p.ssp,
                             K = K, d = d, scale = N)
      dL.sq.ssp <- softmax.dL.sq(X = x.ssp, Y.matrix = y.matrix.ssp,
                                 P = P.ssp, p = p.ssp,
                                 K = K, d = d, scale = N^2)
      Lambda.ssp <- 0
    }
    cov.ssp.full <- solve(ddL.ssp) %*% (dL.sq.ssp ) %*% solve(ddL.ssp)
    cov.ssp <- solve(ddL.ssp) %*% (dL.sq.ssp + Lambda.ssp) %*% solve(ddL.ssp)
  } else if (estimate.method == 'MSCLE'){
    results <- softmax.coef.estimate(x.ssp, y.ssp, offset = offset)
    beta.ssp <- results$beta
    P.ssp <- (results$P1)[, -1] # n.ssp*K matrix
    if (sampling.method == "WithReplacement") {
      stop("The 'MSCLE' estimate method with 'WithReplacement' sampling method
           has not been implemented yet.")
    } else if (sampling.method == 'Poisson') {
      ddL.ssp <- softmax.ddL(X = x.ssp, P = P.ssp, p = p.ssp,
                             K = K, d = d, scale = N)
      dL.sq.ssp <- softmax.dL.sq(X = x.ssp, Y.matrix = y.matrix.ssp,
                                 P = P.ssp, p = p.ssp,
                                 K = K, d = d, scale = N^2)
      Lambda.ssp <- 0
    }
    cov.ssp.full <- solve(ddL.ssp) %*% (dL.sq.ssp ) %*% solve(ddL.ssp)
    cov.ssp <- solve(ddL.ssp) %*% (dL.sq.ssp + Lambda.ssp) %*% solve(ddL.ssp)
  }

  return(list(beta.ssp = beta.ssp,
              P.ssp = P.ssp,
              ddL.ssp = ddL.ssp,
              dL.sq.ssp = dL.sq.ssp,
              Lambda.ssp = Lambda.ssp,
              cov.ssp = cov.ssp,
              cov.ssp.full = cov.ssp.full)
  )
}
###############################################################################
softmax.combine.estimator.2 <- function(X, Y, Y.matrix, index.plt, index.ssp,
                                        P.plt, P.ssp, p.plt, p.ssp, N, K, d,
                                   sampling.method = sampling.method,
                                   estimate.method = estimate.method) {
  x.cmb <- rbind(X[index.plt,], X[index.ssp,])
  y.cmb <- c(Y[index.plt], Y[index.ssp])
  y.matrix.cmb <-rbind(Y.matrix[index.plt, ], Y.matrix[index.ssp, ])

  P.cmb <- rbind(P.plt[index.plt, -1] , P.ssp)

  n.plt <- length(p.plt)
  n.ssp <- length(p.ssp)
  # p.plt <- p.plt * n.plt
  p.cmb <- c(p.plt, p.ssp)
  p.cmb[1:n.plt] <- p.cmb[1:n.plt] / n.plt
  p.cmb[-(1:n.ssp)] <- p.cmb[-(1:n.ssp)] / n.ssp

  n.cmb <- length(p.cmb)

  ddL.cmb <- softmax.ddL(X = x.cmb, P = P.cmb, p = p.cmb,
                         K = K, d = d, scale = N)
  ddL.inv <- solve(ddL.cmb)
  print(ddL.inv)

  if (estimate.method == "Weighted"){
    results <- softmax.coef.estimate(x.cmb, y.cmb, weights = 1 / p.cmb)
    beta.cmb <- results$beta
    P.cmb <- (results$P1)[, -1] # n.ssp*K matrix
    if (sampling.method == "WithReplacement") {
      ddL.cmb <- softmax.ddL(X = x.cmb, P = P.cmb, p = p.cmb,
                             K = K, d = d, scale = N*n.cmb)
      dL.sq.cmb <- softmax.dL.sq(X = x.cmb, Y.matrix = y.matrix.cmb,
                                 P = P.cmb, p = p.cmb,
                                 K = K, d = d, scale = (N^2)*(n.cmb^2))
      c <- n.cmb / N
      Lambda.cmb <- c * softmax.dL.sq(X = x.cmb, Y.matrix = y.matrix.cmb,
                                      P = P.cmb, p = sqrt(p.cmb),
                                      K = K, d = d, scale = N*(n.cmb^2))
    } else if (sampling.method == 'Poisson') {

      ddL.cmb <- softmax.ddL(X = x.cmb, P = P.cmb, p = p.cmb,
                             K = K, d = d, scale = N)
      dL.sq.cmb <- softmax.dL.sq(X = x.cmb, Y.matrix = y.matrix.cmb,
                                 P = P.cmb, p = p.cmb,
                                 K = K, d = d, scale = N^2)
      Lambda.cmb <- 0
    }
    cov.cmb.full <- solve(ddL.cmb) %*% (dL.sq.cmb ) %*% solve(ddL.cmb)
    cov.cmb <- solve(ddL.cmb) %*% (dL.sq.cmb + Lambda.cmb) %*% solve(ddL.cmb)
  }
  P.cmb <- matrix(NA, nrow = N, ncol = K+1)
  P.cmb <- pbeta.multi(X, beta.cmb)
  return(list(beta.cmb = beta.cmb,
              cov.cmb = cov.cmb,
              cov.cmb.full = cov.cmb.full,
              P.cmb = P.cmb)
  )
}
###############################################################################
softmax.combine.sample <- function(X, Y, Y.matrix, index.plt, index.ssp, p.plt, p.ssp, N, K, d,
                                   sampling.method = sampling.method,
                                   estimate.method = estimate.method) {
  x.cmb <- rbind(X[index.plt,], X[index.ssp,])
  y.cmb <- c(Y[index.plt], Y[index.ssp])
  y.matrix.cmb <-rbind(Y.matrix[index.plt, ], Y.matrix[index.ssp, ])
  n.plt <- length(p.plt)
  # p.plt <- p.plt * n.plt
  p.cmb <- c(p.plt, p.ssp)
  n.cmb <- length(p.cmb)
  if (estimate.method == "Weighted"){
    results <- softmax.coef.estimate(x.cmb, y.cmb, weights = 1 / p.cmb)
    beta.cmb <- results$beta
    P.cmb <- (results$P1)[, -1] # n.ssp*K matrix
    if (sampling.method == "WithReplacement") {
      ddL.cmb <- softmax.ddL(X = x.cmb, P = P.cmb, p = p.cmb,
                             K = K, d = d, scale = N*n.cmb)
      dL.sq.cmb <- softmax.dL.sq(X = x.cmb, Y.matrix = y.matrix.cmb,
                                 P = P.cmb, p = p.cmb,
                                 K = K, d = d, scale = (N^2)*(n.cmb^2))
      c <- n.cmb / N
      Lambda.cmb <- c * softmax.dL.sq(X = x.cmb, Y.matrix = y.matrix.cmb,
                                      P = P.cmb, p = sqrt(p.cmb),
                                      K = K, d = d, scale = N*(n.cmb^2))
    } else if (sampling.method == 'Poisson') {

      ddL.cmb <- softmax.ddL(X = x.cmb, P = P.cmb, p = p.cmb,
                             K = K, d = d, scale = N)
      dL.sq.cmb <- softmax.dL.sq(X = x.cmb, Y.matrix = y.matrix.cmb,
                                 P = P.cmb, p = p.cmb,
                                 K = K, d = d, scale = N^2)
      Lambda.cmb <- 0
    }
    cov.cmb.full <- solve(ddL.cmb) %*% (dL.sq.cmb ) %*% solve(ddL.cmb)
    cov.cmb <- solve(ddL.cmb) %*% (dL.sq.cmb + Lambda.cmb) %*% solve(ddL.cmb)
  }
  P.cmb <- matrix(NA, nrow = N, ncol = K+1)
  P.cmb <- pbeta.multi(X, beta.cmb)
  return(list(beta.cmb = beta.cmb,
              cov.cmb = cov.cmb,
              cov.cmb.full = cov.cmb.full,
              P.cmb = P.cmb)
  )
}
###############################################################################

softmax.combining <- function(ddL.plt, ddL.ssp, dL.sq.plt, dL.sq.ssp,
                              Lambda.plt, Lambda.ssp,
                              n.plt, n.ssp, beta.plt, beta.ssp, X, N, K, d) {
  ddL.plt <- n.plt * ddL.plt
  ddL.ssp <- n.ssp * ddL.ssp
  dL.sq.plt <- n.plt^2 * dL.sq.plt
  dL.sq.ssp <- n.ssp^2 * dL.sq.ssp
  Lambda.plt <- n.plt^2 * Lambda.plt
  Lambda.ssp <- n.ssp^2 * Lambda.ssp

  #

  # # print(beta.plt)
  # # print(c(beta.plt))
  # print(ddL.plt)
  # print(ddL.ssp)
  # print(dL.sq.plt)
  # print(dL.sq.ssp)
  ddL.inv <- solve(ddL.plt + ddL.ssp)
  print(ddL.inv)
  beta.cmb <- ddL.inv %*% (ddL.plt %*% c(beta.plt) + ddL.ssp %*% c(beta.ssp))
  beta.1 <- ddL.inv %*% (ddL.plt %*% c(beta.plt))
  beta.2 <- ddL.inv %*% (ddL.ssp %*% c(beta.ssp))
  # print(beta.1)
  # print(beta.2)


  beta.cmb <- matrix(beta.cmb, nrow = d)
  cov.cmb.full <- ddL.inv %*% (dL.sq.plt + dL.sq.ssp + Lambda.ssp) %*% ddL.inv
  cov.cmb <- ddL.inv %*% (dL.sq.plt + Lambda.plt + dL.sq.ssp + Lambda.ssp) %*% ddL.inv
  # print(cov.cmb - cov.cmb.old)
  # cov.cmb <- cov.cmb.full + ddL.inv %*% (Lambda.ssp) %*% ddL.inv
  # var.1 <- diag(ddL.inv %*% dL.sq.plt %*% ddL.inv)
  # var.2 <- diag(ddL.inv %*% Lambda.plt %*% ddL.inv)
  # var.3 <- diag(ddL.inv %*% dL.sq.ssp %*% ddL.inv)
  # var.4 <- diag(ddL.inv %*% Lambda.ssp %*% ddL.inv)
  # print(mean(Lambda.plt/dL.sq.plt )) # should be n.plt / N
  # print(mean(dL.sq.plt/ dL.sq.ssp)) # should be n.plt / n.ssp
  # print(mean(Lambda.ssp/dL.sq.ssp )) # should be n.ssp / N
  # print(mean(Lambda.plt/ Lambda.ssp))# should be (n.plt / n.ssp)^2
  P.cmb <- matrix(NA, nrow = N, ncol = K+1)
  P.cmb <- pbeta.multi(X, beta.cmb)
  # P.cmb[, 1] <- 1 - rowSums(P.cmb[, -1])

  return(list(beta.cmb = beta.cmb,
              cov.cmb = cov.cmb,
              cov.cmb.full = cov.cmb.full,
              P.cmb = P.cmb)
  )
}
###############################################################################
fitting.measure <-  function(beta.pred = NA, beta.true = NA, P.pred = NA,
                             P.true = NA) {
  MSE <- sum((beta.pred - beta.true)^2)
  MSPE <- sum((P.pred - P.true)^2)
  return(list(MSE = MSE,
              MSPE = MSPE))
}
###############################################################################
subsampling.summary <- function(object) {
  coef <- object$beta.cmb
  se <- sqrt(diag(object$cov.cmb))
  N <- object$N
  n.ssp.expect <- object$subsample.size.expect
  n.ssp.actual <- length(object$index)
  n.ssp.unique <- length(unique(object$index))
  subsample.rate.expect <- (n.ssp.expect / N) * 100
  subsample.rate.actual <- (n.ssp.actual / N) * 100
  subsample.rate.unique <- (n.ssp.unique / N) * 100
  cat("Model Summary\n\n")
  cat("\nCall:\n")
  cat("\n")
  print(object$model.call)
  cat("\n")
  cat("Subsample Size:\n")
  size_table <- data.frame(
    'Variable' = c(
      'Total Sample Size',
      'Expected Subsample Size',
      'Actual Subsample Size',
      'Unique Subsample Size',
      'Expected Subample Rate',
      'Actual Subample Rate',
      'Unique Subample Rate'
    ),
    'Value' = c(
      N,
      n.ssp.expect,
      n.ssp.actual,
      n.ssp.unique,
      paste0(subsample.rate.expect, "%"),
      paste0(subsample.rate.actual, "%"),
      paste0(subsample.rate.unique, "%")
    )
  )
  colnames(size_table) <- NULL
  rownames(size_table) <- NULL
  print(size_table)
  cat("\n")
  cat("Coefficients:\n")
  cat("\n")
  coef_table <- data.frame(
    Estimate = round(coef, digits = 4),
    `Std. Error` = round(se, digits = 4),
    `z value` = round(coef / se, digits = 4),
    `Pr(>|z|)` = format.p.values(2 * (1 - pnorm(abs(coef / se))),
                                 threshold = 0.0001),
    check.names = FALSE
  )
  rownames(coef_table) <- names(coef)
  print(coef_table)
  # Add more summary information as needed
}
