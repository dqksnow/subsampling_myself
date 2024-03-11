softmax.optimal.subsampling <- function(X, Y, n.plt, n.ssp,
                                        criterion = c('OptL', 'OptA', 'MSPE'),
                                        sampling.method = c('Poisson', 'WithReplacement'),
                                        estimate.method = c('Weighted', 'MSCLE', 'Uniform'),
                                        constraint = c('baseline', 'summation'),
                                        alpha = 0,
                                        b = 2) {
  model.call <- match.call()
  # family <- match.arg(family)

  # mf <- model.frame(formula, data)
  # Extract the response and predictor variables from the model frame
  # Y <- model.response(mf, "any")
  # if (is.character(Y) && length(unique(Y)) == 2) {
  #   levels <- unique(Y)
  #   Y <- as.integer(Y == levels[2])  # Assuming levels[2] is the 'success' category
  # }
  # X <- model.matrix(formula, mf)
  # colnames(X)[1] <- "intercept"
  # print(X)
  N <- nrow(X)
  d <- ncol(X)
  K <- length(unique(Y)) - 1
  G <- rbind(rep(-1/(K+1), K), diag(K) - 1/(K+1)) %x% diag(d)
  Y.matrix <- matrix(0, nrow = N, ncol = K)
  Y.matrix[cbind(seq_along(Y), Y)] <- 1

  if (estimate.method %in% c("Weighted", "MSCLE")) {
    # pilot step
    plt.estimate.results <- softmax.plt.estimate(X = X, Y = Y, Y.matrix,
                                                 n.plt = n.plt, N, K, d,
                                                 criterion)
    p.plt <- plt.estimate.results$p.plt
    beta.plt <- plt.estimate.results$beta.plt
    P1.plt = plt.estimate.results$P1.plt
    ddL.plt <- plt.estimate.results$ddL.plt
    dL.sq.plt <- plt.estimate.results$dL.sq.plt
    index.plt <- plt.estimate.results$index.plt
    cov.plt.b <- plt.estimate.results$cov.plt
    Omega.plt <- plt.estimate.results$Omega.plt
    Lambda.plt <- plt.estimate.results$Lambda.plt
    # n.plt <- length(index.plt)
    # subsampling step
    ssp.results <- softmax.subsampling(X = X,
                                       Y = Y,
                                       Y.matrix = Y.matrix,
                                       G = G,
                                       n.ssp = n.ssp,
                                       N = N, K = K, d = d,
                                       alpha = alpha,
                                       b = b,
                                       criterion = criterion,
                                       estimate.method = estimate.method,
                                       sampling.method = sampling.method,
                                       constraint = constraint,
                                       p.plt = p.plt,
                                       ddL.plt = ddL.plt,
                                       P1.plt = P1.plt,
                                       Omega.plt = Omega.plt,
                                       index.plt = index.plt)
    index.ssp <- ssp.results$index.ssp
    p.ssp <- ssp.results$p.ssp
    offset <- ssp.results$offset
    n.ssp <- length(index.ssp)
    # subsample estimating step
    ssp.estimate.results <- softmax.subsample.estimate(X[index.ssp, ],
                                                       Y[index.ssp],
                                                       Y.matrix[index.ssp, ],
                                                       n.ssp = n.ssp,
                                                       index.ssp = index.ssp,
                                                       p.ssp = p.ssp[index.ssp],
                                                       offset = offset,
                                                       beta.plt = beta.plt,
                                                       sampling.method = sampling.method,
                                                       estimate.method = estimate.method,
                                                       N, K, d)
    beta.ssp <- ssp.estimate.results$beta.ssp
    P.ssp <- ssp.estimate.results$beta.ssp
    ddL.ssp <- ssp.estimate.results$ddL.ssp
    dL.sq.ssp <- ssp.estimate.results$dL.sq.ssp
    Lambda.ssp <- ssp.estimate.results$Lambda.ssp
    cov.ssp.b <- ssp.estimate.results$cov.ssp
    cov.ssp.full.b <- ssp.estimate.results$cov.ssp.full
    # combining step
    # combining.results <- softmax.combining(ddL.plt = ddL.plt,
    #                                        ddL.ssp = ddL.ssp,
    #                                        dL.sq.plt = dL.sq.plt,
    #                                        dL.sq.ssp = dL.sq.ssp,
    #                                        Lambda.plt = Lambda.plt,
    #                                        Lambda.ssp = Lambda.ssp,
    #                                        n.plt = n.plt,
    #                                        n.ssp = n.ssp,
    #                                        beta.plt = beta.plt,
    #                                        beta.ssp = beta.ssp,
    #                                        X = X, N = N, K = K, d = d)
    # beta.cmb <- combining.results$beta.cmb
    # cov.cmb.b <- combining.results$cov.cmb
    # cov.cmb.full.b <- combining.results$cov.cmb.full
    # P.cmb <- combining.results$P.cmb
    # combining step
    combining.results <- softmax.combine.estimator.2(X, Y, Y.matrix,
                                                index.plt, index.ssp,
                                                P1.plt, P.ssp,
                                                p.plt, p.ssp[index.ssp], N, K, d,
                                                sampling.method = sampling.method,
                                                estimate.method = estimate.method)

    beta.cmb <- combining.results$beta.cmb
    cov.cmb.b <- combining.results$cov.cmb
    cov.cmb.full.b <- combining.results$cov.cmb.full
    P.cmb <- combining.results$P.cmb
    # # combining step
    # combining.results <- softmax.combine.sample(X, Y, Y.matrix,
    #                                             index.plt, index.ssp, p.plt, p.ssp[index.ssp], N, K, d,
    #                                             sampling.method = sampling.method,
    #                                             estimate.method = estimate.method)
    #
    # beta.cmb <- combining.results$beta.cmb
    # cov.cmb.b <- combining.results$cov.cmb
    # cov.cmb.full.b <- combining.results$cov.cmb.full
    # P.cmb <- combining.results$P.cmb
    # results

    beta.plt <- G %*% as.vector(beta.plt)
    beta.ssp <- G %*% as.vector(beta.ssp)
    beta.cmb <- G %*% as.vector(beta.cmb)
    # beta.plt <- as.vector(beta.plt)
    # beta.ssp <- as.vector(beta.ssp)
    # beta.cmb <- as.vector(beta.cmb)
    cov.plt.s <- G %*% cov.plt.b %*% t(G)
    cov.ssp.s <- G %*% cov.ssp.b %*% t(G)
    cov.ssp.full.s <- G %*% cov.ssp.full.b %*% t(G)
    cov.cmb.s <- G %*% cov.cmb.b %*% t(G)
    cov.cmb.full.s <- G %*% cov.cmb.full.b %*% t(G)
    # beta.plt <- matrix(G %*% as.vector(t(beta.plt)), nrow = d)
    # beta.ssp <- matrix(G %*% as.vector(t(beta.ssp)), nrow = d)
    # beta.cmb <- matrix(G %*% as.vector(t(beta.cmb)), nrow = d)

    return(list(
      model.call = model.call,
      beta.plt = beta.plt,
      beta.ssp = beta.ssp,
      beta.cmb = beta.cmb,
      cov.plt = cov.plt.s,
      cov.ssp = cov.ssp.s,
      cov.ssp.full = cov.ssp.full.s,
      cov.cmb = cov.cmb.s,
      cov.cmb.full = cov.cmb.full.s,
      P.cmb = P.cmb,
      index.plt = index.plt,
      index.ssp = index.ssp,
      subsample.size.expect = n.ssp))
  } else if (estimate.method == "Uniform"){
    if (sampling.method == 'WithReplacement') {
      n.uni <- n.plt + n.ssp
      # uni sampling with replacement
      # index.uni <- poisson.index(N, n.uni/N)
      index.uni <- random.index(N, n.uni)
      # index.uni <- sample(N, n.uni, replace = F)
      x.uni <- X[index.uni, ]
      y.uni <- Y[index.uni]
      results <- softmax.coef.estimate(x.uni, y.uni)
      beta.uni <- results$beta
      P.uni <- results$P1
      # cov.uni.nnet <- results$var
      ddL.uni <- softmax.ddL(X = x.uni, P = P.uni[, -1], p = 1, K, d,
                             scale = N*n.uni)
      dL.sq.uni <- softmax.dL.sq(X = x.uni, Y.matrix = Y.matrix[index.uni, ],
                                 P = P.uni[, -1], p = 1, K = K, d = d,
                                 scale=(N^2)*(n.uni^2))
      c <- n.uni / N

      cov.uni.b <- solve(ddL.uni) %*% (dL.sq.uni * (1+c))%*% solve(ddL.uni)
      cov.uni.nnet <- G %*% solve(ddL.uni) %*% (dL.sq.uni) %*% solve(ddL.uni)%*% t(G)
      cov.uni.s <- G %*% cov.uni.b %*% t(G)

      beta <- G %*% as.vector(beta.uni)
    # beta <- as.vector(beta.uni)
    # P.uni <- matrix(NA, nrow = N, ncol = K+1)
    # P.uni <- pbeta.multi(X, beta.uni) # N*(K+1)
#############################################################################
    # uni poisson sampling
    } else if (sampling.method == 'Poisson') {
      n.uni <- n.plt + n.ssp
      index.uni <- poisson.index(N, n.uni/N)
      p.uni <- rep(n.uni / N, length(index.uni))
      x.uni <- X[index.uni, ]
      y.uni <- Y[index.uni]
      results <- softmax.coef.estimate(x.uni, y.uni)
      beta.uni <- results$beta
      P.uni <- results$P1
      ddL.uni <- softmax.ddL(X = x.uni, P = P.uni[, -1], p = p.uni, K, d,
                             scale = N)
      dL.sq.uni <- softmax.dL.sq(X = x.uni, Y.matrix = Y.matrix[index.uni, ],
                                 P = P.uni[, -1], p = p.uni, K = K, d = d,
                                 scale=(N^2))
      cov.uni.b <- solve(ddL.uni) %*% (dL.sq.uni) %*% solve(ddL.uni)
      cov.uni.s <- G %*% cov.uni.b %*% t(G)
      beta <- G %*% as.vector(beta.uni)
    }
    return(list(
      model.call = model.call,
      index = index.uni,
      beta = beta,
      cov = cov.uni.s,
      P = P.uni,
      cov.nnet = cov.uni.nnet))
  }
}
