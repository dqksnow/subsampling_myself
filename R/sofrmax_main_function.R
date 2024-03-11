softmax.optimal.subsampling <- function(X, Y, n.plt, n.ssp,
                                    criterion = c('OptL', 'OptA'),
                                    sampling.method = c('Poisson', 'WithReplacement'),
                                    estimate.method = c('Weighted', 'Uniform'),
                                    constraint = c('baseline', 'summation'),
                                    alpha = 0,
                                    b = 2) {
  # model.call <- match.call()
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
  print(nrow(X))
  N <- nrow(X)
  d <- ncol(X)
  K <- unique(Y) - 1
  Y.matrix <- matrix(0, nrow = N, ncol = K)
  Y.matrix[cbind(seq_along(Y), Y)] <- 1

  if (estimate.method == "Weighted") {
    # pilot step
    plt.estimate.results <- pilot.estimate(X = X, Y = Y, Y.matrix,
                                           n.plt = n.plt, N, K, d)
    p.plt <- plt.estimate.results$p.plt
    beta.plt <- plt.estimate.results$beta.plt
    P.plt <- plt.estimate.results$P.plt
    ddL.plt <- plt.estimate.results$ddL.plt
    dL.sq.plt <- plt.estimate.results$dL.sq.plt
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
                               P.plt = P.plt,
                               index.plt = index.plt)
    index.ssp <- ssp.results$index.ssp
    w.ssp <- ssp.results$w.ssp

    # subsample estimating step
    ssp.estimate.results <- subsample.estimate(X[index.ssp, ],
                                               Y[index.ssp],
                                               Y.matrix[index.ssp, ],
                                               n.ssp = n.ssp,
                                               w.ssp = w.ssp,
                                               beta.plt = beta.plt,
                                               sampling.method = sampling.method,
                                               estimate.method = estimate.method)
    beta.ssp <- ssp.estimate.results$beta.ssp
    ddL.ssp <- ssp.estimate.results$ddL.ssp
    dL.sq.ssp <- ssp.estimate.results$dL.sq.ssp
    var.ssp <- ssp.estimate.results$var.ssp

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
                # var.plt = var.plt,
                var.ssp = var.ssp,
                var.ssp.true = var.ssp.true,
                var.cmb = var.cmb,
                var.cmb.true = var.cmb.true,
                index.plt = index.plt,
                index.ssp = index.ssp))
  } else if (estimate.method == "Uni"){
    n.uni <- n.plt + n.ssp
    index.uni <- random.index(N, n.uni)
    x.uni <- X[index.uni, ]
    y.uni <- Y[index.uni]
    beta.uni <- softmax.coef.estimate(x.uni, y.uni)
    P.uni <- pbeta.multi(x.uni, beta.uni)
    # pbeta.uni <- pbeta(X[index.uni, ], beta.uni)
    # ddL.uni <- MN(X.uni, pbeta.uni, 1 / n.uni) # 2nd derivative of the log likelihood function.
    # Psi.uni <- Psi(X.uni, Y.uni, pbeta.uni, 1 / n.uni^2)
    ddL.uni <- ddL(X = x.uni, P = P.uni[index.uni, ], p = 1 / N, K, d, N, n.uni)
    dL.sq.uni <- dL.sq(X = x.uni, Y.matrix = Y.matrix[index.uni, ],
                       P = P.uni[index.uni, ], p = p.uni, K = K, d = d, N = N,
                       n = n.uni)
    var.uni <- solve(ddL.uni) %*% dL.sq.uni %*% solve(ddL.uni)

    # prediction.results <- in.sample.prediction(X, Y, beta.uni)

    return(list(model.call = model.call,
                index = index.uni,
                beta = beta.uni,
                var = var.uni))
  }
}
