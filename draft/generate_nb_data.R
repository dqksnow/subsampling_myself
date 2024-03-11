generate.nb.data <- function(N, case, beta0, v = 2, seed) { # as same as poisson
  set.seed(seed)
  d <- length(beta0) - 1

  if (case == 1) {
    X <- matrix(runif(N * d), N, d)
  } else if (case == 2) {
    X <- matrix(runif(N * d), N, d)
    epsilon <- runif(N, 0, 0.1)
    X[, 2] <- X[, 1] + epsilon
  } else if (case == 3) {
    X <- matrix(runif(N * d), N, d)
    epsilon <- runif(N)
    X[, 2] <- X[, 1] + epsilon
  } else if (case == 4) {
    X <- matrix(runif(N * d, min = -1, max = 1), N, d)
    epsilon <- runif(N)
    X[, 2] <- X[, 1] + epsilon
  }

  mu <- exp(beta0[1] + X %*% beta0[-1])
  Y <- rnbinom(N, size = v, mu = mu)

  return(list(
    X = X,
    Y = Y
  ))
}
