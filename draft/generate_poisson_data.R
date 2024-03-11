generate.poisson.data <- function(N, case, beta0, seed) {
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

  lambda <- exp(beta0[1] + X %*% beta0[-1])
  Y <- rpois(N, lambda)

  return(list(
    X = X,
    Y = Y
  ))
}

link_function <- function(X, beta0) 1 / (beta0[1] + X %*% beta0[-1])

generate.gamma.data <- function(N, case, beta0, shape, seed) {
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

  scale <- link_function(X, beta0) / shape
  Y <- rgamma(N, shape = shape, scale = scale)
  return(list(
    X = X,
    Y = Y
  ))
}
