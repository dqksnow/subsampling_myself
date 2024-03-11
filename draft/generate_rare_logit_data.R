generate.rare.data <- function(N, case, beta0) {
  if (case == 1) {
    beta0[1] <- -7.65
  } else if (case == 2) {
    beta0[1] <- -0.5
  } else if (case == 3) {
    beta0[1] <- -7
  } else if (case == 4) {
    beta0[1] <- -1.8
  } else if (case == 5) {
    beta0[1] <- -3.2
  }
  d <- length(beta0)
  ds <- d - 1
  corr <- 0.5
  sigmax <- corr ^ abs(outer(1:ds, 1:ds, "-"))
  sigmax <- sigmax / 4

  if (case == 1) { # Normal
    Z <- MASS::mvrnorm(n = N, mu = rep(0, ds), Sigma = sigmax)
  } else if (case == 2) { # lognormal, inbalanced
    Z <- exp(MASS::mvrnorm(n = N, mu = rep(0, ds), Sigma = sigmax))
  } else if (case == 3) { # T3
    Z <- MASS::mvrnorm(n = N, mu = rep(0, ds), Sigma = sigmax)
    df <- 3
    Z <- Z / sqrt(rchisq(n = N ,df = df) / df) / 3
  } else if (case == 4) { # exponential
    Z <- matrix(rexp(N * ds), nrow = N)
  } else if (case == 5) { # uniform
    Z <- matrix(runif(N * ds), nrow = N)
  }

  P <- 1 - 1 / (1 + exp(beta0[1] + Z %*% beta0[-1])) # equal to  P <- plogis(beta0[1] + Z %*% beta0[-1])
  Y <- rbinom(n = N, size = 1, p = P)

  return(list(
    X = Z,
    Y = Y,
    beta0 = beta0
  ))
}
# generate.rare.data <- function(N, case, beta0) {
#   if (case == 1) {
#     beta0[1] <- -7.65
#   } else if (case == 2) {
#     beta0[1] <- -0.5
#   } else if (case == 3) {
#     beta0[1] <- -7
#   } else if (case == 4) {
#     beta0[1] <- -1.8
#   } else if (case == 5) {
#     beta0[1] <- -3.2
#   }
#   d <- length(beta0)
#   ds <- d - 1
#   corr <- 0.5
#   sigmax <- corr ^ abs(outer(1:ds, 1:ds, "-"))
#   sigmax <- sigmax / 4
#
#   if (case == 1) { # Normal
#     Z <- MASS::mvrnorm(n = N, mu = rep(0, ds), Sigma = sigmax)
#   } else if (case == 2) { # lognormal, inbalanced
#     Z <- exp(MASS::mvrnorm(n = N, mu = rep(0, ds), Sigma = sigmax))
#   } else if (case == 3) { # T3
#     Z <- MASS::mvrnorm(n = N, mu = rep(0, ds), Sigma = sigmax)
#     df <- 3
#     Z <- Z / sqrt(rchisq(n = N ,df = df) / df) / 3
#   } else if (case == 4) { # exponential
#     Z <- matrix(rexp(N * ds), nrow = N)
#   } else if (case == 5) { # uniform
#     Z <- matrix(runif(N * ds), nrow = N)
#   }
#
#   P <- 1 - 1 / (1 + exp(beta0[1] + Z %*% beta0[-1])) # equal to  P <- plogis(beta0[1] + Z %*% beta0[-1])
#   Y <- rbinom(n = N, size = 1, p = P)
#
#   return(list(
#     X = cbind(rep(1, N), Z),
#     Y = Y,
#     beta0 = beta0
#   ))
# }
