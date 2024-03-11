#### JASA2018 simu1 multinormal balanced ########
simu_mzNormal <- function(seed, N, beta0, corr){
  d <- length(beta0) - 1
  sigmax  <- matrix(corr, d, d) + diag(1-corr, d)
  set.seed(seed)
  X <- mvtnorm::rmvnorm(N, rep(0, d), sigmax)
  colnames(X) <- paste("X", 1:ncol(X), sep = "")
  P <- 1 - 1 / (1 + exp(beta0[1] + X %*% beta0[-1]))
  Y <- rbinom(N, 1, P)
  # X_numeric <- mvtnorm::rmvnorm(N, rep(0, d), sigmax)
  # predictor <- sample(c("alpha", "beta", 'gamma'), size = N, replace = TRUE)
  # predictor_df <- data.frame(predictor = as.factor(predictor))
  # X_to_model <- as.data.frame(cbind(X_numeric, predictor_df))
  # colnames(X_to_model) <- paste("X", 1:ncol(X_to_model), sep = "")
  # predictor_dummies <- model.matrix(~ predictor - 1, data = predictor_df)
  # X <- cbind(X_numeric, predictor_dummies)
  # beta0 <- rep(0.5, 10)
  # P <- 1 - 1 / (1 + exp(beta0[1] + X %*% beta0[-1]))
  # Y <- rbinom(N, 1, P)



  # P <- 1 - 1 / (1 + exp(beta0[1] + X %*% beta0[-1]))
  # Y <- rbinom(N, 1, P)
  # X[, 1] <- predictor
  return(list(X=X, Y=Y))
}

# #### JASA2018 simu1 multinormal balanced ########
# simu_mzNormal <- function(seed, N, beta0, corr){
#   d <- length(beta0) - 1
#   sigmax  <- matrix(corr, d, d) + diag(1-corr, d)
#   set.seed(seed)
#   X <- mvtnorm::rmvnorm(N, rep(0, d), sigmax)
#   X <- cbind(1, X)
#   P <- 1 - 1 / (1 + exp(X %*% beta0))
#   Y <- rbinom(N, 1, P)
#
#   # beta_full <- logistic_coef_estimate(X, Y, 1, 1:length(Y))
#   return(list(X=X, Y=Y))
# }

#### JASA2018 simu2 multinormal imbalanced ########
simu_nzNormal <- function(seed, N, beta0, corr){
  d <- length(beta0) - 1
  sigmax  <- matrix(corr, d, d) + diag(1-corr, d)
  set.seed(seed)
  X <- mvtnorm::rmvnorm(N, rep(1.5, d), sigmax)
  P <- 1 - 1 / (1 + exp(beta0[1] + X %*% beta0[-1]))
  Y <- rbinom(N, 1, P)
  # beta_full <- logistic_coef_estimate(X, Y, 1, 1:length(Y))
  return(list(X=X, Y=Y))
}
# simu_nzNormal <- function(seed, N, beta0, corr){
#   d <- length(beta0) - 1
#   sigmax  <- matrix(corr, d, d) + diag(1-corr, d)
#   set.seed(seed)
#   X <- mvtnorm::rmvnorm(N, rep(1.5, d), sigmax)
#   X <- cbind(1, X)
#   P <- 1 - 1 / (1 + exp(X %*% beta0))
#   Y <- rbinom(N, 1, P)
#   # beta_full <- logistic_coef_estimate(X, Y, 1, 1:length(Y))
#   return(list(X=X, Y=Y))
# }

#### JASA2018 simu3 hetero multinormal ########
simu_ueNormal <- function(seed, N, beta0, corr){
  d <- length(beta0) - 1
  sigmax <- matrix(0,d,d)
  var <- 1/c(1:d)^2 #
  #var <- c(1:(d)) #^2
  for (i in 1:d){
    for (j in 1:d){
      sigmax[i,j] <- corr^(i!=j) * (sqrt(var[i])*sqrt(var[j]))
    }
    }
  #sigmax <- matrix(corr, d, d)  + diag(as.vector(var)) - diag(corr, d)
  X <- mvtnorm::rmvnorm(N, rep(0, d), sigmax)
  P <- 1 - 1 / (1 + exp(beta0[1] + X %*% beta0[-1]))
  Y <- rbinom(N, 1, P)
  # beta_full <- logistic_coef_estimate(X, Y, 1, 1:length(Y))
  return(list(X=X, Y=Y))
}

#### JASA2018 simu4 biomodal normal ########
simu_mixNormal <- function(seed, N, beta0, corr){
  d <- length(beta0) - 1
  sigmax  <- matrix(corr, d, d) + diag(1-corr, d)
  set.seed(seed)
  X <- 0.5*mvtnorm::rmvnorm(N, rep(1, d), sigmax) + 0.5*rmvnorm(N, rep(-1, d), sigmax)
  P <- 1 - 1 / (1 + exp(beta0[1] + X %*% beta0[-1]))
  Y <- rbinom(N, 1, P)
  # beta_full <- logistic_coef_estimate(X, Y, 1, 1:length(Y))
  return(list(X=X, Y=Y))
}

#### JASA2018 simu5 multi t ########
simu_T3 <- function(seed, N, beta0, corr, df){
  d <- length(beta0) - 1
  sigmax  <- matrix(corr, d, d) + diag(1-corr, d)
  set.seed(seed)
  X <-mvtnorm::rmvt(N, sigmax, df, rep(0, d))/10#, type=c("shifted","Kshirsagar")[1]
  P <- 1 - 1 / (1 + exp(beta0[1] + X %*% beta0[-1]))
  Y <- rbinom(N, 1, P)
  # beta_full <- logistic_coef_estimate(X, Y, 1, 1:length(Y))
  return(list(X=X, Y=Y))
}

#### JASA2018 simu6 exp ########
# simu_EXP <- function(seed, N, beta0, rate){
#   d <- length(beta0)
#   set.seed(seed)
#   X <- matrix(0, N, d)
#   generate_rexp <- function(x) x<-rexp(N, rate=rate)
#   X <- apply(X, 2, generate_rexp)
#   P <- 1 - 1 / (1 + exp(X %*% beta0))
#   Y <- rbinom(N, 1, P)
#   # beta_full <- logistic_coef_estimate(X, Y, 1, 1:length(Y))
#   return(list(X=X, Y=Y))
# }
simu_EXP <- function(seed, N, beta0, rate){
  d <- length(beta0) - 1
  set.seed(seed)
  X <- matrix(0, N, d)
  generate_rexp <- function(x) x<-rexp(N, rate=rate)
  X <- apply(X, 2, generate_rexp)
  P <- 1 - 1 / (1 + exp(beta0[1] + X %*% beta0[-1]))
  Y <- rbinom(N, 1, P)
  # beta_full <- logistic_coef_estimate(X, Y, 1, 1:length(Y))
  return(list(X=X, Y=Y))
}

#
predict_classification <- function(X, Y, beta){
  pbeta_estimate <- pbeta(X, beta)
  Y_pred <- rep(0, length(Y))
  Y_pred[which(pbeta_estimate>0.5)] <- 1
  proportion <- sum(Y_pred==Y)/length(Y)
  return(proportion)
}
###############################################################################
simu_softmax <- function(seed, case, N, d, K, beta){
  set.seed(seed)
  if (case == 1) {
    mu <- rep(0, d)
    # sigma <- matrix(0.5, nrow = d, ncol = d)
    sigma <- matrix(0.5, nrow = d, ncol = d)
    diag(sigma) <- rep(1, d)
    X <- MASS::mvrnorm(N, mu, sigma)
    prob.Kplus2 <- exp( X %*% beta)
    prob.Kplus2 <- prob.Kplus2 / rowSums(prob.Kplus2)
    Y <- apply(prob.Kplus2, 1, function(row) sample(0:K, size = 1, prob = row))
  } else if (case == 2) {
    mu <- rep(1.5, d)
    sigma <- matrix(0.5, nrow = d, ncol = d)
    diag(sigma) <- rep(1, d)
    X <- MASS::mvrnorm(N, mu, sigma)
    prob.Kplus2 <- exp( X %*% beta)
    prob.Kplus2 <- prob.Kplus2 / rowSums(prob.Kplus2)
    Y <- apply(prob.Kplus2, 1, function(row) sample(0:K, size = 1, prob = row))
  } else if (case == 3) {
    mu1 <- rep(1, d)
    mu2 <- rep(-1, d)
    sigma <- matrix(0.5, nrow = d, ncol = d)
    diag(sigma) <- rep(1, d)
    X <- 0.5 * MASS::mvrnorm(N, mu1, sigma) + 0.5 * MASS::mvrnorm(N, mu2, sigma)
    prob.Kplus2 <- exp( X %*% beta)
    prob.Kplus2 <- prob.Kplus2 / rowSums(prob.Kplus2)
    Y <- apply(prob.Kplus2, 1, function(row) sample(0:K, size = 1, prob = row))
  } else if (case == 4) {
    mu <- rep(0, d)
    df <- 3
    sigma <- matrix(0.5, nrow = d, ncol = d)
    diag(sigma) <- rep(1, d)
    X <-mvtnorm::rmvt(N, sigma, df, delta = mu)
    prob.Kplus2 <- exp( X %*% beta)
    prob.Kplus2 <- prob.Kplus2 / rowSums(prob.Kplus2)
    Y <- apply(prob.Kplus2, 1, function(row) sample(0:K, size = 1, prob = row))
  } else if (case == 5) {
    mu <- rep(0, (d-1))
    index_diff <- abs(outer(1:(d-1), 1:(d-1), "-"))
    sigma <- 0.5 ^ index_diff
    # sigma <- matrix(0.5, nrow = (d-1), ncol = (d-1))
    # diag(sigma) <- rep(1, (d-1))
    X <- cbind(1, MASS::mvrnorm(N, mu, sigma))
    prob.Kplus2 <- exp( X %*% beta)
    prob.Kplus2 <- prob.Kplus2 / (rowSums(prob.Kplus2))
    Y <- apply(prob.Kplus2, 1, function(row) sample(0:K, size = 1, prob = row))
  }
  return(list(X=X, Y=Y, P.true=prob.Kplus2))
}
