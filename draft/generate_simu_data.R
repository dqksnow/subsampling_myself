library(mvtnorm)

#### JASA2018 simu1 multinormal balanced ########
simu_mzNormal <- function(seed, N, beta0, corr){
  d <- length(beta0)
  sigmax  <- matrix(corr, d, d) + diag(1-corr, d)
  set.seed(seed)
  X <- rmvnorm(N, rep(0, d), sigmax)
  #X <- cbind(1, X)
  P <- 1 - 1 / (1 + exp(X %*% beta0))
  y <- rbinom(N, 1, P)
  beta_full <- logistic_coef_estimate(X, y, 1, 1:length(y))
  return(list(X=X, y=y, beta_full=beta_full))
}

#### JASA2018 simu2 multinormal imbalanced ########
simu_nzNormal <- function(seed, N, beta0, corr){
  d <- length(beta0)
  sigmax  <- matrix(corr, d, d) + diag(1-corr, d)
  set.seed(seed)
  X <- rmvnorm(N, rep(1.5, d), sigmax)
  #X <- cbind(1, X)
  P <- 1 - 1 / (1 + exp(X %*% beta0))
  y <- rbinom(N, 1, P)
  beta_full <- logistic_coef_estimate(X, y, 1, 1:length(y))
  return(list(X=X, y=y, beta_full=beta_full))
}

#### JASA2018 simu3 hetero multinormal ########
simu_ueNormal <- function(seed, N, beta0, corr){
  d <- length(beta0)
  sigmax <- matrix(0,d,d)
  var <- 1/c(1:d)^2 #
  #var <- c(1:(d)) #^2
  for (i in 1:d){
    for (j in 1:d){
      sigmax[i,j] <- corr^(i!=j) * (sqrt(var[i])*sqrt(var[j]))
    }
    }
  #sigmax <- matrix(corr, d, d)  + diag(as.vector(var)) - diag(corr, d)
  X <- rmvnorm(N, rep(0, d), sigmax)
  #X <- cbind(1, X)
  P <- 1 - 1 / (1 + exp(X %*% beta0))
  y <- rbinom(N, 1, P)
  beta_full <- logistic_coef_estimate(X, y, 1, 1:length(y))
  return(list(X=X, y=y, beta_full=beta_full))
}

#### JASA2018 simu4 biomodal normal ########
simu_mixNormal <- function(seed, N, beta0, corr){
  d <- length(beta0)
  sigmax  <- matrix(corr, d, d) + diag(1-corr, d)
  set.seed(seed)
  X <- 0.5*rmvnorm(N, rep(1, d), sigmax) + 0.5*rmvnorm(N, rep(-1, d), sigmax)
  #X <- cbind(1, X)
  P <- 1 - 1 / (1 + exp(X %*% beta0))
  y <- rbinom(N, 1, P)
  beta_full <- logistic_coef_estimate(X, y, 1, 1:length(y))
  return(list(X=X, y=y, beta_full=beta_full))
}

#### JASA2018 simu5 multi t ########
simu_T3 <- function(seed, N, beta0, corr, df){
  d <- length(beta0)
  sigmax  <- matrix(corr, d, d) + diag(1-corr, d)
  set.seed(seed)
  X <-rmvt(N, sigmax, df, rep(0, d))/10#, type=c("shifted","Kshirsagar")[1]
  #X <- cbind(1, X)
  P <- 1 - 1 / (1 + exp(X %*% beta0))
  y <- rbinom(N, 1, P)
  beta_full <- logistic_coef_estimate(X, y, 1, 1:length(y))
  return(list(X=X, y=y, beta_full=beta_full))
}

#### JASA2018 simu6 exp ########
simu_EXP <- function(seed, N, beta0, rate){
  d <- length(beta0)
  set.seed(seed)
  X <- matrix(0, N, d)
  generate_rexp <- function(x) x<-rexp(N, rate=rate)
  X <- apply(X, 2, generate_rexp)
  #X <- cbind(1, X)
  P <- 1 - 1 / (1 + exp(X %*% beta0))
  y <- rbinom(N, 1, P)
  beta_full <- logistic_coef_estimate(X, y, 1, 1:length(y))
  return(list(X=X, y=y, beta_full=beta_full))
}

#
predict_classification <- function(X, y, beta){
  pbeta_estimate <- pbeta(X, beta)
  y_pred <- rep(0, length(y))
  y_pred[which(pbeta_estimate>0.5)] <- 1
  proportion <- sum(y_pred==y)/length(y)
  return(proportion)
}
