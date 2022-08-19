
library(mvtnorm)
n <- 1e4
beta0  <- c(rep(0.5, 7))
d <- length(beta0)
corr  <- 0.5
sigmax  <- matrix(corr, d-1, d-1) + diag(1-corr, d-1)

set.seed(123)
X <- rmvnorm(n, rep(0, d-1), sigmax)
X <- cbind(1, X)
P <- 1 - 1 / (1 + exp(X %*% beta0))
y <- rbinom(n, 1, P)



## formula
covariate <- ifelse(is.null(colnames(X)),
                    paste(paste0("V", 1:d), collapse = "+"),
                    paste(colnames(X), collapse = "+"))
rsp <- ifelse(is.null(names(y)), "y", names(y))
fmla <- ifelse(intercept,
               paste0(rsp, "~", covariate),
               paste0(rsp, "~", covariate, "-1"))



