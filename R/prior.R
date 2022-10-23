
# calculate the case control subsampling probabilities
# for logistic regression & softmax regression
proptional_ssp <- function(N, k, y){
  prob <- 1/((k+1) * as.data.frame(table(y))$Freq)
  return(prob[y+1])
}

# subsample with replacement
swr_indx <- function(N, r, pi){
  sample(1:N, r, T, pi)
}

# Poisson subsample
poisson_indx <- function(N, r, pi){
  u <- runif(N)
  c(1:N)[u < r*pi]
}


# compute coefficients for logistic regression
# with survey package
logistic_coef_estimate <- function(X, y, weights, indx){

  data <- as.data.frame(cbind(y, X)[indx,])
  formula <- paste(colnames(data)[1], "~",
                   paste(colnames(data)[-1], collapse = "+"), "-1")
  design <- survey::svydesign(ids=~1,
                              weights=~weights,
                              data = data)
  beta <- survey::svyglm(as.formula(formula),
                         design=design,
                         family=quasibinomial())$coefficients
  beta
}

softmax_coef_estimate <- function(X, y, weights, indx){
  fit <- nnet::multinom(y[indx] ~ X[indx,] - 1, weights = weights)
  t(coef(fit))
}





