
# calculate the case control subsampling probabilities
# for logistic regression & softmax regression
proptional_ssp <- function(N, k, y){
  prob <- 1/((k+1) * as.data.frame(table(y))$Freq)
  return(prob[y+1])
}



# for logistic regression
# pilot sample subsampling without replacement
# what if sum(y == 1) < r/2 ???
plt_indx <- function(N, y, r){
  c(sample(c(1:N)[y == 0], r/2, replace = F),
    sample(c(1:N)[y == 1], r/2, replace = F))
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
# logistic.coef.estimate <-
#   function(X,
#            Y,
#            offset = rep(0, length(Y)),
#            start = rep(0, ncol(X)),
#            weights = rep(1, length(Y))) {
#     if (length(weights) == 1)
#       weights <- rep(weights, length(Y))
#     if (offset != rep(0, length(Y)) & weights != rep(1, length(Y))) {
#       stop("offset will only be needed in the unweighted likelihood function.")
#     }
#     data <- as.data.frame(cbind(Y, X))
#     formula <- paste(colnames(data)[1], "~",
#                      paste(colnames(data)[-1], collapse = "+"), "-1")
#     design <- survey::svydesign(ids =  ~ 1,
#                                 weights =  ~ weights,
#                                 data = data)
#     fit <- try(
#       beta <- survey::svyglm( as.formula(formula),
#                               design = design,
#                               offset = offset,
#                               start = start,
#                               family = quasibinomial(link = "logit"))$coefficients,
#       silent=TRUE
#     )
#     if ("try-error" %in% class(fit)){
#       message("Warning: an error occurred while calling 'glm.fit': ", geterrmessage(),
#               "This is probably due to the iteratively reweighted least squares algorithm called by 'glm.fit' not getting converged. Tring another function 'getMSLE' to replace 'glm.fit'.")
#       beta <- getMSLE(X, Y, offset = offset, start = start)$beta
#     }
#     return(list(beta = beta))
#   }

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





