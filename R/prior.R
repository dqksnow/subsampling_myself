
# calculate the case control subsampling probabilities
# for logistic regression & softmax regression
proptional_ssp <- function(N, k, y){
  I <- matrix(0, nrow = N, ncol = k)
  I[cbind(seq_along(y), y)] <- 1

  prob <- 1/((k+1) * as.data.frame(table(y))$Freq)
  return(prob[y+1])
}

# sample observations with replacement with
# case control subsampling probabilities
# for logistic regression & softmax regression
swr_indx <- function(N, r, prop){
  return(sample(1:N, r, T, prop))
}


# compute coefficients for logistic regression
# with survey package
logistic_coef_estimate <- function(X, y, pinv, indx){

  data <- as.data.frame(cbind(y, X)[indx,])
  formula <- paste(colnames(data)[1], "~",
                   paste(colnames(data)[-1], collapse = "+"))
  design <- survey::svydesign(ids=~1,
                              weights=~pinv,
                              data = data)
  beta <- survey::svyglm(as.formula(formula),
                         design=design,
                         family=quasibinomial())$coefficients
  return(beta)
}






