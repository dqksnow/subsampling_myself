rm(list=ls())
# Newton method for poisson regression
getMLE <- function(x, y, w) {
  beta <- rep(0, d)
  loop  <- 1
  Loop  <- 100
  msg <- "NA"
  while (loop <= Loop) {
    pr <- c(exp(x %*% beta))
    H <- t(x) %*% (pr * w * x)
    S <- colSums((y - pr) * w * x)
    tryCatch(
      {shs <- NA
      shs <- solve(H, S) }, 
      error=function(e){
        cat("\n ERROR :", loop, conditionMessage(e), "\n")})
    if (is.na(shs[1])) {
      msg <- "Not converge"
      beta <- loop <- NA
      break
    }
    beta.new <- beta + shs
    tlr  <- sum((beta.new - beta)^2)
    beta  <- beta.new
    if(tlr < 0.000001) {
      msg <- "Successful convergence"
      break
    }
    if (loop == Loop)
      warning("Maximum iteration reached")
    loop  <- loop + 1
  }
  list(par=beta, message=msg, iter=loop, Infor = H)
}

library(foreach)
## sample code
#foreach(a=list(X[1:5000,],X[5001:10000,]),b=list(Y[1:5000],Y[5001:10000]))%dopar%getMLE(a,b,1)


## r1 <- r1+r; r2=0
AlgTwoStpsimp <- function(X,Y,r1=r1,rho=0.01,beta=0,Psi=0){
  n <- nrow(X)
  idx <- 1:n
  pi <- rep((r1/n),n)
  decision <- rbinom(n,rep(1,n),prob=pi)
  idx.simp <- idx[decision==1]
  X.simp <- X[idx.simp,]
  Y.simp <- Y[idx.simp]
  pinv.simp <- n/r1
  fit.simp <- getMLE(x=X.simp, y=Y.simp, w=pinv.simp)
  msg <- fit.simp$message
  if (msg != "Successful convergence") {
    warning(paste("not converge", msg))
    beta.simp <- rep(NA,ncol(X))
  } 
  beta.simp <- fit.simp$par
  hatpsidot <- exp(X.simp %*% beta.simp)
  V_c <- t(X.simp)%*%diag(c((Y.simp-hatpsidot)^2*(pinv.simp^2)))%*%X.simp-t(X.simp)%*%diag(c((Y.simp-hatpsidot)^2*(pinv.simp)))%*%X.simp
  return(list(simp=beta.simp, msg=msg, infor = fit.simp$Infor, V_c=V_c))
}


AlgTwoStpmse <- function(X,Y,r1=r1,rho=0.01,beta,W.prop,Psi){
  n <- nrow(X)
  psi.dot  <-  ( exp(X %*% beta))
  idx <- 1:n
  PI.mMSE <- abs(Y - psi.dot) * rowSums((X%*%W.prop)^2)
  PI.mMSE <- PI.mMSE / (n*Psi)
  PI.mMSE1 <- abs((1-rho)*PI.mMSE+rho/n+1/(r1+1e-6))/2-abs((1-rho)*PI.mMSE+rho/n-1/(r1+1e-6))/2
  decision <- rbinom(n,rep(1,n),prob=r1*PI.mMSE1)
  idx.mMSE <- idx[decision==1]     
  X.mMSE <- X[idx.mMSE,]
  Y.mMSE <- Y[idx.mMSE]
  pinv.mMSE <- 1 / (r1*PI.mMSE1[idx.mMSE])
  fit.mMSE <- getMLE(x=X.mMSE, y=Y.mMSE, w=pinv.mMSE)
  msg <- fit.mMSE$message
  if (msg != "Successful convergence") {
    warning(paste("not converge", msg))
    beta.mMSE <- rep(NA,ncol(X))
  } 
  beta.mMSE <- fit.mMSE$par
  hatpsidot <- exp(X.mMSE %*% beta.mMSE)
  V_c <- t(X.mMSE)%*%diag(c((Y.mMSE-hatpsidot)^2*(pinv.mMSE^2)))%*%X.mMSE- t(X.mMSE)%*%diag(c((Y.mMSE-hatpsidot)^2*(pinv.mMSE)))%*%X.mMSE
  return(list(mMSE=beta.mMSE, msg=msg, infor = fit.mMSE$Infor, V_c = V_c))
}


AlgTwoStpmvc <- function(X,Y,r1=r1,rho=0.01,beta,Psi){
  n <- nrow(X)
  psi.dot  <-  ( exp(X %*% beta))
  idx <- 1:n
  PI.mVc <- abs(Y - psi.dot) * rowSums(X^2)
  PI.mVc <- PI.mVc / (n*Psi)
  PI.mVc1 <- abs((1-rho)*PI.mVc+rho/n+1/(r1+1e-6))/2-abs((1-rho)*PI.mVc+rho/n-1/(r1+1e-6))/2
  decision <- rbinom(n,rep(1,n),prob=r1*PI.mVc1)
  idx.mVc <- idx[decision==1]     
  X.mVc <- X[idx.mVc,]
  Y.mVc <- Y[idx.mVc]
  pinv.mVc <- 1 / (r1*PI.mVc1[idx.mVc])
  fit.mVc <- getMLE(x=X.mVc, y=Y.mVc, w=pinv.mVc)
  msg <- fit.mVc$message
  if (msg != "Successful convergence") {
    warning(paste("not converge", msg))
    beta.mVc <- rep(NA,ncol(X))
  } 
  beta.mVc <- fit.mVc$par
  hatpsidot <- exp(X.mVc %*% beta.mVc)
  V_c <- t(X.mVc)%*%diag(c((Y.mVc-hatpsidot)^2*(pinv.mVc^2)))%*%X.mVc-t(X.mVc)%*%diag(c((Y.mVc-hatpsidot)^2*(pinv.mVc)))%*%X.mVc
  return(list(mVc=beta.mVc, msg=msg, infor = fit.mVc$Infor, V_c=V_c))
}


AlgTwoStp <- function(X,Y,r1=r1, r2=r2,rho=0.01,K=1,method="simp") {
  if (method == "simp") {
    pilot <- AlgTwoStpsimp(X,Y,r1)
    n <- nrow(X)
    m <- floor(n/K)
    a <- b <- c()
    for(i in 1:K){
      if(i==K){
        a <- c(a,list(X[((K-1)*m+1):nrow(X),]))
        b <- c(b,list(Y[((K-1)*m+1):nrow(X)]))
      }
      else{
        a <- c(a,list(X[((i-1)*m+1):(i*m),]))
        b <- c(b,list(Y[((i-1)*m+1):(i*m)]))
      }
    }
    result <- foreach(a=a,b=b)%do%AlgTwoStpsimp(a,b,r2/K)
    H <- pilot$infor
    B <- pilot$infor%*%pilot$simp
    V_c <- pilot$V_c
    for(i in 1:K){
      H <- H + result[[i]]$infor
      B <- B + result[[i]]$infor%*%result[[i]]$simp
      V_c <- V_c + result[[i]]$V_c
    }
    if(any(is.na(B))){
      msg <- "not converge"
      beta.simp <- rep(NA,ncol(X))
      sigma_simp <- diag(solve(H)%*%V_c%*%solve(H))[2]
      return(list(simp=beta.simp, msg=msg,sigma=sigma_simp))
    }
    beta.simp <- solve(H)%*%B
    msg <- "Successful convergence"
    sigma_simp <- diag(solve(H)%*%V_c%*%solve(H))[2]
    return(list(simp=beta.simp, msg=msg,sigma=sigma_simp))
  }
  if (method != "simp") {
    n <- nrow(X)
    idx <- 1:n
    idx.prop <- sample(1:n, r1, T)
    x.prop <- X[idx.prop,]
    y.prop <- Y[idx.prop]
    pinv.prop <- rep(n/r1,r1)
    fit.prop <- getMLE(x=x.prop, y=y.prop, w=pinv.prop)
    beta.prop <- fit.prop$par
    if (is.na(beta.prop[1]))
      return(list(opt=NA, msg="first stage not converge"))
    psi.dot  <-  ( exp(X %*% beta.prop))
    hatpsidot <- exp(x.prop%*%fit.prop$par)
    V_c0 <- t(x.prop)%*%diag(c((y.prop-hatpsidot)^2*(pinv.prop^2)))%*%x.prop
    
    ## mMSE
    psi.ddot <- psi.dot[idx.prop]
    W.prop <- solve(t(x.prop) %*% (x.prop * psi.ddot * pinv.prop))
    Psi <- sum(abs(Y[idx.prop] - psi.ddot) * rowSums((X[idx.prop,]%*%W.prop)^2) )/r1
    m <- floor(n/K)
    a <- b <- c()
    for(i in 1:K){
      if(i==K){
        a <- c(a,list(X[((K-1)*m+1):nrow(X),]))
        b <- c(b,list(Y[((K-1)*m+1):nrow(X)]))
      }
      else{
        a <- c(a,list(X[((i-1)*m+1):(i*m),]))
        b <- c(b,list(Y[((i-1)*m+1):(i*m)]))
      }
    }
    result <- foreach(a=a,b=b)%do%AlgTwoStpmse(a,b,r2/K,rho,beta.prop,W.prop,Psi)
    H <- fit.prop$Infor
    B <- fit.prop$Infor%*%fit.prop$par
    V_cmse <- V_c0
    for(i in 1:K){
      H <- H + result[[i]]$infor
      B <- B + result[[i]]$infor%*%result[[i]]$mMSE
      V_cmse <- V_cmse + result[[i]]$V_c
    }
    if(any(is.na(B))){
      msg.mMSE <- "not converge"
      beta.mMSE <- rep(NA,ncol(X))
      sigma_mMSE <- NA
    }
    else{
      beta.mMSE <- solve(H)%*%B
      sigma_mMSE <- diag(solve(H)%*%V_cmse%*%solve(H))[2]
      msg.mMSE <- "Successful convergence"
    }
    
    
    ## mVc
    PI.mVc <- abs(Y - psi.dot) * rowSums(X^2)
    Psi <- sum(abs(Y[idx.prop] - psi.ddot) * rowSums((X[idx.prop,])^2))/r1
    m <- floor(n/K)
    a <- b <- c()
    for(i in 1:K){
      if(i==K){
        a <- c(a,list(X[((K-1)*m+1):nrow(X),]))
        b <- c(b,list(Y[((K-1)*m+1):nrow(X)]))
      }
      else{
        a <- c(a,list(X[((i-1)*m+1):(i*m),]))
        b <- c(b,list(Y[((i-1)*m+1):(i*m)]))
      }
    }
    result <- foreach(a=a,b=b)%do%AlgTwoStpmvc(a,b,r2/K,rho,beta.prop,Psi)
    H <- fit.prop$Infor
    B <- fit.prop$Infor%*%fit.prop$par
    V_cmv <- V_c0
    for(i in 1:K){
      H <- H + result[[i]]$infor
      B <- B + result[[i]]$infor%*%result[[i]]$mVc
      V_cmv <- V_cmv + result[[i]]$V_c
    }
    if(any(is.na(B))){
      msg.mVc <- "not converge"
      beta.mVc <- rep(NA,ncol(X))
      sigma_mvc <- NA
    }
    else{
      beta.mVc <- solve(H)%*%B
      sigma_mvc <- diag(solve(H)%*%V_cmv%*%solve(H))[2]
      msg.mVc <- "Successful convergence" 
    }
    
    opt <- cbind(beta.mMSE,  beta.mVc)
    msg <- c(msg.mMSE,  msg.mVc)
    sigma <- c(sigma_mMSE, sigma_mvc)
    
    return(list(opt=opt, msg=msg, sigma =sigma))
  }
}

### test code ###
set.seed(1234)
n<-5e5;
beta0  <- c(rep(1/2, 7))
d <- length(beta0)
X  <- matrix(runif(n*d, 0, 1),n,d)
#X[,2] <- X[,1]+runif(n,0,0.1)
#X[,6] <- runif(n,-1,1)
#X[,7] <- runif(n,-1,1)
lambda  <- exp(X %*% beta0)
Y  <- rpois(n, lambda)
print(mean(Y))

fit.full <- getMLE(X, Y, 1)

#AlgTwoStp(X,Y,200,1000)
#AlgTwoStp(X,Y,200,500,rho=0.2,K=25,method="opt")

rpt <- 1000
lrs <- 1
Beta.simp <- matrix(NA, d, rpt*lrs)
nm <- 2
Beta.opt <- matrix(NA, d*nm, rpt*lrs)
Sigma.est <- matrix(NA,nm,rpt)
Sigma.simp <- matrix(NA,1,rpt)
coverrate <- rep(NA,nm+1)
conflen <- matrix(NA, nm+1, rpt)


r1<-200


for(k in c(1,5)){
cat("\n K=",k,"-")
  for(r in c(1000,1500)*k){
cat(r,":\n")

set.seed(0)
for(i in 1:rpt){
  if (i%/%100 == i/100) cat(i, "-")
  fit.alg <- AlgTwoStp(X,Y,r1,r,rho=0.2,K=k,method="opt")
  Beta.opt[,i] <- c(fit.alg$opt)
  Sigma.est[,i] <- c(fit.alg$sigma)
}

set.seed(0)
for(i in 1:rpt){
  if (i%/%100 == i/100) cat(i, "-")
  fit.alg <- AlgTwoStp(X,Y,r1,r,rho=0.2,K=k,method="simp")
  Beta.simp[,i] <- c(fit.alg$simp)
  Sigma.simp[,i] <- c(fit.alg$sigma)
}


for(j in 1:nm){
  cover <- 0
  for (i in 1:rpt){
    interval1 <- Beta.opt[7*j-5,i]-1.96*sqrt(Sigma.est[j,i])
    interval2 <- Beta.opt[7*j-5,i]+1.96*sqrt(Sigma.est[j,i])
    if(0.5>=interval1 & 0.5<=interval2){
      cover <- cover + 1
    }
    conflen[j,i] <- interval2-interval1
  }
  coverrate[j] <- cover/rpt
}

cover <- 0
for (i in 1:rpt){
  interval1 <- Beta.simp[2,i]-1.96*sqrt(Sigma.simp[1,i])
  interval2 <- Beta.simp[2,i]+1.96*sqrt(Sigma.simp[1,i])
  if(0.5>=interval1 & 0.5<=interval2){
    cover <- cover + 1
  }
  conflen[nm+1,i] <- interval2-interval1
}
coverrate[nm+1] <- cover/rpt

cat("coverrate:\n")
cat(coverrate)
cat("conflen:\n")
cat(rowMeans(conflen))
cat("\n")
  }
}

