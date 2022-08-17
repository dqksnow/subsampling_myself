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
  list(par=beta, message=msg, iter=loop)
}

cal <- function(i,PI1,r){
  return((r-i+1)*PI1[i]/sum(PI1[-c(1:i)]))
}
findM <- function(PI,r){
  PI1 <- sort(PI,decreasing = T)
  if(r*PI1[1]/sum(PI1)<=1){
    return(0)
  }
  index <- 2:r
  rough <- ((r-1):1)*PI1[2:r]/rev(cumsum(rev(PI1)))[2:r]
  return(min(index[rough<=1]))
}




## r1 <- r1+r; r2=0

AlgTwoStp <- function(r1=r1, r2=r2) {
  if (r2 == 0) {
    idx <- 1:n
    pi <- rep((r1/n),n)
    decision <- rbinom(n,rep(1,n),prob=pi)
    idx.simp <- idx[decision==1]
    x.simp <- X[idx.simp,]
    y.simp <- Y[idx.simp]
    pinv.simp <- n
    fit.simp <- getMLE(x=x.simp, y=y.simp, w=pinv.simp)
    beta.simp <- fit.simp$par
    msg <- fit.simp$message
    return(list(simp=beta.simp, msg=msg))
  }
  if (r2 != 0) {
    idx <- 1:n
    idx.prop <- sample(1:n, r1, T)
    x.prop <- X[idx.prop,]
    y.prop <- Y[idx.prop]
    pinv.prop <- rep(n,r1)
    fit.prop <- getMLE(x=x.prop, y=y.prop, w=pinv.prop)
    beta.prop <- fit.prop$par
    if (is.na(beta.prop[1]))
      return(list(opt=NA, msg="first stage not converge"))
    psi.dot  <-  ( exp(X %*% beta.prop))
    psi.ddot <- psi.dot[idx.prop]
    
    
    ## mVc,rho=0.001,M opt
    rho = 0.001
    PI.mVc <-PI.mVc2 <- abs(Y - psi.dot) * rowSums(X^2)
    M <- findM(PI.mVc,r2)
    if(M==0){
      PI.mVc2 <- r2*PI.mVc/sum(PI.mVc)
    }
    else{
      P1 <- sort(PI.mVc,decreasing = T)[M] 
      PI.mVc2[PI.mVc<=P1] <- (r2-M+1)*PI.mVc2[PI.mVc<=P1]/sum(PI.mVc2[PI.mVc<=P1])
      PI.mVc2[PI.mVc>P1] <- 1
    }
    PI.mVc <- (1-rho)*PI.mVc2 +rho*r2/n
    PI.mVc1 <- abs(PI.mVc+1-1e-6)/2-abs(PI.mVc-1+1e-6)/2
    decision <- rbinom(n,rep(1,n),prob=PI.mVc1)
    idx.mVc <- idx[decision==1]
    x.mVc <- X[c(idx.mVc, idx.prop),]
    y.mVc <- Y[c(idx.mVc, idx.prop)]
    pinv.mVc <- c(r2 / PI.mVc1[idx.mVc], pinv.prop)
    fit.mVc <- getMLE(x=x.mVc, y=y.mVc, w=pinv.mVc)
    
    ## mVc,rho=0.2
    rho = 0.2
    PI.mVc <- abs(Y - psi.dot) * rowSums(X^2)
    Psi <- sum(abs(Y[idx.prop] - psi.ddot) * rowSums((X[idx.prop,])^2))/r1
    PI.mVc <- (1-rho)*PI.mVc / (n*Psi)+rho/n
    PI.mVc1 <- abs(r2*PI.mVc+1-1e-6)/2-abs(r2*PI.mVc-1+1e-6)/2
    decision <- rbinom(n,rep(1,n),prob=PI.mVc1)
    idx.mVc <- idx[decision==1]
    x.mVc <- X[c(idx.mVc, idx.prop),]
    y.mVc <- Y[c(idx.mVc, idx.prop)]
    pinv.mVc <- c(r2 / PI.mVc1[idx.mVc], pinv.prop)
    fit.mVc1 <- getMLE(x=x.mVc, y=y.mVc, w=pinv.mVc)
    
    ## mVc, hard thresholding, rho=0.2
    rho = 0.2
    PI.mVc <- abs(Y - psi.dot) * rowSums(X^2)
    PI.prop <- abs(Y[idx.prop] - psi.ddot) * rowSums((X[idx.prop,])^2)
    M <- quantile(PI.prop,(1-r/n*0.5))
    PI.prop <- abs(PI.prop+M)/2-abs(PI.prop-M)/2
    Psi <- sum(PI.prop)/r1
    PI.mVc <- abs(PI.mVc+M)/2-abs(PI.mVc-M)/2
    PI.mVc <- (1-rho)*PI.mVc / (n*Psi)+rho/n
    PI.mVc1 <- abs(r2*PI.mVc+1-1e-6)/2-abs(r2*PI.mVc-1+1e-6)/2
    decision <- rbinom(n,rep(1,n),prob=PI.mVc1)
    idx.mVc <- idx[decision==1]
    x.mVc <- X[c(idx.mVc, idx.prop),]
    y.mVc <- Y[c(idx.mVc, idx.prop)]
    pinv.mVc <- c(r2 / PI.mVc1[idx.mVc], pinv.prop)
    fit.mVc2 <- getMLE(x=x.mVc, y=y.mVc, w=pinv.mVc)
    
    ## mVc, optimal M, rho=0.2
    rho = 0.2
    PI.mVc <-PI.mVc2 <- abs(Y - psi.dot) * rowSums(X^2)
    M <- findM(PI.mVc,r2)
    if(M==0){
      PI.mVc2 <- r2*PI.mVc/sum(PI.mVc)
    }
    else{
      P1 <- sort(PI.mVc,decreasing = T)[M] 
      PI.mVc2[PI.mVc<=P1] <- (r2-M+1)*PI.mVc2[PI.mVc<=P1]/sum(PI.mVc2[PI.mVc<=P1])
      PI.mVc2[PI.mVc>P1] <- 1
    }
    PI.mVc <- (1-rho)*PI.mVc2 +rho*r2/n
    PI.mVc1 <- abs(PI.mVc+1-1e-6)/2-abs(PI.mVc-1+1e-6)/2
    decision <- rbinom(n,rep(1,n),prob=PI.mVc1)
    idx.mVc <- idx[decision==1]
    x.mVc <- X[c(idx.mVc, idx.prop),]
    y.mVc <- Y[c(idx.mVc, idx.prop)]
    pinv.mVc <- c(r2 / PI.mVc1[idx.mVc], pinv.prop)
    fit.mVc3 <- getMLE(x=x.mVc, y=y.mVc, w=pinv.mVc)
    
    ## mMSE,rho = 0.001,M opt
    rho = 0.001
    W.prop <- solve(t(x.prop) %*% (x.prop * psi.ddot * pinv.prop))
    PI.mMSE <- PI.mMSE2 <- abs(Y - psi.dot) * rowSums((X%*%W.prop)^2)
    M <- findM(PI.mMSE,r2)
    if(M==0){
      PI.mMSE2 <- r2*PI.mMSE/sum(PI.mMSE)
    }
    else{
      P1 <- sort(PI.mMSE,decreasing = T)[M] 
      PI.mMSE2[PI.mMSE<=P1] <- (r2-M+1)*PI.mMSE2[PI.mMSE<=P1]/sum(PI.mMSE2[PI.mMSE<=P1])
      PI.mMSE2[PI.mMSE>P1] <- 1
    }
    PI.mMSE <- (1-rho)*PI.mMSE2  + rho*r2/n
    PI.mMSE1 <- abs(PI.mMSE+1-1e-6)/2-abs(PI.mMSE-(1-1e-6))/2
    decision <- rbinom(n,rep(1,n),prob=PI.mMSE1)
    idx.mMSE <- idx[decision==1]    
    x.mMSE <- X[c(idx.mMSE, idx.prop),]
    y.mMSE <- Y[c(idx.mMSE, idx.prop)]
    pinv.mMSE <- c(r2 / PI.mMSE1[idx.mMSE], pinv.prop)
    fit.mMSE <- getMLE(x=x.mMSE, y=y.mMSE, w=pinv.mMSE)
    
    
    ## mMSE,rho = 0.2
    rho = 0.2
    W.prop <- solve(t(x.prop) %*% (x.prop * psi.ddot * pinv.prop))
    PI.mMSE <- abs(Y - psi.dot) * rowSums((X%*%W.prop)^2)
    psi.ddot <- psi.dot[idx.prop]
    W.prop <- solve(t(x.prop) %*% (x.prop * psi.ddot * pinv.prop))
    Psi <- sum(abs(Y[idx.prop] - psi.ddot) * rowSums((X[idx.prop,]%*%W.prop)^2) )/r1
    PI.mMSE <- (1-rho)*PI.mMSE / (n*Psi) + rho/n
    PI.mMSE1 <- abs(r2*PI.mMSE+1-1e-6)/2-abs(r2*PI.mMSE-(1-1e-6))/2
    decision <- rbinom(n,rep(1,n),prob=PI.mMSE1)
    idx.mMSE <- idx[decision==1]    
    x.mMSE <- X[c(idx.mMSE, idx.prop),]
    y.mMSE <- Y[c(idx.mMSE, idx.prop)]
    pinv.mMSE <- c(r2 / PI.mMSE1[idx.mMSE], pinv.prop)
    fit.mMSE1 <- getMLE(x=x.mMSE, y=y.mMSE, w=pinv.mMSE)    
    
    ## mMSE, hard thresholding, rho=0.2
    rho = 0.2
    W.prop <- solve(t(x.prop) %*% (x.prop * psi.ddot * pinv.prop))
    PI.mMSE <- abs(Y - psi.dot) * rowSums((X%*%W.prop)^2)
    psi.ddot <- psi.dot[idx.prop]
    W.prop <- solve(t(x.prop) %*% (x.prop * psi.ddot * pinv.prop))
    PI.prop <- abs(Y[idx.prop] - psi.ddot) * rowSums((X[idx.prop,]%*%W.prop)^2)
    M <- quantile(PI.prop,(1-r/n*0.5))
    PI.prop <- abs(PI.prop+M)/2-abs(PI.prop-M)/2
    Psi <- sum(PI.prop)/r1
    PI.mMSE <- abs(PI.mMSE+M)/2-abs(PI.mMSE-M)/2
    PI.mMSE <- (1-rho)*PI.mMSE / (n*Psi) + rho/n
    PI.mMSE1 <- abs(r2*PI.mMSE+1-1e-6)/2-abs(r2*PI.mMSE-(1-1e-6))/2
    decision <- rbinom(n,rep(1,n),prob=PI.mMSE1)
    idx.mMSE <- idx[decision==1]    
    x.mMSE <- X[c(idx.mMSE, idx.prop),]
    y.mMSE <- Y[c(idx.mMSE, idx.prop)]
    pinv.mMSE <- c(r2 / PI.mMSE1[idx.mMSE], pinv.prop)
    fit.mMSE2 <- getMLE(x=x.mMSE, y=y.mMSE, w=pinv.mMSE) 
    
    ## mMSE, optimal M, rho=0.2
    rho = 0.2
    W.prop <- solve(t(x.prop) %*% (x.prop * psi.ddot * pinv.prop))
    PI.mMSE <- PI.mMSE2 <- abs(Y - psi.dot) * rowSums((X%*%W.prop)^2)
    M <- findM(PI.mMSE,r2)
    if(M==0){
      PI.mMSE2 <- r2*PI.mMSE/sum(PI.mMSE)
    }
    else{
      P1 <- sort(PI.mMSE,decreasing = T)[M] 
      PI.mMSE2[PI.mMSE<=P1] <- (r2-M+1)*PI.mMSE2[PI.mMSE<=P1]/sum(PI.mMSE2[PI.mMSE<=P1])
      PI.mMSE2[PI.mMSE>P1] <- 1
    }
    PI.mMSE <- (1-rho)*PI.mMSE2  + rho*r2/n
    PI.mMSE1 <- abs(PI.mMSE+1-1e-6)/2-abs(PI.mMSE-(1-1e-6))/2
    decision <- rbinom(n,rep(1,n),prob=PI.mMSE1)
    idx.mMSE <- idx[decision==1]    
    x.mMSE <- X[c(idx.mMSE, idx.prop),]
    y.mMSE <- Y[c(idx.mMSE, idx.prop)]
    pinv.mMSE <- c(r2 / PI.mMSE1[idx.mMSE], pinv.prop)
    fit.mMSE3 <- getMLE(x=x.mMSE, y=y.mMSE, w=pinv.mMSE)  
    
    
    opt <- cbind(fit.mMSE$par,fit.mMSE1$par,fit.mMSE2$par, fit.mMSE3$par, 
                 fit.mVc$par,fit.mVc1$par,fit.mVc2$par,fit.mVc3$par)
    msg <- c(fit.mMSE$message,fit.mMSE1$message,fit.mMSE2$message, fit.mMSE3$message, 
             fit.mVc$message,fit.mVc1$message,fit.mVc2$message,fit.mVc3$message)
    
    return(list(opt=opt))
  }
}

### test code ###
set.seed(1234)
n<-5e5;
beta0  <- c(rep(1/2, 7))
d <- length(beta0)
X  <- matrix(runif(n*d, 0, 1),n,d)
X[,2] <- X[,1]+runif(n,0,1)
X[,6] <- runif(n,-1,1)
X[,7] <- runif(n,-1,1)
lambda  <- exp(X %*% beta0)
Y  <- rpois(n, lambda)
print(mean(Y))

fit.full <- getMLE(X, Y, 1)

r1<-200
r.ss <- c(0.01,0.1,0.3,0.5,0.7)*n

set.seed(0)
rpt <- 1000
#r.ss <- c( 300, 500, 700, 1000, 1200, 1400)## [5:6]
lrs <- length(r.ss)
Beta.simp <- matrix(NA, d, rpt*lrs)
nm <- 8
Beta.opt <- matrix(NA, d*nm, rpt*lrs)
itr <- 0

for (r in r.ss) {
  cat(r, " ");
  set.seed(0)
  for (rr in 1:rpt) {
    if (rr%/%100 == rr/100) cat(rr, "-")
    itr <- itr+1
    fit.alg.simp <- AlgTwoStp(r1+r, 0)
    Beta.simp[,itr] <- fit.alg.simp$simp
    fit.alg <- AlgTwoStp(r1, r)
    Beta.opt[,itr] <- c(fit.alg$opt)
  }
  cat("\n\n")
}
l2.simp <- (Beta.simp - fit.full$par)^2
l2.opt  <- (Beta.opt  - fit.full$par)^2
mse.simp <- mse.opt <- c()
for (i in 1:lrs) {
  loc <- ((i-1)*rpt+1):(i*rpt)
  mse.simp <- c(mse.simp, sum(rowMeans(l2.simp[,loc], na.rm=TRUE)))
  mse.opt <- rbind(mse.opt, colSums(matrix(
    rowMeans(l2.opt[,loc], na.rm=TRUE), d, nm)))
}

rec <- cbind(mse.simp, mse.opt)


write.csv(Beta.simp,"D:/report5/data/largeR-case4-compare-betasimp.csv")
write.csv(Beta.opt,"D:/report5/data/largeR-case4-compare-betaopt.csv")
write.csv(rec,"D:/report5/data/largeR-case4-compare-rec.csv")