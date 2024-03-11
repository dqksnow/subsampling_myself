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
  pinv.simp <- n
  fit.simp <- getMLE(x=X.simp, y=Y.simp, w=pinv.simp)
  msg <- fit.simp$message
  if (msg != "Successful convergence") {
    warning(paste("not converge", msg))
    beta.simp <- rep(NA,ncol(X))
  } 
  beta.simp <- fit.simp$par
  return(list(simp=beta.simp, msg=msg, infor = fit.simp$Infor))
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
  pinv.mMSE <- 1 / PI.mMSE1[idx.mMSE]
  fit.mMSE <- getMLE(x=X.mMSE, y=Y.mMSE, w=pinv.mMSE)
  msg <- fit.mMSE$message
  if (msg != "Successful convergence") {
    warning(paste("not converge", msg))
    beta.mMSE <- rep(NA,ncol(X))
  } 
  beta.mMSE <- fit.mMSE$par
  return(list(mMSE=beta.mMSE, msg=msg, infor = fit.mMSE$Infor))
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
  pinv.mVc <- 1 / PI.mVc1[idx.mVc]
  fit.mVc <- getMLE(x=X.mVc, y=Y.mVc, w=pinv.mVc)
  msg <- fit.mVc$message
  if (msg != "Successful convergence") {
    warning(paste("not converge", msg))
    beta.mVc <- rep(NA,ncol(X))
  } 
  beta.mVc <- fit.mVc$par
  return(list(mVc=beta.mVc, msg=msg, infor = fit.mVc$Infor))
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
    for(i in 1:K){
      H <- H + result[[i]]$infor
      B <- B + result[[i]]$infor%*%result[[i]]$simp
    }
    if(any(is.na(B))){
      msg <- "not converge"
      beta.simp <- rep(NA,ncol(X))
      return(list(simp=beta.simp, msg=msg))
    }
    beta.simp <- solve(H)%*%B
    msg <- "Successful convergence"
    return(list(simp=beta.simp, msg=msg))
  }
  if (method != "simp") {
    n <- nrow(X)
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
    for(i in 1:K){
      H <- H + result[[i]]$infor
      B <- B + result[[i]]$infor%*%result[[i]]$mMSE
    }
    if(any(is.na(B))){
      msg.mMSE <- "not converge"
      beta.mMSE <- rep(NA,ncol(X))
    }
    else{
      beta.mMSE <- solve(H)%*%B
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
    for(i in 1:K){
      H <- H + result[[i]]$infor
      B <- B + result[[i]]$infor%*%result[[i]]$mVc
    }
    if(any(is.na(B))){
      msg.mVc <- "not converge"
      beta.mVc <- rep(NA,ncol(X))
    }
    else{
      beta.mVc <- solve(H)%*%B
      msg.mVc <- "Successful convergence" 
    }
    
    opt <- cbind(beta.mMSE,  beta.mVc)
    msg <- c(msg.mMSE,  msg.mVc)
    
    return(list(opt=opt, msg=msg))
  }
}

### test code ###
set.seed(1234)
n<-5e5;
beta0  <- c(rep(1/2, 7))
d <- length(beta0)
X  <- matrix(runif(n*d, 0, 1),n,d)
# X[,2] <- X[,1]+runif(n,0,1)
# X[,6] <- runif(n,-1,1)
# X[,7] <- runif(n,-1,1)
lambda  <- exp(X %*% beta0)
Y  <- rpois(n, lambda)
print(mean(Y))

fit.full <- getMLE(X, Y, 1)

#AlgTwoStp(X,Y,200,1000)
#AlgTwoStp(X,Y,200,500,rho=0.2,K=25,method="opt")

set.seed(0)
rpt <- 1000
r.ss <- c(300,500,700,1000,1200,1500,1700,2000)## [5:6]
lrs <- length(r.ss)
Beta.simp <- matrix(NA, d, rpt*lrs)
nm <- 2
Beta.opt <- matrix(NA, d*nm, rpt*lrs)


k <- 1
r1<-200
itr <- 0
for (r in r.ss) {
  cat(r, " ",k," ");
  set.seed(0)
  for (rr in 1:rpt) {
    if (rr%/%100 == rr/100) cat(rr, "-")
    itr <- itr+1
    fit.alg.simp <- AlgTwoStp(X,Y,r1,r,rho=0.2,K=k,method="simp")
    Beta.simp[,itr] <- fit.alg.simp$simp
    if (any(fit.alg.simp$msg != "Successful convergence")) {
      warning(paste(rr, "not converge", fit.alg.simp$msg))
      Beta.simp[,itr] <- NA
    }
    fit.alg <- AlgTwoStp(X,Y,r1,r,rho=0.2,K=k,method="opt")
    if (any(fit.alg$msg != "Successful convergence")) {
      warning(paste(rr, "not converge:", fit.alg$msg))
      next
    }
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


write.csv(rec,"D:/report5/data/case1-new-k1.csv")

plot(r.ss, log(rec[,1]), xlab="r", ylab="MSE", pch="1",lwd=2,
     ylim=c(min(log(rec)), max(log(rec))), type="o")
pp <- dim(rec)[2]
for(i in 2:pp)
  lines(r.ss, log(rec[,i]), type="o", pch=paste(i), lty=i, col=i,lwd=2)
legend("topright", lty=1:pp, pch=paste(1:pp), col=1:pp,
       legend=c("uniform", "MV", "MVC"),cex=1.5)



Deff <- data.frame(r.ss=rep(r.ss,3),mse=log(c(rec)),method=rep(c("UNIF","MV","MVc"),
                                                               each=length(r.ss)))

library(ggplot2)#275
theme_set(theme_bw())
b <- qplot(x=r.ss, y=mse, data=Deff,  geom=c("line", "point"), shape=method, main=" ",
           xlab="K",ylab="log(MSE)",ylim = c(min(Deff$mse), max(Deff$mse)),color=method)+geom_point(size = 3.5)
b+theme(legend.title = element_text(colour="blue", size=16),
        legend.text = element_text(size = 16, face = "bold"))




# layout(cbind(1,2), widths =c(6,1.5))  # put legend on bottom 1/8th of the chart
# par(mar=c(4.5, 3, 2, 1))
# plot(rho.ss, rec[,1], xlab=expression(rho), ylab="MSE", pch="1",lwd=2,
#      ylim=c(min(rec), max(rec)), type="o")
# pp <- dim(rec)[2]
# for(i in 2:pp)
#   lines(rho.ss, rec[,i], type="o", pch=paste(i), lty=i, col=i,lwd=2)
# 
# # setup for no margins on the legend
# par(mar=c(0, 0, 0, 0))
# # c(bottom, left, top, right)
# plot.new()
# legend("center",'groups', lty=1:pp, pch=paste(1:pp), col=1:pp,
#        legend=c("UNIF", "MV", "MVC"),cex=1.3,bty='n',horiz = F)


# plot(rho.ss, rec[,1], xlab=expression(rho), ylab="MSE", pch="1",lwd=2,
#      ylim=c(min(rec), max(rec)), type="o")
# pp <- dim(rec)[2]
# for(i in 2:pp)
#   lines(rho.ss, rec[,i], type="o", pch=paste(i), lty=i, col=i,lwd=2)
# legend("right", lty=1:pp, pch=paste(1:pp), col=1:pp,
#        legend=c("uniform", "MV", "MVC"),cex=1.3)

