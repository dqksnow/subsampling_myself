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





## r1 <- r1+r; r2=0

AlgTwoStpsimp <- function(r1=r1, r2=r2) {
  r1 <- r1+r2
  idx.simp <- sample(1:n, r1, T)
  x.simp <- X[idx.simp,]
  y.simp <- Y[idx.simp]
  pinv.simp <- n
  fit.simp <- getMLE(x=x.simp, y=y.simp, w=pinv.simp)
  beta.simp <- fit.simp$par
  msg <- fit.simp$message
  return(list(simp=beta.simp, msg=msg))
}

AlgTwoStpmse <- function(r1=r1, r2=r2) {
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
  W.prop <- solve(t(x.prop) %*% (x.prop * psi.ddot * pinv.prop))
  PI.mMSE <- sqrt((abs(Y - psi.dot-1e-6)+abs(Y - psi.dot+1e-6))^2/4 * rowSums((X%*%W.prop)^2))
  PI.mMSE <- PI.mMSE / sum(PI.mMSE)
  idx.mMSE <- sample(1:n, r2, T, PI.mMSE)
  x.mMSE <- X[c(idx.mMSE, idx.prop),]
  y.mMSE <- Y[c(idx.mMSE, idx.prop)]
  pinv.mMSE <- c(1 / PI.mMSE[idx.mMSE], pinv.prop)
  fit.mMSE <- getMLE(x=x.mMSE, y=y.mMSE, w=pinv.mMSE)
  beta.mMSE <- fit.mMSE$par
  msg <- fit.mMSE$message
  return(list(simp=beta.mMSE, msg=msg))    
}


AlgTwoStpmvc <- function(r1=r1, r2=r2) {
  idx.prop <- sample(1:n, r1, T)
  x.prop <- X[idx.prop,]
  y.prop <- Y[idx.prop]
  pinv.prop <- rep(n,r1)
  fit.prop <- getMLE(x=x.prop, y=y.prop, w=pinv.prop)
  beta.prop <- fit.prop$par
  if (is.na(beta.prop[1]))
    return(list(opt=NA, msg="first stage not converge"))
  psi.dot  <-  ( exp(X %*% beta.prop))
  
  
  PI.mVc <- sqrt((abs(Y - psi.dot-1e-6)+abs(Y - psi.dot+1e-6))^2/4 * rowSums(X^2))
  PI.mVc <- PI.mVc / sum(PI.mVc)
  idx.mVc <- sample(1:n, r2, T, PI.mVc)
  x.mVc <- X[c(idx.mVc, idx.prop),]
  y.mVc <- Y[c(idx.mVc, idx.prop)]
  pinv.mVc <- c(1 / PI.mVc[idx.mVc], pinv.prop)
  fit.mVc <- getMLE(x=x.mVc, y=y.mVc, w=pinv.mVc)
  beta.mVc <- fit.mVc$par
  msg <- fit.mVc$message
  return(list(simp=beta.mVc, msg=msg))    
}


leveragescore <- function(X){
  p<-d;
  if(d <= 10){
    temp <- svd(X)
    PI.Lva <- apply(temp$u^2,1,sum)
    p.slev <- PI.Lva/sum(PI.Lva)
  }
  else{
  r1<-p;r2<-20;
  SX<-matrix(0,r1,p);
  rn<-rmultinom(1,n,prob=rep(1,r1));
  for(i in 1:r1){
    index.pm<-sample(1:n,rn[i],F);
    Dd<-2*rbinom(rn[i],1,0.5)-1;
    SX[i,]<-Dd%*%(matrix(X[index.pm,],,p));
  }
  sv.cw<-svd(SX);
  V<-sv.cw$v;D<-sv.cw$d;
  Rinv<-V%*%diag(D);
  SRinv<-matrix(0,r2,p);
  rn<-rmultinom(1,p,prob=rep(1,r2));
  Rinv<-t(Rinv);
  for(i in 1:r2){
    index.pm<-sample(1:p,rn[i],F);
    Dd<-2*rbinom(rn[i],1,0.5)-1;
    SRinv[i,]<-Dd%*%(matrix(Rinv[index.pm,],,p));
  }
  B<-X%*%t(SRinv);
  p.slev<-rep(0,n);
  for(i in 1:n) {
    p.slev[i] <- sqrt(B[i,]%*%B[i,]);
  }
  p.slev <- p.slev/sum(p.slev)
  }
  return(p.slev)
}



AlgTwoStpLv <- function(r1=r1, r2=r2) {    
  ## leverage score
  PI.Lv <- leveragescore(X)
  idx.Lv <- sample(1:n, r1+r2, T, PI.Lv)
  x.Lv <- X[c(idx.Lv),]
  y.Lv <- Y[c(idx.Lv)]
  pinv.Lv <- c(1 / PI.Lv[idx.Lv])
  fit.Lv <- getMLE(x=x.Lv, y=y.Lv, w=pinv.Lv)
  beta.Lv <- fit.Lv$par
  msg <- fit.Lv$message
  return(list(simp=beta.Lv, msg=msg))    
}


AlgTwoStpLva <- function(r1=r1, r2=r2) {
  idx.prop <- sample(1:n, r1, T)
  x.prop <- X[idx.prop,]
  y.prop <- Y[idx.prop]
  pinv.prop <- rep(n,r1)
  fit.prop <- getMLE(x=x.prop, y=y.prop, w=pinv.prop)
  beta.prop <- fit.prop$par
  if (is.na(beta.prop[1]))
    return(list(opt=NA, msg="first stage not converge"))
  psi.dot  <-  ( exp(X %*% beta.prop))    
  
  ## leverage score:adjust
  lambdahat  <- sqrt(exp(X %*% beta.prop))
  temp <- sweep(X, 1, lambdahat[,1],FUN = '*')
  #temp <- svd(temp)
  #PI.Lva <- apply(temp$u^2,1,sum)
  #PI.Lva <- PI.Lva/sum(PI.Lva)
  PI.Lva <- leveragescore(temp)
  idx.Lva <- sample(1:n, r2, T, PI.Lva)
  x.Lva <- X[c(idx.Lva, idx.prop),]
  y.Lva <- Y[c(idx.Lva, idx.prop)]
  pinv.Lva <- c(1 / PI.Lva[idx.Lva], pinv.prop)
  fit.Lva <- getMLE(x=x.Lva, y=y.Lva, w=pinv.Lva)
  beta.Lva <- fit.Lva$par
  msg <- fit.Lva$message
  return(list(simp=beta.Lva, msg=msg))        
}

### test code ###

beta0  <- c(rep(1/2, 10),rep(0.1,20),rep(-0.2,10),rep(-0.1,20),rep(-0.5,20))
d <- length(beta0)
for(d in c(7,40,80)){
for(n in c(1e5,5e5)){
set.seed(1234)
beta1 <- beta0[1:d]
X  <- matrix(runif(n*d, 0, 1),n,d)
X[,2] <- X[,1]+runif(n,0,1)
X[,6] <- runif(n,-1,1)
X[,7] <- runif(n,-1,1)
lambda  <- exp(X %*% beta1)
Y  <- rpois(n, lambda)

r1 <- ifelse(d==7,800,ifelse(d==40,800,800))

for(r in c(1000,1500,2000,2500)){

library(microbenchmark)
tm <- microbenchmark(getMLE(X, Y, 1),
                     AlgTwoStpsimp(r1+r,0),
                     AlgTwoStpmse(r1, r),
                     AlgTwoStpmvc(r1, r),
                     AlgTwoStpLv(r1, r),
                     AlgTwoStpLva(r1, r),times=50L)
cat("d=",d,"n=",n,"r=",r,"\n")
print(tm)
write.csv(tm,paste0("D:/usetime_", as.character(d), "_",as.character(n),"_",as.character(r),"e1.csv"))
}
}
}