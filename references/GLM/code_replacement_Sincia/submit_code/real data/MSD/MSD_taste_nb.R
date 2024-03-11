rm(list=ls())
# Newton method for poisson regression
getMLE <- function(x, y, w, size) {
  beta <- c(0.92771138,0.04061366,-0.01477433,-0.05080527,3.98471833,2.10953425,1.24296298)
  loop  <- 1
  Loop  <- 400
  msg <- "NA"
  while (loop <= Loop) {
    pr <- c(exp(x %*% beta))
    u <- size/(size+pr)
    H <- t(x) %*% (pr * u * w * x)
    S <- colSums((y - pr) * u * w * x)
    tryCatch(
      {shs <- NA
      shs <- solve(H, S) }, 
      error=function(e){
        cat("\n ERROR :", loop, conditionMessage(e), "\n")})
    if (is.na(shs[1])) {
      msg <- "Not converge"
      beta <- loop <- NA
      #break
      return(list(par=beta, message=msg, iter=loop))
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
  return(list(par=beta, message=msg, iter=loop))
}





## r1 <- r1+r; r2=0
library(MASS)
AlgTwoStp <- function(r1=r1, r2=r2,size=size) {
  if (r2 == 0) {
    idx.simp <- sample(1:n, r1, T)
    x.simp <- X[idx.simp,]
    y.simp <- Y[idx.simp]
    pinv.simp <- n
    # fit.pilot <- glm.nb(y.simp[1:200]~-1+x.simp[1:200])
    # size <- fit.pilot$theta
    fit.simp <- getMLE(x=x.simp, y=y.simp, w=pinv.simp, size = size)
    beta.simp <- fit.simp$par
    msg <- fit.simp$message
    return(list(simp=beta.simp, msg=msg))
  }
  if (r2 != 0) {
    idx.prop <- sample(1:n, r1, T)
    x.prop <- X[idx.prop,]
    y.prop <- Y[idx.prop]
    pinv.prop <- rep(n,r1)
    # fit.pilot <- glm.nb(y.prop~-1+x.prop)
    # size <- fit.pilot$theta
    fit.prop <- getMLE(x=x.prop, y=y.prop, w=pinv.prop, size = size)
    beta.prop <- fit.prop$par
    if (is.na(beta.prop[1]))
      return(list(opt=NA, msg="first stage not converge"))
    psi.dot  <-  ( exp(X %*% beta.prop))
    u <- size/(size+psi.dot)
    
    ## mVc
    PI.mVc <- sqrt((abs(Y - psi.dot-1e-6)+abs(Y - psi.dot+1e-6))^2/4 *abs(u^2)* rowSums(X^2))
    PI.mVc <- PI.mVc / sum(PI.mVc)
    if(!any(is.na(PI.mVc))){
      idx.mVc <- sample(1:n, r2, T, PI.mVc)
      x.mVc <- X[c(idx.mVc, idx.prop),]
      y.mVc <- Y[c(idx.mVc, idx.prop)]
      pinv.mVc <- c(1 / PI.mVc[idx.mVc], pinv.prop)
      fit.mVc <- getMLE(x=x.mVc, y=y.mVc, w=pinv.mVc, size = size)
    }
    else{
      fit.mVc <- list(par=rep(NA,d),message="NA in Prob")
    }

    
    
    
    ## mMSE
    psi.ddot <- psi.dot[idx.prop]* size * (size+Y[idx.prop])/((size+psi.dot[idx.prop])^2)
    tryCatch(
      {W.prop <- NA
       W.prop <- solve(t(x.prop) %*% (x.prop * psi.ddot * pinv.prop)) }, 
      error=function(e){
        cat("\n ERROR :",  conditionMessage(e), "\n")})
    if (any(is.na(W.prop))) {
      fit.mMSE <- list(par=rep(NA,d),message="NA in J_X
                       li")
    }
    else{
    PI.mMSE <- sqrt((abs(Y - psi.dot-1e-6)+abs(Y - psi.dot+1e-6))^2/4 * abs(u^2) * rowSums((X%*%W.prop)^2))
    PI.mMSE <- PI.mMSE / sum(PI.mMSE)
    if(!any(is.na(PI.mMSE))){
      idx.mMSE <- sample(1:n, r2, T, PI.mMSE)
      x.mMSE <- X[c(idx.mMSE, idx.prop),]
      y.mMSE <- Y[c(idx.mMSE, idx.prop)]
      pinv.mMSE <- c(1 / PI.mMSE[idx.mMSE], pinv.prop)
      fit.mMSE <- getMLE(x=x.mMSE, y=y.mMSE, w=pinv.mMSE, size = size)
    }
    else{
      fit.mMSE <- list(par=rep(NA,d),message="NA in Prob")
    }
    }

    
    
    
    ## leverage score
    temp <- svd(X)
    PI.Lv <- apply(temp$u^2,1,sum)
    PI.Lv <- PI.Lv/sum(PI.Lv)
    idx.Lv <- sample(1:n, r1+r2, T, PI.Lv)
    x.Lv <- X[c(idx.Lv),]
    y.Lv <- Y[c(idx.Lv)]
    pinv.Lv <- c(1 / PI.Lv[idx.Lv])
    fit.Lv <- getMLE(x=x.Lv, y=y.Lv, w=pinv.Lv, size = size)
    
    ## leverage score:adjust
    lambdahat  <- sqrt(u * psi.dot)
    temp <- sweep(X, 1, lambdahat[,1],FUN = '*')
    tryCatch(
      {temp1 <- NA
      temp1 <- svd(temp)},
      error=function(e){
        cat("\n ERROR :",  conditionMessage(e), "\n")})
    if (any(is.na(temp1))) {
      fit.Lva <- list(par=rep(NA,d),message="NA in SVD")
    }
    else{
      PI.Lva <- apply(temp1$u^2,1,sum)
      PI.Lva <- PI.Lva/sum(PI.Lva) 
      if(!any(is.na(PI.Lva))){
      idx.Lva <- sample(1:n, r2, T, PI.Lva)
      x.Lva <- X[c(idx.Lva, idx.prop),]
      y.Lva <- Y[c(idx.Lva, idx.prop)]
      pinv.Lva <- c(1 / PI.Lva[idx.Lva], pinv.prop)
      fit.Lva <- getMLE(x=x.Lva, y=y.Lva, w=pinv.Lva, size = size)
    }
    else{
      fit.Lva <- list(par=rep(NA,d),message="NA in Prob")
    }
    }
    

    
    
    opt <- cbind(fit.mMSE$par,  fit.mVc$par,fit.Lv$par,fit.Lva$par)
    msg <- c( fit.mMSE$message,  fit.mVc$message,fit.Lv$message,fit.Lva$message)
    
    return(list(opt=opt, msg=msg))
  }
}

### test code ###
df <- read.csv("D:/ÌÖÂÛ°à/report4/data/MSDtaste1.csv")
df1 <- df[complete.cases(df),c(3,8:18)]
X <- df1[,c(2:3,7:9)]
X <- scale(X)
data <- cbind(X,df1[,c(4:6,10:12,1)])
data[,9] <- ifelse(df1[,10]>2000,1,0)

data <- data[data[,8]==1,-c(2,5,8,9,11)]
head(data)
X <- cbind(1,data[,1:6])
X <- as.matrix(X)
Y <- data[,7]
n<- nrow(data)
d <- 7

fit.full <- getMLE(X, Y, 1, size = 1.4)


### select theta
r1<-400
set.seed(0)
idx.simp <- sample(1:n, r1, T)
x.simp <- X[idx.simp,]
y.simp <- Y[idx.simp]
glm.nb(Y[Y<247]~-1+X[Y<247,])

winsorVar<-function(x,probs=c(0.25,0.75)) {
  xq<-quantile(x,probs=probs)
  x[x < xq[1]]<-xq[1]
  x[x > xq[2]]<-xq[2]
  return(var(x,na.rm = T))
}
winsorMean<-function(x,probs=c(0.25,0.75)) {
  xq<-quantile(x,probs=probs)
  x[x < xq[1]]<-xq[1]
  x[x > xq[2]]<-xq[2]
  return(mean(x,na.rm = T))
}
winsorMean(Y)^2/(winsorVar(Y)-winsorMean(Y))

### test code ###
r1 <- 400
r <- 1000

set.seed(0)
rpt <- 1000
r.ss <- c( 300, 500, 700, 1000, 1200)*2## [5:6]
lrs <- length(r.ss)
Beta.simp <- matrix(NA, d, rpt*lrs)
nm <- 4
Beta.opt <- matrix(NA, d*nm, rpt*lrs)
itr <- 0

for (r in r.ss) {
  cat(r, " ");
  set.seed(0)
  for (rr in 1:rpt) {
    if (rr%/%100 == rr/100) cat(rr, "-")
    itr <- itr+1
    fit.alg.simp <- AlgTwoStp(r1+r, 0, size = 0.8)
    Beta.simp[,itr] <- fit.alg.simp$simp
    if (any(fit.alg.simp$msg != "Successful convergence")) {
      warning(paste(rr, "not converge", fit.alg.simp$msg))
      Beta.simp[,itr] <- NA
    }
    fit.alg <- AlgTwoStp(r1, r, size = 0.8)
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
plot(r.ss, rec[,1], xlab="r", ylab="MSE", pch="1",lwd=2,
     ylim=c(min(rec), max(rec)), type="o")
pp <- dim(rec)[2]
for(i in 2:pp)
  lines(r.ss, rec[,i], type="o", pch=paste(i), lty=i, col=i,lwd=2)
legend("topright", lty=1:pp, pch=paste(1:pp), col=1:pp,
       legend=c("UNIF", "mV","mVc", "Lev","Lev-A"),cex=1.5,ncol=1)

write.csv(rec,"D:/report4/caserd-nb-theta-08.csv")