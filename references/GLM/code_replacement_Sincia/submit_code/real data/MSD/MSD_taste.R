rm(list=ls())
# Newton method for poisson regression
getMLE <- function(x, y, w) {
  beta <- c(-1.36834,0.02895,-0.02203,-0.03506,5.30977,4.45400,1.72954)
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
      beta <-  NA
      print(loop)
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


### data ###
df <- read.csv("D:/讨论班/report4/data/MSDtaste1.csv")
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

fit.full <- getMLE(X, Y, 1)
##seems right


## r1 <- r1+r; r2=0

AlgTwoStp <- function(r1=r1, r2=r2) {
  if (r2 == 0) {
    idx.simp <- sample(1:n, r1, T)
    x.simp <- X[idx.simp,]
    y.simp <- Y[idx.simp]
    pinv.simp <- n
    fit.simp <- getMLE(x=x.simp, y=y.simp, w=pinv.simp)
    beta.simp <- fit.simp$par
    msg <- fit.simp$message
    return(list(simp=beta.simp, msg=msg))
  }
  if (r2 != 0) {
    idx.prop <- sample(1:n, r1, T)
    x.prop <- X[idx.prop,]
    y.prop <- Y[idx.prop]
    pinv.prop <- rep(n,r1)
    fit.prop <- getMLE(x=x.prop, y=y.prop, w=pinv.prop)
    beta.prop <- fit.prop$par
    if (is.na(beta.prop[1]))
      return(list(opt=NA, msg="first stage not converge"))
    psi.dot  <-  ( exp(X %*% beta.prop))
    
    ## mVc
    PI.mVc <- sqrt((abs(Y - psi.dot-1e-6)+abs(Y - psi.dot+1e-6))^2/4 * rowSums(X^2))
    PI.mVc <- PI.mVc / sum(PI.mVc)
    idx.mVc <- sample(1:n, r2, T, PI.mVc)
    x.mVc <- X[c(idx.mVc, idx.prop),]
    y.mVc <- Y[c(idx.mVc, idx.prop)]
    pinv.mVc <- c(1 / PI.mVc[idx.mVc], pinv.prop)
    fit.mVc <- getMLE(x=x.mVc, y=y.mVc, w=pinv.mVc)
    
    
    
    ## mMSE
    psi.ddot <- psi.dot[idx.prop]
    W.prop <- solve(t(x.prop) %*% (x.prop * psi.ddot * pinv.prop))
    PI.mMSE <- sqrt((abs(Y - psi.dot-1e-6)+abs(Y - psi.dot+1e-6))^2/4 * rowSums((X%*%W.prop)^2))
    PI.mMSE <- PI.mMSE / sum(PI.mMSE)
    idx.mMSE <- sample(1:n, r2, T, PI.mMSE)
    x.mMSE <- X[c(idx.mMSE, idx.prop),]
    y.mMSE <- Y[c(idx.mMSE, idx.prop)]
    pinv.mMSE <- c(1 / PI.mMSE[idx.mMSE], pinv.prop)
    fit.mMSE <- getMLE(x=x.mMSE, y=y.mMSE, w=pinv.mMSE)
    
    
    
    ## leverage score
    temp <- svd(X)
    PI.Lv <- apply(temp$u^2,1,sum)
    PI.Lv <- PI.Lv/sum(PI.Lv)
    idx.Lv <- sample(1:n, r1+r2, T, PI.Lv)
    x.Lv <- X[c(idx.Lv),]
    y.Lv <- Y[c(idx.Lv)]
    pinv.Lv <- c(1 / PI.Lv[idx.Lv])
    fit.Lv <- getMLE(x=x.Lv, y=y.Lv, w=pinv.Lv)
    
    ## leverage score:adjust
    lambdahat  <- sqrt(exp(X %*% beta.prop))
    temp <- sweep(X, 1, lambdahat[,1],FUN = '*')
    temp <- svd(temp)
    PI.Lva <- apply(temp$u^2,1,sum)
    PI.Lva <- PI.Lva/sum(PI.Lva)
    idx.Lva <- sample(1:n, r2, T, PI.Lva)
    x.Lva <- X[c(idx.Lva, idx.prop),]
    y.Lva <- Y[c(idx.Lva, idx.prop)]
    pinv.Lva <- c(1 / PI.Lva[idx.Lva], pinv.prop)
    fit.Lva <- getMLE(x=x.Lva, y=y.Lva, w=pinv.Lva)
    
    
    opt <- cbind(fit.mMSE$par, fit.mVc$par, fit.Lv$par, fit.Lva$par)
    msg <- c( fit.mMSE$message,  fit.mVc$message, fit.Lv$message, fit.Lva$message)
    
    return(list(opt=opt, msg=msg))
  }
}

### test code ###
r1 <- 400
r <- 2000

set.seed(0)
rpt <- 1000
r.ss <- c(300, 500, 700,1000,1200,1400)*2## [5:6]
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
    fit.alg.simp <- AlgTwoStp(r1+r, 0)
    Beta.simp[,itr] <- fit.alg.simp$simp
    if (any(fit.alg.simp$msg != "Successful convergence")) {
      warning(paste(rr, "not converge", fit.alg.simp$msg))
      Beta.simp[,itr] <- NA
    }
    fit.alg <- AlgTwoStp(r1, r)
    if (any(fit.alg$msg != "Successful convergence")) {
      warning(paste(rr, "not converge:", fit.alg$msg))
      next
    }
    Beta.opt[,itr] <- c(fit.alg$opt)
  }
  cat("\n\n")
}

beta0 <- fit.full$par



# pred <- function(beta,x,y){
#   l2.loss <- (y-exp((x)%*%(beta)))^2
#   mean(l2.loss)
# }
# pred.simp <- apply(Beta.simp,2,pred, x=as.matrix(X),y=Y)
# 
# pred.optmse <- apply(Beta.opt[1:8,],2,pred, x=as.matrix(X),y=Y)
# pred.optmvc <- apply(Beta.opt[9:16,],2,pred, x=as.matrix(X),y=Y)
# pred.optlev <- apply(Beta.opt[17:24,],2,pred, x=as.matrix(X),y=Y)
# pred.optleva <- apply(Beta.opt[25:32,],2,pred, x=as.matrix(X),y=Y)
# 
# 
# 
# pred.opt <- rbind(pred.optmse,pred.optmvc,pred.optlev,pred.optleva)
# 
# mspe.simp <- mspe.opt <- c()
# for (i in 1:lrs) {
#   loc <- ((i-1)*rpt+1):(i*rpt)
#   mspe.simp <- c(mspe.simp, (mean(pred.simp[loc], na.rm=TRUE)))
#   mspe.opt <- rbind(mspe.opt, apply(pred.opt[,loc],1,mean, na.rm=TRUE))
# }
# par(mar=c(4,5,4,3))
# pp <- dim(repc)[2]
# repc <- cbind(mspe.simp, mspe.opt)
# plot(r.ss, log10(repc[,1]), xlab="r", ylab=expression(log[10](MSPE)), pch="1",lwd=2,
#      ylim=c(min(log10(repc)), max(log10(repc))), type="o")
# for(i in 2:pp)
#   lines(r.ss, log10(repc[,i]), type="o", pch=paste(i), lty=i, col=i,lwd=2)
# legend("topright", lty=1:pp, pch=paste(1:pp), col=1:pp,
#        legend=c("UNIF", "mMV", "mMVc", "Lev","Lev-A"),cex=1.5)
# 
# 
l2.simp <- (Beta.simp - beta0)^2
l2.opt  <- (Beta.opt  - beta0)^2

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
       legend=c("UNIF", "mMV", "mMVc", "Lev","Lev-A"),cex=1.5)