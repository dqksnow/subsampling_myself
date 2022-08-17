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





## two-step-alg

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
    

    
    
    
    opt <- cbind(fit.mMSE$par, fit.mVc$par)
    msg <- c(fit.prop$message, fit.mMSE$message,  fit.mVc$message)
    
    return(list(opt=opt, msg=msg))
  }
}



### test code ###
set.seed(1234)
n<-1e5;
beta0  <- c(rep(1/2, 7))
d <- length(beta0)
X  <- matrix(runif(n*d, 0, 1),n,d)
X[,2] <- X[,1]+runif(n,0,1)
X[,6] <- runif(n,-1,1)
X[,7] <- runif(n,-1,1)
lambda  <- exp(X %*% beta0)
Y  <- rpois(n, lambda)
print(mean(Y))

lambda  <- exp(X %*% beta0)
Y  <- rpois(n, lambda)
print(mean(Y))

fit.full <- getMLE(X, Y, 1)

rpt <- 1000
r <-  1200

r.ss <- round(seq(100, r-50, length=20))

lrs <- length(r.ss)
Beta.simp <- matrix(NA, d, rpt*lrs)
nm <- 2
Beta.opt <- matrix(NA, d*nm, rpt*lrs)
itr <- 0

for (r1 in r.ss) {
  r2 <- r - r1
  cat(r1, r2, r, " ");
  set.seed(10)
  for (rr in 1:rpt) {
    if (rr%/%100 == rr/100) cat(rr, "-")
    itr <- itr+1
    fit.alg.simp <- AlgTwoStp(r1+r2, 0)
    Beta.simp[,itr] <- fit.alg.simp$simp
    if (any(fit.alg.simp$msg != "Successful convergence")) {
      warning(paste(rr, "not converge", fit.alg$msg.simp))
      Beta.simp[,itr] <- NA
    }
    fit.alg <- AlgTwoStp(r1, r2)
    if (any(fit.alg$msg != "Successful convergence")) {
      warning(paste(rr, "not converge", fit.alg$msg))
      next
    }
    Beta.opt[,itr] <- c(fit.alg$opt)
  }
  cat("\n")
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


plot(r.ss/r, rec[,1], xlab=expression("r"[0]/("r"[0]+r)), ylab="MSE", pch="1",
     ylim=c(min(rec), max(rec)), type="o")
pp <- dim(rec)[2]
for(i in 2:pp)
  lines(r.ss/r, rec[,i], type="o", pch=paste(i), lty=i, col=i)
legend("bottomright", lty=1:pp, pch=paste(1:pp), col=1:pp,
       legend=c("uniform", "MV", "MV-LD", "MVC", "MVC-LD"))


plot(r.ss/r, rec[,1], xlab=expression("r"[0]/("r"[0]+r)), ylab="MSE", pch="1",lwd=2,
     ylim=c(min(rec), max(rec)), type="o")
pp <- dim(rec)[2]
for(i in 2:pp)
  lines(r.ss/r, rec[,i], type="o", pch=paste(i), lty=i, col=i,lwd=2)
legend("bottomright", lty=1:pp, pch=paste(1:pp), col=1:pp,
       legend=c("UNIF", "mMV ","mMVc"),cex=1.5,ncol=1)


