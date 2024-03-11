rm(list=ls())
# Newton method for poisson regression
getMLE <- function(x, y, w, size) {
  beta <- rep(0, d)
  loop  <- 1
  Loop  <- 100
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

AlgTwoStp <- function(r1=r1, r2=r2, size) {
  if (r2 == 0) {
    idx.simp <- sample(1:n, r1, T)
    x.simp <- X[idx.simp,]
    y.simp <- Y[idx.simp]
    pinv.simp <- n
    fit.simp <- getMLE(x=x.simp, y=y.simp, w=pinv.simp, size = size)
    beta.simp <- fit.simp$par
    msg <- fit.simp$message
    psi.dot <- exp(x.simp%*%beta.simp)
    hatpsidot <- exp(x.simp%*%beta.simp)* size * (size+Y[idx.simp])/((size+psi.dot)^2)
    J_Xsimp <- 1/(n*(r1+r2))*t(x.simp)%*%diag(c(hatpsidot*pinv.simp))%*%x.simp
    V_simp <- 1/(n^2*(r1+r2)^2)*t(x.simp)%*%diag(c((y.simp-psi.dot)^2*(size)^2/((size+psi.dot)^2)*(pinv.simp^2)))%*%x.simp
    sigma_simp <- diag(solve(J_Xsimp)%*%V_simp%*%solve(J_Xsimp))[2]
    return(list(simp=beta.simp, msg=msg,sigma_simp=sigma_simp))
  }
  if (r2 != 0) {
    idx.prop <- sample(1:n, r1, T)
    x.prop <- X[idx.prop,]
    y.prop <- Y[idx.prop]
    pinv.prop <- rep(n,r1)
    fit.prop <- getMLE(x=x.prop, y=y.prop, w=pinv.prop, size = size)
    beta.prop <- fit.prop$par
    if (is.na(beta.prop[1]))
      return(list(opt=NA, msg="first stage not converge"))
    psi.dot  <-  ( exp(X %*% beta.prop))
    u <- size/(size+psi.dot)
    
    ## mVc
    PI.mVc <- sqrt((abs(Y - psi.dot-1e-6)+abs(Y - psi.dot+1e-6))^2/4 *abs(u^2)* rowSums(X^2))
    PI.mVc <- PI.mVc / sum(PI.mVc)
    idx.mVc <- sample(1:n, r2, T, PI.mVc)
    x.mVc <- X[c(idx.mVc, idx.prop),]
    y.mVc <- Y[c(idx.mVc, idx.prop)]
    pinv.mVc <- c(1 / PI.mVc[idx.mVc], pinv.prop)
    fit.mVc <- getMLE(x=x.mVc, y=y.mVc, w=pinv.mVc, size = size)
    psi.dot1 <- exp(x.mVc%*%fit.mVc$par)
    hatpsidot <- psi.dot1 * size * (size+y.mVc)/((size+psi.dot1)^2)
    J_XmVc <- 1/(n*(r1+r2))*t(x.mVc)%*%diag(c(hatpsidot*pinv.mVc))%*%x.mVc
    V_mVc <- 1/(n^2*(r1+r2)^2)*t(x.mVc)%*%diag(c((y.mVc-psi.dot1)^2*(size)^2/((size+psi.dot1)^2)*(pinv.mVc^2)))%*%x.mVc
    sigma_mVc <- diag(solve(J_XmVc)%*%V_mVc%*%solve(J_XmVc))[2]
    
    
    ## mMSE
    psi.ddot <- psi.dot[idx.prop]* size * (size+Y[idx.prop])/((size+psi.dot[idx.prop])^2)
    W.prop <- solve(t(x.prop) %*% (x.prop * psi.ddot * pinv.prop))
    PI.mMSE <- sqrt((abs(Y - psi.dot-1e-6)+abs(Y - psi.dot+1e-6))^2/4 * abs(u^2) * rowSums((X%*%W.prop)^2))
    PI.mMSE <- PI.mMSE / sum(PI.mMSE)
    idx.mMSE <- sample(1:n, r2, T, PI.mMSE)
    x.mMSE <- X[c(idx.mMSE, idx.prop),]
    y.mMSE <- Y[c(idx.mMSE, idx.prop)]
    pinv.mMSE <- c(1 / PI.mMSE[idx.mMSE], pinv.prop)
    fit.mMSE <- getMLE(x=x.mMSE, y=y.mMSE, w=pinv.mMSE, size = size)
    psi.dot1 <- exp(x.mMSE%*%fit.mMSE$par)
    hatpsidot <- psi.dot1 * size * (size+y.mMSE)/((size+psi.dot1)^2)
    J_XmMSE <- 1/(n*(r1+r2))*t(x.mMSE)%*%diag(c(hatpsidot*pinv.mMSE))%*%x.mMSE
    V_mMSE <- 1/(n^2*(r1+r2)^2)*t(x.mMSE)%*%diag(c((y.mMSE-psi.dot1)^2*(size)^2/((size+psi.dot1)^2)*(pinv.mMSE^2)))%*%x.mMSE
    sigma_mMSE <- diag(solve(J_XmMSE)%*%V_mMSE%*%solve(J_XmMSE))[2]
    
    

    
    
    opt <- cbind(fit.mMSE$par,  fit.mVc$par)
    msg <- c( fit.mMSE$message,  fit.mVc$message)
    sigma_est <- c(sigma_mMSE,sigma_mVc)
    return(list(opt=opt, msg=msg, sigma_est=sigma_est))
  }
}

### test code ###
set.seed(1234)
n<-1e5;
beta0  <- c(rep(1/2, 7))
d <- length(beta0)
X  <- matrix(runif(n*d, 0, 1),n,d)
#X[,2] <- X[,1]+runif(n,0,0.1)
# X[,6] <- runif(n,-1,1)
# X[,7] <- runif(n,-1,1)


lambda  <- exp(X %*% beta0)
Y  <- rpois(n, lambda)
print(mean(Y))

fit.full <- getMLE(X, Y, 1,2)

rpt <- 1000
r1 <- 200
r <- 300
nm <- 2
Beta.opt <- matrix(NA, d*nm, rpt)
Sigma.est <- matrix(NA,nm,rpt)
coverrate <- rep(NA,nm+1)
Beta.simp <- matrix(NA, d, rpt)
Sigma.simp <- matrix(NA,1,rpt)
conflen <- matrix(NA, nm+1, rpt)

set.seed(0)
for(i in 1:rpt){
  if (i%/%100 == i/100) cat(i, "-")
  fit.alg <- AlgTwoStp(r1, r,2)
  Beta.opt[,i] <- c(fit.alg$opt)
  Sigma.est[,i] <- c(fit.alg$sigma_est)
}

set.seed(0)
for(i in 1:rpt){
  if (i%/%100 == i/100) cat(i, "-")
  fit.alg <- AlgTwoStp(r1+r,0,2)
  Beta.simp[,i] <- c(fit.alg$simp)
  Sigma.simp[,i] <- c(fit.alg$sigma_simp)
}

for(j in 1:nm){
  # plot(NULL
  #      ,xlim = c(0,100)
  #      ,ylim = c(0.3,0.7)
  #      ,xaxt = 'n'
  #      ,ylab = '95% Confidence level'
  #      ,xlab = 'index'
  #      #,main = "Confidence Intervals of 1000 Subsamples"
  # )
  # 
  # abline(h = fit.full$par[2], col = 'black')
  # axis(1,at=c(1,10,20,30,40,50,60,70,80,90,100))
  # 
  # for (i in 1:100){
  #   interval1 <- Beta.opt[7*j-5,i]-1.96*sqrt(Sigma.est[j,i])
  #   interval2 <- Beta.opt[7*j-5,i]+1.96*sqrt(Sigma.est[j,i])
  #   if(0.5>=interval1 & 0.5<=interval2){
  #     lines(c(i,i),c(interval1,interval2), lwd=2,col='black')
  #   }
  #   else{
  #     lines(c(i,i),c(interval1,interval2), lwd=2,col='red' )
  #   } 
  # }
  
  
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





# plot(NULL
#      ,xlim = c(0,100)
#      ,ylim = c(0.3,0.7)
#      ,xaxt = 'n'
#      ,ylab = '95% Confidence level'
#      ,xlab = 'index'
#      #,main = "Confidence Intervals of 1000 Subsamples"
# )
# 
# abline(h = 0.5, col = 'black')
# axis(1,at=c(1,10,20,30,40,50,60,70,80,90,100))
# 
# for (i in 1:100){
#   interval1 <- Beta.simp[2,i]-1.96*sqrt(Sigma.simp[1,i])
#   interval2 <- Beta.simp[2,i]+1.96*sqrt(Sigma.simp[1,i])
#   if(0.5>=interval1 & fit.full$par[2]<=interval2){
#     lines(c(i,i),c(interval1,interval2), lwd=2,col='black')
#   }
#   else{
#     lines(c(i,i),c(interval1,interval2), lwd=2,col='red' )
#   } 
# }

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

coverrate
rowMeans(conflen)