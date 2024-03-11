rm(list=ls())
# Newton method for logistic regression
getMLE <- function(x, y, w) {
  beta <- c(1.375e-01,4.013e-04,2.945e-05,2.155e-05,-6.708e-04)
  loop  <- 1
  Loop  <- 100
  msg <- "NA"
  while (loop <= Loop) {
    pr <- c(1/(x %*% beta))
    Jx <- c(1/((x %*% beta)^2))
    H <- t(x) %*% ((Jx) * w * x) 
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
    tlr  <- sum((beta.new - beta)^2)/sum((beta)^2)
    beta  <- beta.new
    if(tlr < 0.001) {
      msg <- "Successful convergence"
      break
    }
    if (loop == Loop)
      warning("Maximum iteration reached")
    loop  <- loop + 1
  }
  list(par=beta, message=msg, iter=loop, Infor = H)
}

#fit1 <- getMLE(X1,Y1,1)


## r1 <- r1+r; r2=0

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
  psi.dot  <-   1/(X %*% beta)
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
  psi.dot  <-   1/(X %*% beta)
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


AlgTwoStp <- function(X,Y,r1=r1, r2=r2,rho=0.01,K=1,method="simp",a1,b1) {
  if (method == "simp") {
    n <- nrow(X)
    idx <- 1:n
    idx.prop <- sample(1:n, r1, T)
    x.prop <- X[idx.prop,]
    y.prop <- Y[idx.prop]
    pinv.prop <- rep(n,r1)
    pilot <- getMLE(x=x.prop, y=y.prop, w=pinv.prop)
    #beta.prop <- fit.prop$par
    # pilot <- AlgTwoStpsimp(X,Y,r1)
    #n <- nrow(X)
    # <- floor(n/K)

    result <- foreach(a=a1,b=b1)%do%AlgTwoStpsimp(a,b,r2/K)
    H <- pilot$Infor
    B <- pilot$Infor%*%pilot$par
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
  if (method == "mse") {
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
    psi.dot  <-   c(1/(X[idx.prop,] %*% beta.prop))
    
    ## mMSE
    p.prop <- psi.dot
    psi.ddot <- p.prop^2 
    W.prop <- solve(t(x.prop) %*% (x.prop * psi.ddot * pinv.prop))
    Psi <- sum(abs(Y[idx.prop] - p.prop) * rowSums((X[idx.prop,]%*%W.prop)^2) )/r1
    # m <- floor(n/K)
    # a <- b <- c()
    # for(i in 1:K){
    #   if(i==K){
    #     a <- c(a,list(X[((K-1)*m+1):nrow(X),]))
    #     b <- c(b,list(Y[((K-1)*m+1):nrow(X)]))
    #   }
    #   else{
    #     a <- c(a,list(X[((i-1)*m+1):(i*m),]))
    #     b <- c(b,list(Y[((i-1)*m+1):(i*m)]))
    #   }
    # }
    result <- foreach(a=a1,b=b1)%do%AlgTwoStpmse(a,b,r2/K,rho,beta.prop,W.prop,Psi)
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
    opt <- cbind(beta.mMSE)
    msg <- c(msg.mMSE)
    
    return(list(opt=opt, msg=msg))
  }
  if (method == "mvc") { 
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
    psi.dot  <-   c(1/(X[idx.prop,] %*% beta.prop))
    #p.prop <- psi.dot
    ## mVc
    #PI.mVc <- abs(Y - psi.dot) * rowSums(X^2)
    Psi <- sum(abs(Y[idx.prop] - p.prop) * rowSums((X[idx.prop,])^2))/r1
    # m <- floor(n/K)
    # a <- b <- c()
    # for(i in 1:K){
    #   if(i==K){
    #     a <- c(a,list(X[((K-1)*m+1):nrow(X),]))
    #     b <- c(b,list(Y[((K-1)*m+1):nrow(X)]))
    #   }
    #   else{
    #     a <- c(a,list(X[((i-1)*m+1):(i*m),]))
    #     b <- c(b,list(Y[((i-1)*m+1):(i*m)]))
    #   }
    # }
    result <- foreach(a=a1,b=b1)%do%AlgTwoStpmvc(a,b,r2/K,rho,beta.prop,Psi)
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
    
    opt <- cbind( beta.mVc)
    msg <- c(msg.mVc)
    
    return(list(opt=opt, msg=msg))
  }
}

### test code ###
library(data.table)
realdata_1 <- fread("D:/讨论班/report5/data/airlinedata/2008_clean.csv", header = T, sep = ',')
realdata_2 <- fread("D:/讨论班/report5/data/airlinedata/2007_clean.csv", header = T, sep = ',')
realdata_3 <- fread("D:/讨论班/report5/data/airlinedata/2006_clean.csv", header = T, sep = ',')
realdata_4 <- fread("D:/讨论班/report5/data/airlinedata/2005_clean.csv", header = T, sep = ',')
realdata_5 <- fread("D:/讨论班/report5/data/airlinedata/2004_clean.csv", header = T, sep = ',')
realdata_6 <- fread("D:/讨论班/report5/data/airlinedata/2003_clean.csv", header = T, sep = ',')
realdata_7 <- fread("D:/讨论班/report5/data/airlinedata/2002_clean.csv", header = T, sep = ',')
realdata_8 <- fread("D:/讨论班/report5/data/airlinedata/2001_clean.csv", header = T, sep = ',')
realdata_9 <- fread("D:/讨论班/report5/data/airlinedata/2000_clean.csv", header = T, sep = ',')
realdata_10 <- fread("D:/讨论班/report5/data/airlinedata/1999_clean.csv", header = T, sep = ',')
realdata_11 <- fread("D:/讨论班/report5/data/airlinedata/1998_clean.csv", header = T, sep = ',')
realdata_12 <- fread("D:/讨论班/report5/data/airlinedata/1997_clean.csv", header = T, sep = ',')
realdata_13 <- fread("D:/讨论班/report5/data/airlinedata/1996_clean.csv", header = T, sep = ',')
realdata_14 <- fread("D:/讨论班/report5/data/airlinedata/1995_clean.csv", header = T, sep = ',')
realdata_15 <- fread("D:/讨论班/report5/data/airlinedata/1994_clean.csv", header = T, sep = ',')
realdata_16 <- fread("D:/讨论班/report5/data/airlinedata/1993_clean.csv", header = T, sep = ',')
realdata_17 <- fread("D:/讨论班/report5/data/airlinedata/1992_clean.csv", header = T, sep = ',')
realdata_18 <- fread("D:/讨论班/report5/data/airlinedata/1991_clean.csv", header = T, sep = ',')
realdata_19 <- fread("D:/讨论班/report5/data/airlinedata/1990_clean.csv", header = T, sep = ',')
realdata_20 <- fread("D:/讨论班/report5/data/airlinedata/1989_clean.csv", header = T, sep = ',')
realdata_21 <- fread("D:/讨论班/report5/data/airlinedata/1988_clean.csv", header = T, sep = ',')
realdata_22 <- fread("D:/讨论班/report5/data/airlinedata/1987_clean.csv", header = T, sep = ',')


X <- rbind(realdata_1,realdata_2,realdata_3,realdata_4,realdata_5,realdata_6,realdata_7,realdata_8,
           realdata_9,realdata_10,realdata_11,realdata_12,realdata_13,realdata_14,realdata_15,realdata_16,
           realdata_17,realdata_18,realdata_19,realdata_20,realdata_21,realdata_22)

#rec1[,4] <- scale(rec1[,4])



Y <- (X$ArrDelay+1440)
X1 <- cbind(1,X[,3:6])
X1 <- as.matrix(X1)
X1[,2] <- scale(X1[,2])

Y1 <- log(Y)
# index <- sample(1:nrow(X1),5e6,replace = F)
# model1 <- glm(Y1[index]~-1+X1[index,],family = Gamma(link = "inverse"))
d <- dim(X1)[2]
fit.full <- getMLE(X1, Y1, 1)

r1 <- 800
k <- 5
r.ss <- c(500,800,1000,1200,1500,1700,2000,2200,2500)*2

set.seed(0)
rpt <- 500
#r.ss <- c( 300, 500, 700, 1000, 1200, 1400)## [5:6]
lrs <- length(r.ss)
Beta.simp <- matrix(NA, d, rpt*lrs)
nm <- 2
rho <- 0.2
Beta.opt <- matrix(NA, d*nm, rpt*lrs)
itr <- 0

n <- nrow(X1) 
m <- floor(n/k)
acovariate <- bresponse <- c()
K <- k
for(i in 1:K){
  if(i==K){
    acovariate <- c(acovariate,list(X1[((K-1)*m+1):nrow(X),]))
    bresponse <- c(bresponse,list(Y1[((K-1)*m+1):nrow(X)]))
  }
  else{
    acovariate <- c(acovariate,list(X1[((i-1)*m+1):(i*m),]))
    bresponse <- c(bresponse,list(Y1[((i-1)*m+1):(i*m)]))
  }
}


rm(X)
rm(Y)
rm(realdata_1,realdata_2,realdata_3,realdata_4,realdata_5,realdata_6,realdata_7,realdata_8,
          realdata_9,realdata_10,realdata_11,realdata_12,realdata_13,realdata_14,realdata_15,realdata_16,
          realdata_17,realdata_18,realdata_19,realdata_20,realdata_21,realdata_22)

gc()

set.seed(0)
for (r in r.ss) {
  cat(r, " ",k," ");
  set.seed(0)
  for (rr in 1:rpt) {
    if (rr%/%10 == rr/10) cat(rr, "-")
    itr <- itr+1
    tryCatch(
      {flag <- 1
      fit.alg.simp <- AlgTwoStp(X1,Y1,r1,r*k,rho=0.2,K=k,method="simp",acovariate,bresponse)
      Beta.simp[,itr] <- fit.alg.simp$simp
      gc()
      flag <- 0},
      error=function(e){
        cat("\n ERROR :",  conditionMessage(e), "\n")})
    if(flag != 0){
      Beta.simp[,itr] <- NA
    }
    if (any(fit.alg.simp$msg != "Successful convergence")) {
      warning(paste(rr, "not converge", fit.alg.simp$msg))
      Beta.simp[,itr] <- NA
    }
    tryCatch(
      {flag <- 1
      fit.alg1 <- AlgTwoStp(X1,Y1,r1,r*k,rho=0.2,K=k,method="mse",acovariate,bresponse)
      fit.alg2 <- AlgTwoStp(X1,Y1,r1,r*k,rho=0.2,K=k,method="mvc",acovariate,bresponse)
      gc()
      flag <- 0},
      error=function(e){
        cat("\n ERROR :",  conditionMessage(e), "\n")})
    if(flag != 0){
      Beta.opt[,itr] <- NA
      next
    }
    if (any(c(fit.alg1$msg,fit.alg2$msg) != "Successful convergence")) {
      warning(paste(rr, "not converge:", fit.alg$msg))
      next
    }
    Beta.opt[,itr] <- c(fit.alg1$opt,fit.alg2$opt)
  }
  cat("\n\n")
}

write.csv(Beta.simp,"D:/report5/data/Airlinefull_gamma_betasimp_K5.csv")
write.csv(Beta.opt,"D:/report5/data/Airlinefull_gamma_betaopt_K5.csv")


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
write.csv(rec,"D:/report5/data/Airlinefull_gamma_rec_K5.csv")

Deff <- data.frame(rho.ss=rep(r.ss,3),mse=c(log(rec)),method=rep(c("UNIF","MV","MVc"),
                                                                 each=length(r.ss)))


library(ggplot2)#275
theme_set(theme_bw())
b <- qplot(x=rho.ss, y=mse, data=Deff,  geom=c("line", "point"), shape=method, main=" ",
           xlab="r",ylab="log(MSE)",ylim = c(min(log(rec)), max(log(rec))),color=method)+geom_point(size = 3.5)
b+theme(legend.title = element_text(colour="blue", size=16),
        legend.text = element_text(size = 16, face = "bold"))




# 
# 
# pred <- function(beta,x,y){
#   l2.loss <- abs(y-((x)%*%(beta))>0)
#   (1-mean(l2.loss))
# }
# pred.simp <- apply(Beta.simp,2,pred, x=as.matrix(X.t),y=Y.t)
# pred.optmse <- apply(Beta.opt[1:d,],2,pred, x=(X.t),y=Y.t)
# pred.optmvc <- apply(Beta.opt[(d+1):(2*d),],2,pred, x=(X.t),y=Y.t)
# 
# pred.opt <- rbind(pred.optmse,pred.optmvc)
# write.csv(pred.simp,"D:/report5/data/SUSY_predsimp_K5.csv")
# write.csv(pred.opt,"D:/report5/data/SUSY_predopt_K5.csv")
# 
# 
# mspe.simp <- mspe.opt <- c()
# for (i in 1:lrs) {
#   loc <- ((i-1)*rpt+1):(i*rpt)
#   mspe.simp <- c(mspe.simp, (mean(pred.simp[loc], na.rm=TRUE)))
#   mspe.opt <- rbind(mspe.opt, apply(pred.opt[,loc],1,mean, na.rm=TRUE))
# }
# repc <- cbind(mspe.simp, mspe.opt)
# 
# write.csv(repc,paste0("D:/report5/data/SUSY_mspe_k_",as.character(k),".csv"))
# 
# 
# 
# 
# Deff <- data.frame(r.ss=rep(r.ss,3),mse=(c(repc)),method=rep(c("UNIF","MV","MVc"),
#                                                              each=length(r.ss)))
# 
# library(ggplot2)#275
# theme_set(theme_bw())
# b <- qplot(x=r.ss, y=mse, data=Deff,  geom=c("line", "point"), shape=method, main=" ",
#            xlab="K",ylab="proportion",ylim = c(min(Deff$mse), max(Deff$mse)),color=method)+geom_point(size = 3.5)
# b+theme(legend.title = element_text(colour="blue", size=16),
#         legend.text = element_text(size = 16, face = "bold"))+ geom_hline(yintercept = full.pred, size=2, color="grey",linetype=2)
# 
