rm(list = ls())
getMLE <- function(x, y, w) {
  beta <- c(2.59444145,-0.02597514,0.11718343,0.34156773,0.84544824,0.60636395,-0.22549084,0.49362886,0.11556359 )
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
    if(tlr < 0.0001) {
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
library(data.table)

realdata_article1 <- fread("D:/讨论班/report5/data/citationjoint0.csv", header = T)
realdata_article2 <- fread("D:/讨论班/report5/data/citationjoint1.csv", header = T)
realdata_article3 <- fread("D:/讨论班/report5/data/citationjoint2.csv", header = T)
realdata_article4 <- fread("D:/讨论班/report5/data/citationjoint3.csv", header = T)

article_info <- rbind(realdata_article1,realdata_article2,realdata_article3,realdata_article4)
rm(realdata_article1,realdata_article2,realdata_article3,realdata_article4)
gc()
#article_info <- fread("D:/讨论班/report5/data/citationjoint3.csv", header = T)
article_info <- article_info[complete.cases(article_info),-1]
article_info <- unique( article_info )
Y <- article_info$citation
Year <- 2018-as.numeric(article_info$year)
abs <- ifelse(article_info$abs_words>100,1,0)
abs1 <- ifelse(article_info$abs_words<=100&article_info$abs_words>0,1,0)
titleinfo <- ifelse(article_info$num_title>9,1,0)
#teaminfo <- ifelse(article_info$number_authors>3,1,0)
quality  <-  ifelse(article_info$`SJR Best Quartile`%in%c('Q1'),1,0)
isjournal <- ifelse(article_info$Type=='journal',1,0)
X <- cbind(article_info$SJR,article_info$authorpaper2,Year)
X1 <- cbind(1,scale(X),abs,abs1,titleinfo,quality,isjournal)

#model1 <- glm(Y~-1+X1,family = poisson())




# nrow(rec1)
# [1] 3139874

d<- dim(X1)[2]
X <- as.matrix(X1)
fit.full <- getMLE(X, Y, 1)
#model1 <- glm(Y~-1+X,family = poisson())
### Sampling procedure
r1 <- 800
k <- 5
r.ss <- c(1000,1200,1500,1700,2000,2200)*2

set.seed(0)
rpt <- 1000
#r.ss <- c( 300, 500, 700, 1000, 1200, 1400)## [5:6]
lrs <- length(r.ss)
Beta.simp <- matrix(NA, d, rpt*lrs)
nm <- 2
rho <- 0.2
Beta.opt <- matrix(NA, d*nm, rpt*lrs)
itr <- 0

for (r in r.ss) {
  cat(r, " ",k," ");
  set.seed(0)
  for (rr in 1:rpt) {
    if (rr%/%10 == rr/10) cat(rr, "-")
    itr <- itr+1
    tryCatch(
      {flag <- 1
      fit.alg.simp <- AlgTwoStp(X,Y,r1,r*k,rho=0.2,K=k,method="simp")
      Beta.simp[,itr] <- fit.alg.simp$simp 
      flag <- 0}, 
      error=function(e){
        cat("\n ERROR :", itr, conditionMessage(e), "\n")})
    if(flag != 0){
      Beta.simp[,itr] <- NA
    }
    if (any(fit.alg.simp$msg != "Successful convergence")) {
      warning(paste(rr, "not converge", fit.alg.simp$msg))
      Beta.simp[,itr] <- NA
    }
    tryCatch(
      {flag <- 1
      fit.alg <- AlgTwoStp(X,Y,r1,r*k,rho=0.2,K=k,method="opt")
      flag <- 0}, 
      error=function(e){
        cat("\n ERROR :", itr, conditionMessage(e), "\n")})
    if(flag != 0){
      Beta.opt[,itr] <- NA
      next
    }
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

write.csv(Beta.simp,"D:/report5/data/cite_betasimp_K5_withoutimpact.csv")
write.csv(Beta.opt,"D:/report5/data/cite_betaopt_K5_withoutimpact.csv")

mse.simp <- mse.opt <- c()
for (i in 1:lrs) {
  loc <- ((i-1)*rpt+1):(i*rpt)
  mse.simp <- c(mse.simp, sum(rowMeans(l2.simp[,loc], na.rm=TRUE)))
  mse.opt <- rbind(mse.opt, colSums(matrix(
    rowMeans(l2.opt[,loc], na.rm=TRUE), d, nm)))
}



rec <- cbind(mse.simp, mse.opt)
write.csv(rec,"D:/report5/data/cite_rec_K5_withoutimpact.csv")

Deff <- data.frame(rho.ss=rep(r.ss,3),mse=c(log(rec)),method=rep(c("UNIF","MV","MVc"),
                                                                 each=length(r.ss)))

library(ggplot2)#275
theme_set(theme_bw())
b <- qplot(x=rho.ss, y=mse, data=Deff,  geom=c("line", "point"), shape=method, main=" ",
           xlab="r",ylab="log(MSE)",ylim = c(min(log(rec)), max(log(rec))),color=method)+geom_point(size = 3.5)
b+theme(legend.title = element_text(colour="blue", size=16),
        legend.text = element_text(size = 16, face = "bold"))
