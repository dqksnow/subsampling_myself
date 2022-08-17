require('mvtnorm')
require('MASS')
library(Rcpp)
library(inline)
rtools <- "C:\\Rtools\\bin"
gcc <- "C:\\Rtools\\mingw_64\\bin"
path <- strsplit(Sys.getenv("PATH"), ";")[[1]]
new_path <- c(rtools, gcc, path)
new_path <- new_path[!duplicated(tolower(new_path))]
Sys.setenv(PATH = paste(new_path, collapse = ";"))
#x(Rx) should be like renaming the numeric vector
src <- '
#include <algorithm>
#include <iostream>
Rcpp::NumericVector x(Rx);
Rcpp::NumericVector f(Rf);
Rcpp::NumericVector r(Rr);//subsample size
Rcpp::NumericVector max(Rmax);
int n = x.size(), nf = f.size();//so essentially x and f are same thing for different use
int k = r[0], j=0, mx = max[0], loc[mx];
double y[n];
for (int i = 0; i < n; i++) y[i] = x[i];
//nth_element put smaller front, larger back
std::nth_element(y, y + k - 1, y + sizeof(y)/sizeof(*y));//then y become partial sorted ver of x
double  kl = y[k-1];//the kth smallest number, used as threshold
for (int i = 0; i < n; i++) y[i] = -x[i];
std::nth_element(y, y + k - 1, y + sizeof(y)/sizeof(*y));//then y become partial sorted ver of -x


double  ku = -y[k-1];//the kth largest number, used as threshold
for (int i = 0; i < nf; i++) {
if (f[i] <= kl || f[i] >= ku)
loc[j++] = i + 1;//keep those smaller than lower thresh and those larger than upper thresh
}
Rcpp::NumericVector rt(j);
for (int l = 0; l < j; l++) {
rt[l] = loc[l];//turn int array to numbericVector?
}
return rt;
'
getIdx <- cxxfunction(signature(Rx="numeric", Rf="numeric",
                                Rr="numeric", Rmax="numeric"),
                      src, plugin = "Rcpp")




getMLE <- function(x, y, w) {
  d <- ncol(x)
  beta <- rep(0, d)
  loop  <- 1
  Loop  <- 100
  msg <- "NA"
  while (loop <= Loop) {
    pr <- c(1 - 1 / (1 + exp(x %*% beta)))
    H <- t(x) %*% (pr * (1 - pr) * w * x)
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

###generate data for logistic regression
simulate.data=function(num.covariate,num.datalines,dist)
{
  #beta = runif(num.covariate,-1,1)
  
  beta=rep(0.5,num.covariate)
  #x=rnorm(num.datalines*(num.covariate-1),0,5)
  #x.matrix=matrix(x,num.datalines,num.covariate-1)
  #x.matrix=rmvt(num.datalines,diag(num.covariate-1)*5,df=8)
  covariance=matrix(0.5,(num.covariate),(num.covariate))+diag(num.covariate)*0.5
  #x.matrix1=mvrnorm(num.datalines,rep(1,(num.covariate)),covariance)
  #x.matrix2=mvrnorm(num.datalines,rep(-1,(num.covariate)),covariance)
  #x.matrix=0.5*x.matrix1+0.5*x.matrix2
  #covariance=matrix(0.5,(num.covariate),(num.covariate))+diag(num.covariate)*0.5
  #x.matrix=mvrnorm(num.datalines,rep(1,(num.covariate)),covariance)
  # x.matrix=mvrnorm(num.datalines,rep(0,(num.covariate)),diag(num.covariate)*5)
  #x.matrix=rmvt(num.datalines,diag(num.covariate)*5,df=8)
  if(dist == 1)
  {
    x.matrix=rmvt(num.datalines,covariance,df=3)/10
  }
  if(dist == 2)
  {
    x.matrix1=mvrnorm(num.datalines,rep(1,(num.covariate)),covariance)
    x.matrix2=mvrnorm(num.datalines,rep(-1,(num.covariate)),covariance)
    x.matrix=0.5*x.matrix1+0.5*x.matrix2
  }
  if(dist == 3)
  {
    x.matrix=mvrnorm(num.datalines,rep(0,(num.covariate)),covariance)
  }
  if(dist==4)
  {
    x.matrix=mvrnorm(num.datalines,rep(1,(num.covariate)),covariance)
   
  }
  #x.matrix=cbind(rep(1,num.datalines),x.matrix)
  c=1/(1+exp(-x.matrix%*%beta))
  y=rbinom(num.datalines,1,c)
  return(list(y=y,x.matrix=x.matrix,beta=beta))
  
}


min.fun=function(x,num.covariate)
{
  value=1/(x^2*(exp(x)/((1+exp(x))^2))^(num.covariate+1))
  return(value)  
}
min.fun2=function(x,num.covariate,coef,beta)
{
  
  value = 1/(x^2*(exp(x)/(1+exp(x))^2))*(beta[num.covariate])^2+1/(exp(x)/(1+exp(x))^2)/(beta[num.covariate])^2*coef
  #value=1/(x^2*(exp(x)/((1+exp(x))^2))^num.covariate)
  return(value)  
}

second.samp=function(c,y,x.matrix,beta,num.dim,index)
{
  c_index=as.matrix(cbind(rep(1,nrow(x.matrix)),x.matrix))%*%beta
  distance1=abs(c_index-c)
  distance2=abs(c_index+c)
  a.index=(distance1>distance2)
  distance1[a.index]=0
  distance2[!a.index]=0
  distance=distance1+distance2
  initial=order(distance)[1:num.dim]
  y=y[initial]
  x.matrix=x.matrix[initial,]
  index=index[initial]
  
  return(list(y=y,x.matrix=x.matrix,index=index))
}


###optimal sample approach for logistic case; delta is prespicified boundary condition

samp.data=function(data,beta,num.covariate,num.datalines,num.samp,delta)
{
  dim=num.covariate
  x.matrix=data$x.matrix
  y=data$y
  samp.y=c()
  samp.x=data.frame()
  dim.samp=ceiling(num.samp/(2*(dim-1)))
  op.result=optimize(min.fun,num.covariate=num.covariate,interval=c(-10,10))
  c=op.result$minimum
  #print(c)
  # min.v=op.result$objective
  
  c_index=as.matrix(x.matrix)%*%beta
  distance1=abs(c_index-c)
  distance2=abs(c_index+c)
  index=(distance1>distance2)
  distance1[index]=0
  distance2[!index]=0
  distance=distance1+distance2
  
  index.update=(distance<delta)
  distance=distance[index.update]
  c_index=c_index[index.update];
  # print(c(max(c_index),min(c_index)))
  # ratio=min.v/min.fun(c_index[index.update[num]],num.covariate)
  samp.x=x.matrix[index.update,]
  # x.matrix=x.matrix[-index.update,]
  
  samp.y=y[index.update]
  x.matrix=samp.x
  y=samp.y
  
  u=length(y)
  print(u)
  gc()
  #print(u)
  samp.y=c()
  samp.x=c()
  #plot(x.matrix[,1],x.matrix[,2]);
  #for(i in 1:(dim-1))
  # {
  
  # index.update=order(x.matrix[,i])[c(1:dim.samp,(u-dim.samp+1):u)]
  #u=u-2*dim.samp
  # samp.y=c(samp.y,y[index.update])
  #samp.x=rbind(samp.x,x.matrix[index.update,])
  #y=y[-index.update]
  
  # x.matrix=x.matrix[-index.update,]
  #x.matrix=x.matrix[-]
  
  # }
  index <- c()
  use_column <- x.matrix[,1]
  # x.matrix_temp <- x.matrix
  index_raw<-c(1:nrow(x.matrix))
  for(i in 1:(dim-1))
  {
    index_pick <- getIdx(use_column, use_column, dim.samp, num.samp)
    index<- c(index,index_raw[index_pick])
    index_raw<- index_raw[-index_pick]
    use_column <- x.matrix[-index,i+1]
    #x.matrix <- x.matrix[-index_pick,]
  }
  
  #print(c(length(index),length(unique(index)),length(index)!=length(unique(index))))
  samp.y<-y[index]
  samp.x <- x.matrix[index,]
  
  print(length(samp.y))
  
  #plot(samp.x[,1],samp.x[,2])
  
  return(list(y=samp.y,x.matrix=samp.x))
}




samp.data.Aoptimal=function(data,beta,num.covariate,num.datalines,num.samp,delta)
{
  dim=num.covariate
  x.matrix=data$x.matrix
  y=data$y
  #print(ncol(x.matrix))
  max = apply(x.matrix,1,max)
  min = apply(x.matrix,1,min)
  coef = (1+sum(4/(max[-num.covariate]-min[-num.covariate])^2))
  samp.y=c()
  samp.x=data.frame()
  dim.samp=ceiling(num.samp/(2*(dim-1)))
  op.result=optimize(min.fun2,num.covariate=num.covariate,coef=coef,beta=beta,interval=c(0,1))
  c=op.result$minimum
  #print(c)
  # min.v=op.result$objective
  
  c_index=as.matrix(x.matrix)%*%beta
  distance1=abs(c_index-c)
  distance2=abs(c_index+c)
  index=(distance1>distance2)
  distance1[index]=0
  distance2[!index]=0
  distance=distance1+distance2
  
  index.update=(distance<delta)
  distance=distance[index.update]
  c_index=c_index[index.update];
  # print(c(max(c_index),min(c_index)))
  # ratio=min.v/min.fun(c_index[index.update[num]],num.covariate)
  samp.x=x.matrix[index.update,]
  # x.matrix=x.matrix[-index.update,]
  
  samp.y=y[index.update]
  x.matrix=samp.x
  y=samp.y
  
  u=length(y)
  print(u)
  #print(u)
  samp.y=c()
  samp.x=c()
  #plot(x.matrix[,1],x.matrix[,2]);
  for(i in 1:(dim-1))
  {
    
    index.update=order(x.matrix[,i])[c(1:dim.samp,(u-dim.samp+1):u)]
    u=u-2*dim.samp
    samp.y=c(samp.y,y[index.update])
    samp.x=rbind(samp.x,x.matrix[index.update,])
    y=y[-index.update]
    
    x.matrix=x.matrix[-index.update,]
    #x.matrix=x.matrix[-]
    
  }
  
  #plot(samp.x[,1],samp.x[,2])
  
  return(list(y=samp.y,x.matrix=samp.x))
}




###optimal sample approach for logistic case; delta is prespicified boundary condition

samp.data.LEV=function(data,beta,num.covariate,num.datalines,num.samp)
{
  #dim=num.covariate-1
  x.matrix=data$x.matrix
  #beta = beta[-1]
  
  y=data$y
  samp.y=c()
  samp.x=data.frame()
  x.norm=apply(x.matrix^2,1,sum)
  x.norm = sqrt(x.norm)
  coef = abs(y-1/(1+exp(-x.matrix%*%beta)))
  weight = coef*x.norm/(sum(coef*x.norm));
  
  index = sample(1:num.datalines,num.samp,prob = weight,replace = TRUE)
  #u=length(y)
  #print(u)
  #print(u)
  samp.y=c()
  samp.x=c()
  samp.x=x.matrix[index,]
  samp.y = y[index]
  samp.weight = weight[index]
  #plot(samp.x[,1],samp.x[,2])
  
  return(list(y=samp.y,x.matrix=samp.x,weight =samp.weight))
}
###simulations to test the performance of optimal sampling , sequential sampling and random sampling
full.regression=function(iter,num.covariate,num.datalines,num.samp,delta,initial,dist)
{
  MSE.samp3=rep(0,iter)
  MSE.samp=rep(0,iter)
  MSE.samp2=rep(0,iter)
  MSE = rep(0,iter)
  user_time_new = rep(0,iter)
  user_time_srs = rep(0,iter)
  user_time_lev =rep(0,iter)
  user_time_full =rep(0,iter)
  sys_time_new = rep(0,iter)
  sys_time_srs = rep(0,iter)
  sys_time_lev =rep(0,iter)
  sys_time_full =rep(0,iter)
  # H_matrix4= cbind(c(1,-1,-1,1),c(1,1,-1,-1),c(1,-1,1,-1),c(1,1,1,1));
  #H_matrix8 = rbind(cbind(H_matrix4,H_matrix4),cbind(-H_matrix4,H_matrix4));
  #op.result=optimize(min.fun,num.covariate=num.covariate,interval=c(-10,10))
  #c=op.result$minimum
  #max.distance=c()
  i=1
  count_new =count_lev=count_srs=count_full=0
  true_iter <- iter
  while(i <= iter)
  {
    data.simu=simulate.data(num.covariate,num.datalines,dist)
    #data.reg=data.frame(y=data.simu$y,data.simu$x.matrix[,-1])
    data.reg=data.frame(y=data.simu$y,data.simu$x.matrix)
    beta.true=data.simu$beta
    #x.matrix=data.simu$x.matrix[,-1]
    x.matrix=data.simu$x.matrix
    y=data.simu$y
    n1<- sum(y)
    n0 <- length(y)-sum(y)
    weight_initial <- rep(1/n0,n0+n1)
    weight_initial[y==1]<- 1/n1
    pick=sample(1:num.datalines,initial,prob = weight_initial,replace = TRUE)
    #new.data=data.frame(y=data$y,x.matrix))
    
    random.samp1=data.frame(ky=y[pick],x.matrix[pick,])
    #fit.random.samp1=glm(ky~.-1,random.samp1,family = "binomial")
    #beta=fit.random.samp1$coefficients
    fit.random.samp1=getMLE(as.matrix(random.samp1[,-1]),random.samp1$ky,w=1/weight_initial[pick])
    beta=fit.random.samp1$par
    
    
    
    #data.sample.new=samp.data(data.simu,beta,num.covariate,num.datalines,num.samp2)
    # plot_c (data.sample.new, beta.true);
    
    start1 <- proc.time()
    #data.sample=samp.data4(data.simu,num.covariate,num.datalines,num.samp,multiplyer)
    data.sample=samp.data(data.simu,beta,num.covariate,num.datalines,num.samp,delta)
    #plot_c (data.sample, beta.true);
    
    #data.sample=samp.data(data.simu,num.covariate,num.datalines,num.samp)
    data.reg.samp=data.frame(y=data.sample$y,data.sample$x.matrix)
    # fit=glm(y~.,data.reg,family = "binomial")
    fit.samp=getMLE(as.matrix(data.reg.samp[,-1]),data.reg.samp$y,w= 1)
    #fit.samp2=glm(y~-1+.,data.reg.samp2,family = "binomial",weights = 1/data.sample2$weight/num.samp)
    #fit.samp2=glm(y~-1+.,data.reg.samp2,family = "binomial")
    MSE.samp[i]=sum((fit.samp$par-beta.true)^2)
    time_new<- proc.time()-start1
    sys_time_new[i] <-time_new[2]
    user_time_new[i] <- time_new[1]
    #MSE[i]=sum((fit$coefficients[-1]-beta.true[-1])^2)
   
    #data.sample=samp.data(data.simu,num.covariate,num.datalines,num.samp)
    #data.reg.samp=data.frame(y=data.sample$y,data.sample$x.matrix)
    #fit.samp=glm(y~.,data.reg.samp,family = "binomial")
    start2 <- proc.time()
    data.sample2=samp.data.LEV(data.simu,beta,num.covariate,num.datalines,num.samp)
    
    #data.sample=samp.data(data.simu,num.covariate,num.datalines,num.samp)
    data.reg.samp2=data.frame(y=data.sample2$y,data.sample2$x.matrix)
    # fit=glm(y~.,data.reg,family = "binomial")
    fit.samp2=getMLE(as.matrix(data.reg.samp2[,-1]),data.reg.samp2$y,w= 1/data.sample2$weight)
    #fit.samp2=glm(y~-1+.,data.reg.samp2,family = "binomial",weights = 1/data.sample2$weight/num.samp)
    #fit.samp2=glm(y~-1+.,data.reg.samp2,family = "binomial")
    MSE.samp2[i]=sum((fit.samp2$par-beta.true)^2)
    time_lev <- (proc.time()-start2)
    sys_time_lev[i] <- time_lev[2]
    user_time_lev[i] <- time_lev[1] 
    #plot_c (data.sample3, beta.true);
    #plot_c (data.sample3, beta.true);
    
    start3 <- proc.time()
    data.reg.samp3=data.reg[sample(1:num.datalines,num.samp+initial,replace = TRUE),]
    
    #data.frame(y=data.sample2$y,data.sample2$x.matrix)
    fit.samp3<-getMLE(as.matrix(data.reg.samp3[,-1]),data.reg.samp3[,1],w=1)
    MSE.samp3[i]=sum((fit.samp3$par-beta.true)^2)
    time_srs <-  proc.time()-start3
    sys_time_srs[i] <-  time_srs[2]
    user_time_srs[i] <-  time_srs[1]
    
    
    start4 <- proc.time()
    
    #data.frame(y=data.sample2$y,data.sample2$x.matrix)
    #fit<-getMLE(as.matrix(data.reg[,-1]),data.reg[,1],w=1)
    #MSE[i]=sum((fit$par-beta.true)^2)
    time_full <-  proc.time()-start4
    sys_time_full[i] <-  time_full[2]
    user_time_full[i] <-  time_full[1]
    #plot_c2 (data.reg.samp2, beta.true);
    #print(c(nrow(data.reg.samp),nrow(data.reg.samp2),nrow(data.reg.samp3)))
    #print(c(log(MSE.samp[i]),log(MSE.samp2[i]),log(MSE.samp3[i])))
    print(c(log(MSE.samp[i]),log(MSE.samp2[i]),log(MSE.samp3[i])))
    print(rbind(time_srs=time_srs,time_lev=time_lev,time_full=time_full,time_new=time_new))
    if(max(c(MSE.samp[i],MSE.samp2[i],MSE.samp3[i]))>20)
    {
      
      if(MSE.samp[i]>20)
      {
        count_new =count_new+1
      }
      if(MSE.samp2[i]>20)
      {
        count_lev =count_lev+1
      }
      if(MSE.samp3[i]>20)
      {
        count_srs =count_srs+1
      }
      if(MSE[i]>20)
      {
        count_full =count_full+1
      }
      true_iter<- true_iter+1
      i=i-1;
      
    }
    i=i+1
    print(i)
    # if(MSE.samp < (MSE.random.samp/30))
    #{
    #  print(MSE.samp)
    # print(MSE.random.samp)
    #c_index_marked=cbind(rep(1,nrow(data.sample$x.matrix)),data.sample$x.matrix)%*%beta.true
    #distance1=abs(c_index_marked-c)
    #distance2=abs(c_index_marked+c)
    #index=(distance1>distance2)
    #distance1[index]=0
    #distance2[!index]=0
    #  distance=distance1+distance2
    # max.distance=c(max.distance,max(distance))
    #i=i+1
    
    #}
  } 
  #print(max.distance)
  return(list(MSE=list(MSE_full=mean(MSE),MSE_lev=mean(MSE.samp2),MSE_srs=mean(MSE.samp3),MSE_new=mean(MSE.samp)),
              time=list(sys_time_full=mean(sys_time_full),sys_time_lev=mean(sys_time_lev),sys_time_srs = mean(sys_time_srs),sys_time_new=mean(sys_time_new),
                        user_time_full=mean(user_time_full),user_time_lev=mean(user_time_lev),user_time_srs = mean(user_time_srs),user_time_new=mean(user_time_new)),
              count = list(count_full =count_full,count_lev =count_lev,count_srs =count_srs,
                           count_new =count_new),true_iter=true_iter
  ))
  
  
}






iter=1000
set.seed(10)
num.covariate=c(7)
num.datalines=c(500000)
initial=1000
delta=c(0.5,0.5,0.5,2.5)
dist =c(1:4)
distribution = c("T3 Distribution", "MixNormal DSistribution", "MzNormal Distribution","NzNormal Distribution")
num.samp =c(1000,2000,5000,8000)

MSE_full=MSE_srs=MSE_lev=MSE_new=user_time_full=user_time_srs=user_time_lev=user_time_new=
  sys_time_full=sys_time_srs=sys_time_lev=sys_time_new=
  count_full=count_srs=count_lev=count_new=true_iter=rep(0,length(num.samp)*length(dist))
for(i in 1:length(num.samp)){
  for(j in 1:length(dist))
  {
    a=full.regression(iter,num.covariate,num.datalines,num.samp[i],delta[j],initial,dist[j])
    MSE_full[(i-1)*length(dist)+j] <- a$MSE$MSE_full
    MSE_lev[(i-1)*length(dist)+j] <- a$MSE$MSE_lev
    MSE_srs[(i-1)*length(dist)+j] <- a$MSE$MSE_srs
    MSE_new[(i-1)*length(dist)+j] <- a$MSE$MSE_new
    
    
    sys_time_full[(i-1)*length(dist)+j] <- a$time$sys_time_full
    sys_time_lev[(i-1)*length(dist)+j] <- a$time$sys_time_lev
    sys_time_srs[(i-1)*length(dist)+j] <- a$time$sys_time_srs
    sys_time_new[(i-1)*length(dist)+j] <- a$time$sys_time_new
    user_time_full[(i-1)*length(dist)+j] <- a$time$user_time_full
    user_time_lev[(i-1)*length(dist)+j] <- a$time$user_time_lev
    user_time_srs[(i-1)*length(dist)+j] <- a$time$user_time_srs
    user_time_new[(i-1)*length(dist)+j] <- a$time$user_time_new
    
    
    count_full[(i-1)*length(dist)+j] <- a$count$count_full
    count_lev[(i-1)*length(dist)+j] <- a$count$count_lev
    count_srs[(i-1)*length(dist)+j] <- a$count$count_srs
    count_new[(i-1)*length(dist)+j] <- a$count$count_new
    
    true_iter[(i-1)*length(dist)+j] <- a$true_iter
  }
  
}

result <- data.frame(delta = delta,num.datalines = num.datalines,num.samp=kronecker(num.samp,rep(1,4)),initial=initial,
                     MSE_full=MSE_full, MSE_lev=MSE_lev, MSE_srs = MSE_srs, MSE_new=MSE_new,
                     sys_time_full=sys_time_full, sys_time_lev=sys_time_lev, sys_time_srs=sys_time_srs, sys_time_new=sys_time_new,
                     user_time_full=user_time_full, user_time_lev=user_time_lev, user_time_srs=user_time_srs, user_time_new=user_time_new,
                     count_full=count_full, count_lev=count_lev, count_srs=count_srs, count_new=count_new,true_iter=true_iter,Distribution = rep(distribution,4))
setwd("E:/DESKTOP/big data/DIFFERENT DISTRIBUTION")
write.csv(result,file = "result04_05_2018 large sample size different distribution.txt")



result <- read.table("result09_09_2017 large sample size different distribution.txt",sep = ",",header = TRUE)
library(ggplot2)

result <- rbind(result,result,result)
result$MSE <- result$MSE_lev
result$MSE[33:48] <- result $ MSE_new[1:16]
result$MSE[17:32] <- result $ MSE_srs[1:16]
result$Algorithm <- "Leveraging Algorithm (mVc)"
result$Algorithm[33:48] <- "New Algorithm"
result$Algorithm[17:32] <- "Simple Random Sampling"
result$Log_MSE <- log(result$MSE)
ggplot(result,aes(x = num.samp,y =Log_MSE,by = Algorithm, color = Algorithm))+
  geom_point(aes( size=2,shape = Algorithm))+
  geom_line()+facet_wrap(~Distribution,nrow = 2,scales = "free_y")+
  theme(legend.title = element_text(colour="blue", size=15, face="bold"),
        legend.text = element_text(colour="blue", size = 10, face = "bold"))+
  labs(x = "Size of Subsample")+ggtitle("Comparison of subsampling algorithms under different distribution")

ggsave("big size.png",dpi = 600)
