# Note:

# If you want to get the results of all methods and all criteria for certain
# full-data at once, please use this code.

# Please load all function(including 'subsampling_simu' in this code),
# and choose a type of data to generate,
# and then run the simu code(line 138 in this code).

# For example, run line 11 to line 50, and then run line 138 to the end.
source("R/prior.R")
source('R/logistic_optimal_subsampling.R')
source('R/ossp_log.R')
library(mvtnorm)
#install.packages('parallel')
library(parallel)

subsampling_simu <- function(j){ #this function is used for parLapply
  source("R/prior.R")
  source('R/logistic_optimal_subsampling.R')
  source('R/ossp_log.R')
  beta_temp <- matrix(0, S, d)
  r <- num_r[j]
  for(i in 1:S){
    set.seed(i)
    result <- logistic_optimal_subsampling(X, y, r0, r,
                                           criteria = criteria,
                                           method = method,
                                           unweighted.estimator = unweighted.estimator,
                                           b = 2)
    beta_temp[i,] <- result$beta
  }
  mse <- sum((beta_temp - beta_full)^2)/S # I just collect mse, not beta_temp
  return(mse)
}

#### JASA2018 simu1 multinormal balanced ########
rm(list = setdiff(ls(), lsf.str())) # remove all variables except functions
title <- 'simu1 multinormal balanced '
N <- 1e4
beta0  <- c(rep(0.5, 7))
d <- length(beta0)
corr  <- 0.5
sigmax  <- matrix(corr, d-1, d-1) + diag(1-corr, d-1)
set.seed(123)
X <- rmvnorm(N, rep(0, d-1), sigmax)
X <- cbind(1, X)
P <- 1 - 1 / (1 + exp(X %*% beta0))
y <- rbinom(N, 1, P)
beta_full <- logistic_coef_estimate(X, y, 1, 1:length(y))

#### JASA2018 simu2 multinormal imbalanced ########
rm(list = setdiff(ls(), lsf.str())) # remove all variables except functions
title <- 'simu2 multinormal imbalanced '
N <- 1e4
beta0  <- c(rep(0.5, 7))
d <- length(beta0)
corr  <- 0.5
sigmax  <- matrix(corr, d, d) + diag(1-corr, d)
set.seed(123)
X <- rmvnorm(N, rep(1.5, d), sigmax)
#X <- cbind(1, X)
P <- 1 - 1 / (1 + exp(X %*% beta0))
y <- rbinom(N, 1, P)
#table(y)
beta_full <- logistic_coef_estimate(X, y, 1, 1:length(y))

#### JASA2018 simu3 hetero multinormal ########
rm(list = setdiff(ls(), lsf.str())) # remove all variables except functions
title <- 'simu3 hetero multinormal'
N <- 1e4
beta0  <- c(rep(0.5, 7))
d <- length(beta0)
#var <- 1/c(1:d)^2 #
var <- c(1:(d)) #^2
corr <- 0.5
sigmax <- matrix(corr, d, d)  + diag(as.vector(var)) - diag(corr, d)
set.seed(123)
X <- rmvnorm(N, rep(0, d), sigmax)
# note: if var <- 1/c(1:(d-1))^2 as in JASA2018,
# there will be a warning: sigma is numerically not positive semidefinite
# and subsampling will have an error.
# so I change var to var <- c(1:(d-1))
#X <- cbind(1, X)
P <- 1 - 1 / (1 + exp(X %*% beta0))
y <- rbinom(N, 1, P)
table(y)
beta_full <- logistic_coef_estimate(X, y, 1, 1:length(y))

#### JASA2018 simu4 biomodal normal ########
rm(list = setdiff(ls(), lsf.str())) # remove all variables except functions
title <- 'simu4 biomodal normal'
N <- 1e4
beta0  <- c(rep(0.5, 7))
d <- length(beta0)
corr  <- 0.5
sigmax  <- matrix(corr, d, d) + diag(1-corr, d)
set.seed(123)
X <- 0.5*rmvnorm(N, rep(1, d), sigmax) + 0.5*rmvnorm(N, rep(-1, d), sigmax)
# X <- cbind(1, X)
P <- 1 - 1 / (1 + exp(X %*% beta0))
y <- rbinom(N, 1, P)
table(y)
beta_full <- logistic_coef_estimate(X, y, 1, 1:length(y))

#### JASA2018 simu5 multi t ########
rm(list = setdiff(ls(), lsf.str())) # remove all variables except functions
title <- "simu5 multi t"
N <- 1e4
beta0  <- c(rep(0.5, 7))
d <- length(beta0)
corr  <- 0.5
sigmax  <- matrix(corr, d, d) + diag(1-corr, d)
set.seed(123)
X <-rmvt(N, sigmax, df=3, rep(0, d), type=c("shifted","Kshirsagar")[1])#in JASA2018, it is "X/10". maybe that is a typo.
X <- cbind(1, X)
P <- 1 - 1 / (1 + exp(X %*% beta0))
y <- rbinom(N, 1, P)
table(y)
beta_full <- logistic_coef_estimate(X, y, 1, 1:length(y))

#### JASA2018 simu6 exp ########
rm(list = setdiff(ls(), lsf.str())) # remove all variables except functions
title <- "simu6 exp"
N <- 1e4
beta0  <- c(rep(0.5, 7))
d <- length(beta0)
set.seed(123)
X <- matrix(0, N, d)
generate_rexp <- function(x) x<-rexp(N, rate=2)
X <- apply(X, 2, generate_rexp)
#X <- cbind(1, X)
P <- 1 - 1 / (1 + exp(X %*% beta0))
y <- rbinom(N, 1, P)
table(y)
beta_full <- logistic_coef_estimate(X, y, 1, 1:length(y))

#### simu code ######
#### please load all function, generate a type of simu data and then run below
S <- 1000
r0 <- 200
r_min <- 200
r_max <- 1000
num_r <- c(100, seq(r_min, r_max, 200)) # different r to be considered
weight_parameter <-  c(F, T)
criteria_all <- c("optA", "optL", "LCC")
method_all <- c("SWR", "Poisson")
mse <- c()
cl <- makeCluster(length(num_r)) # how many CPU cores are called
t1=proc.time()
for (w in 1:length(weight_parameter)){
  unweighted.estimator <- weight_parameter[w]
  for (p in 1:length(method_all)){
    method <- method_all[p]
    mse_criteria <- matrix(0, length(criteria_all), length(num_r))  # collect mse for different criteria in one method
    for (q in 1:length(criteria_all)){
      criteria <- criteria_all[q]
      clusterExport(cl=cl,
                    varlist=c("criteria", "method", 'unweighted.estimator',
                              "X", "y", "r0", 'S', 'd', "beta_full", "num_r"),
                    envir=environment()) #import environment variables into the function 'subsampling_simu'
      results <- parLapply(cl, 1:length(num_r), subsampling_simu) #only output the mse
      for(i in 1:length(num_r)){
        mse_criteria[q,i] <- results[[i]]
      }
      t2=proc.time()
      time.cost=t2-t1
      print(paste0('unweighted = ', unweighted.estimator, ', ', method, ' + ',
                   criteria, ', total time: ',time.cost[3][[1]],'s'))
    }
    mse <- rbind(mse, mse_criteria)
  }
}
mse_dataframe <- as.data.frame(mse,
                               row.names = c('we_SWR_optA', 'we_SWR_optL', 'we_SWR_LCC',
                                           'we_Poi_optA', 'we_Poi_optL', 'we_Poi_LCC',
                                           'unwe_SWR_optA', 'unwe_SWR_optL', 'unwe_SWR_LCC',
                                           'unwe_Poi_optA', 'unwe_Poi_optL', 'unwe_Poi_LCC'))
colnames(mse_dataframe) <- as.character(num_r)
# row names are different method, col names are different r to be considered
print(title) # print current simu data type
#mse_dataframe
#log(mse_dataframe)

####  plot  ##################
library(ggplot2)
library(reshape2)
mse_plot <- data.frame(r = num_r
                       , we_SWR_optA = mse[1,]
                       , we_SWR_optL = mse[2,]
                       , we_SWR_LCC = mse[3,]

                       , we_Poi_optA = mse[4,]
                       , we_Poi_optL = mse[5,]
                       , we_Poi_LCC = mse[6,]

                       , unwe_SWR_optA = mse[7,]
                       , unwe_SWR_optL = mse[8,]
                       , unwe_SWR_LCC = mse[9,]

                       , unwe_Poi_optA = mse[10,]
                       , unwe_Poi_optL = mse[11,]
                       , unwe_Poi_LCC = mse[12,]
                       )
# If you don't want to draw so many methods on the same plot,
# you can comment out certain lines.
mse_plot<-melt(mse_plot, id.vars = 'r', variable.name="method",
  value.name="mse") # convert data structure to fit ggplot
# mse_plot$mse <- log(mse_plot$mse) # take log(mse)
mse_figure <- ggplot(mse_plot, aes(x=r, y=mse)) +
  geom_line(aes(colour = method)) + geom_point() +
  scale_x_continuous(breaks = num_r)
mse_figure
