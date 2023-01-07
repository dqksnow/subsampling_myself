# Note:

# If you want to get the results of all methods and all criteria at once, please use this code.

# choose a type of simulation data to generate,and then run the simu code.

source("R/prior.R")
source('R/logistic_optimal_subsampling.R')
source('R/ossp_log.R')
source('draft/generate_simu_data.R')
library(mvtnorm)
#install.packages('parallel')
library(parallel)
library(ggplot2)
library(reshape2)

subsampling_simu <- function(j){ #this function is used for parLapply
  source("R/prior.R")
  source('R/logistic_optimal_subsampling.R')
  source('R/ossp_log.R')
  source('draft/generate_simu_data.R')
  beta_temp <- matrix(0, S, d)
  correct_proportion <- rep(0, S)
  r <- num_r[j] # for given criteria and method, run 'logistic_optimal_subsampling' for the same 'r' for S replications
  for(i in 1:S){
    set.seed(i)
    result <- logistic_optimal_subsampling(X, y, r0, r,
                                           criteria = criteria,
                                           method = method,
                                           unweighted.estimator = unweighted.estimator,
                                           b = 2)
    beta_temp[i,] <- result$beta
    correct_proportion[i] <- predict_classification(X_test, y_test, result$beta)
  }
  mse <- sum((beta_temp - beta_full)^2)/S
  correct_proportion_mean <- mean(correct_proportion)
  return(list(mse=mse, proportion=correct_proportion_mean))
}

#### JASA2018 simu1 multinormal balanced ########
rm(list = setdiff(ls(), lsf.str())) # remove all variables except functions
dev.off() # remove the plot
title <- 'simu1_multinormal_balanced'
simu_data <- simu_mzNormal(seed=123, N=1e4, beta0=c(rep(0.5, 7)), corr=0.5)
X <- simu_data$X
y <- simu_data$y
d <- ncol(X)
beta_full <- simu_data$beta_full
test_data <- simu_mzNormal(seed=124, N=1e4, beta0=c(rep(0.5, 7)), corr=0.5)
X_test <- test_data$X
y_test <- test_data$y
correct_proportion_full <- predict_classification(X_test, y_test, beta_full)

#### JASA2018 simu2 multinormal imbalanced ########
rm(list = setdiff(ls(), lsf.str()))
dev.off() # remove the plot
title <- 'simu2_multinormal_imbalanced'
simu_data <- simu_nzNormal(seed=123, N=1e4, beta0=c(rep(0.5, 7)), corr=0.5)
X <- simu_data$X
y <- simu_data$y; table(y)
d <- ncol(X)
beta_full <- simu_data$beta_full
test_data <- simu_nzNormal(seed=124, N=1e4, beta0=c(rep(0.5, 7)), corr=0.5)
X_test <- test_data$X
y_test <- test_data$y
correct_proportion_full <- predict_classification(X_test, y_test, beta_full)

#### JASA2018 simu3 hetero multinormal ########
rm(list = setdiff(ls(), lsf.str()))
dev.off() # remove the plot
title <- 'simu3_hetero_multinormal'
simu_data <- simu_ueNormal(seed=123, N=1e4, beta0=c(rep(0.5, 7)), corr=0.5)
X <- simu_data$X
y <- simu_data$y; table(y)
d <- ncol(X)
beta_full <- simu_data$beta_full
test_data <- simu_ueNormal(seed=124, N=1e4, beta0=c(rep(0.5, 7)), corr=0.5) #generate a new full data
X_test <- test_data$X
y_test <- test_data$y
correct_proportion_full <- predict_classification(X_test, y_test, beta_full)
# note: if var <- 1/c(1:(d-1))^2 as in JASA2018,
# there will be a warning: sigma is numerically not positive semidefinite
# and subsampling will have an error.
# so I change var to var <- c(1:(d-1))


#### JASA2018 simu4 biomodal normal ########
rm(list = setdiff(ls(), lsf.str()))
dev.off() # remove the plot
title <- 'simu4_biomodal_normal'
simu_data <- simu_mixNormal(seed=123, N=1e4, beta0=c(rep(0.5, 7)), corr=0.5)
X <- simu_data$X
y <- simu_data$y; table(y)
d <- ncol(X)
beta_full <- simu_data$beta_full
test_data <- simu_mixNormal(seed=124, N=1e4, beta0=c(rep(0.5, 7)), corr=0.5)
X_test <- test_data$X
y_test <- test_data$y
correct_proportion_full <- predict_classification(X_test, y_test, beta_full)

#### JASA2018 simu5 multi t ########
rm(list = setdiff(ls(), lsf.str()))
dev.off() # remove the plot
title <- "simu5_multi_t"
simu_data <- simu_T3(seed=123, N=1e4, beta0=c(rep(0.5, 7)), corr=0.5, df=3)
X <- simu_data$X
y <- simu_data$y; table(y)
d <- ncol(X)
beta_full <- simu_data$beta_full
test_data <- simu_T3(seed=124, N=1e4, beta0=c(rep(0.5, 7)), corr=0.5, df=3)
X_test <- test_data$X
y_test <- test_data$y
correct_proportion_full <- predict_classification(X_test, y_test, beta_full)

#### JASA2018 simu6 exp ########
rm(list = setdiff(ls(), lsf.str())) # remove all variables except functions
dev.off() # remove the plot
title <- "simu6_exp"
simu_data <- simu_EXP(seed=123, N=1e4, beta0=c(rep(0.5, 7)), rate=2)
X <- simu_data$X
y <- simu_data$y
d <- ncol(X)
beta_full <- simu_data$beta_full
test_data <- simu_EXP(seed=124, N=1e4, beta0=c(rep(0.5, 7)), rate=2)
X_test <- test_data$X
y_test <- test_data$y
correct_proportion_full <- predict_classification(X_test, y_test, beta_full)

#### simu code ######
#### please load all function, generate a type of simu data and then run below
S <- 1000
r0 <- 200
r_min <- 200
r_max <- 1000
num_r <- c(100, seq(r_min, r_max, 200)) # different r to be considered
weight_parameter <-  c(F)
#weight_parameter <-  c(F, T)
criteria_all <- c("optA", "optL", "LCC")
method_all <- c("SWR")
#method_all <- c("SWR", "Poisson")
mse <- proportion <- c()

cl <- makeCluster(length(num_r)) # how many CPU cores are called
t1=proc.time()
for (w in 1:length(weight_parameter)){
  unweighted.estimator <- weight_parameter[w]
  for (p in 1:length(method_all)){
    method <- method_all[p]
    mse_criteria <- proportion_criteria<- matrix(0, length(criteria_all), length(num_r))  # collect mse for different criteria under the same method
    for (q in 1:length(criteria_all)){
      criteria <- criteria_all[q]
      clusterExport(cl=cl,
                    varlist=c("criteria", "method", 'unweighted.estimator',
                              "X", "y", "r0", 'S', 'd', "beta_full", "num_r",
                              'X_test', 'y_test'),
                    envir=environment()) #import environment variables into the function 'subsampling_simu'
      results <- parLapply(cl, 1:length(num_r), subsampling_simu)
      for(i in 1:length(num_r)){
        mse_criteria[q,i] <- results[[i]][[1]]
        proportion_criteria[q,i] <- results[[i]][[2]]
      }
      t2=proc.time()
      time.cost=t2-t1
      print(paste0('unweighted = ', unweighted.estimator, ', ', method, ' + ',
                   criteria, ', total time: ',time.cost[3][[1]],'s'))
    }
    mse <- rbind(mse, mse_criteria)
    proportion <- rbind(proportion, proportion_criteria)
  }
}
mse_dataframe <- as.data.frame(mse,
                               row.names = c('we_SWR_optA', 'we_SWR_optL', 'we_SWR_LCC'))#,
                                           #'we_Poi_optA', 'we_Poi_optL', 'we_Poi_LCC',
                                           #'unwe_SWR_optA', 'unwe_SWR_optL', 'unwe_SWR_LCC',
                                           #'unwe_Poi_optA', 'unwe_Poi_optL', 'unwe_Poi_LCC'))
colnames(mse_dataframe) <- as.character(num_r)
proportion_dataframe <- as.data.frame(proportion,
                               row.names = c('we_SWR_optA', 'we_SWR_optL', 'we_SWR_LCC'))#,
                                            #'we_Poi_optA', 'we_Poi_optL', 'we_Poi_LCC',
                                            #'unwe_SWR_optA', 'unwe_SWR_optL', 'unwe_SWR_LCC',
                                            #'unwe_Poi_optA', 'unwe_Poi_optL', 'unwe_Poi_LCC'))
colnames(proportion_dataframe) <- as.character(num_r)
# row names are different methods, col names are different r to be considered
print(title) # print current simu data type
#mse_dataframe
#proportion_dataframe
#log(mse_dataframe)

####  plot  ##################
# plot for mse
# If you don't want to draw so many methods on the same plot,
# you can comment out certain methods.
mse_plot <- data.frame(r = num_r
                       , we_SWR_optA = mse[1,]
                       , we_SWR_optL = mse[2,]
                       , we_SWR_LCC = mse[3,]

                       # , we_Poi_optA = mse[4,]
                       # , we_Poi_optL = mse[5,]
                       # , we_Poi_LCC = mse[6,]
                       #
                       # , unwe_SWR_optA = mse[7,]
                       # , unwe_SWR_optL = mse[8,]
                       # , unwe_SWR_LCC = mse[9,]
                       #
                       # , unwe_Poi_optA = mse[10,]
                       # , unwe_Poi_optL = mse[11,]
                       # , unwe_Poi_LCC = mse[12,]
                       )

mse_plot<-melt(mse_plot, id.vars = 'r', variable.name="method",
  value.name="mse") # convert data structure to fit ggplot
# mse_plot$mse <- log(mse_plot$mse) # take log(mse)
mse_figure <- ggplot(mse_plot, aes(x=r, y=mse)) +
  geom_line(aes(colour = method)) + geom_point() +
  scale_x_continuous(breaks = num_r) + ggtitle(title) +
  theme_set(theme_bw()) +
  scale_color_manual(values=c("we_SWR_optA" = "red", "we_SWR_optL" = "blue", "we_SWR_LCC" = "green")) #If you draw more than 3 methods, you should define more colors.
#mse_figure
ggsave(paste('draft/mse_', title, '.png'), width = 6, height = 3, dpi = 300)
# plot for proportion
# If you don't want to draw so many methods on the same plot,
# you can comment out certain methods.
proportion_plot <- data.frame(r = num_r
                       , we_SWR_optA = proportion[1,]
                       , we_SWR_optL = proportion[2,]
                       , we_SWR_LCC = proportion[3,]

                       # , we_Poi_optA = proportion[4,]
                       # , we_Poi_optL = proportion[5,]
                       # , we_Poi_LCC = proportion[6,]
                       #
                       # , unwe_SWR_optA = proportion[7,]
                       # , unwe_SWR_optL = proportion[8,]
                       # , unwe_SWR_LCC = proportion[9,]
                       #
                       # , unwe_Poi_optA = proportion[10,]
                       # , unwe_Poi_optL = proportion[11,]
                       # , unwe_Poi_LCC = proportion[12,]
)
proportion_plot<-melt(proportion_plot, id.vars = 'r', variable.name="method",
               value.name="proportion") # convert data structure to fit ggplot
# proportion_plot$proportion <- log(proportion_plot$proportion) # take log(proportion)
proportion_figure <- ggplot(proportion_plot, aes(x=r, y=proportion)) +
  geom_line(aes(colour = method)) + geom_point() +
  scale_x_continuous(breaks = num_r) + ggtitle(title) +
  geom_hline(aes(yintercept=correct_proportion_full), linetype="dashed") +
  theme_set(theme_bw()) +
  scale_color_manual(values=c("we_SWR_optA" = "red", "we_SWR_optL" = "blue", "we_SWR_LCC" = "green"))
proportion_figure
ggsave(paste('draft/proportion_', title, '.png'), width = 6, height = 3, dpi = 300)
