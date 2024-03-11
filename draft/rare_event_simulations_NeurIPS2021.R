library(parallel)
library(ggplot2)
library(reshape2)
# source('R/0327_rare_event_functions.R')
# source('R/0403_rare_event_functions.R')
source('R/rare_event_functions_combine.R')
source('draft/generate_rare_logit_data.R')
rm(list = setdiff(ls(), lsf.str())) # remove all variables except functions

rare_event_simu <- function(j){ #this function is used for parLapply
  # source('R/0327_rare_event_functions.R')
  # source('R/0403_rare_event_functions.R')
  source('R/rare_event_functions_combine.R')
  source('draft/generate_rare_logit_data.R')
  set.seed(j)
  beta.1rep <- array(data = NA, dim = c(length(nss), num.method, d))
  FULL.DATA <- generate.rare.data(N, case = case, beta0)
  X.full <- FULL.DATA$X
  Y.full <- FULL.DATA$Y
  beta0 <- FULL.DATA$beta0
  beta.full <- logistic.coef.estimate(X=X.full, Y=Y.full)$beta #beta.full
  for(i in seq_along(nss)){
    n1 <- nss[i]
    beta.1rep[i, 1, ] <- RareLogistic(X.full, Y.full, n0, n1, method = 'UniW')$beta.ssp
    beta.1rep[i, 2, ] <- RareLogistic(X.full, Y.full, n0, n1, method = 'Uni')$beta.ssp
    beta.1rep[i, 3, ] <- RareLogistic(X.full, Y.full, n0, n1, method = 'Weighted', criterion = 'optA')$beta.ssp
    beta.1rep[i, 4, ] <- RareLogistic(X.full, Y.full, n0, n1, method = 'LogOddsCorrection', criterion = 'optA')$beta.ssp
    beta.1rep[i, 5, ] <- RareLogistic(X.full, Y.full, n0, n1, method = 'LogOddsCorrection', criterion = 'LCC')$beta.ssp
    beta.1rep[i, 6, ] <- RareLogistic(X.full, Y.full, n0, n1, method = 'Weighted', criterion = 'optL')$beta.ssp
    beta.1rep[i, 7, ] <- RareLogistic(X.full, Y.full, n0, n1, method = 'LogOddsCorrection', criterion = 'optL')$beta.ssp
  }
  return(list(beta.1rep = beta.1rep,
              beta.full = beta.full))
}

rpt <- 10
N <-  5 * 1e5
beta0 <- c(NA, -rep(1, 6))
case <- 1
beta0 <- generate.rare.data(N, case = case, beta0)$beta0
d <- length(beta0)
n0 <- 200
nss <- c(1000, 2000, 3000, 5000, 10000)#)#
method.all <- c('UniW', 'Uni', 'Weighted+A', 'LogOddsCorrection+A', 'LCC', 'Weighted+L', 'LogOddsCorrection+L')#)
num.method <- length(method.all)
Betas <- array(data = NA, dim = c(rpt, length(nss), num.method, d))
beta.full <- matrix(NA, rpt, d)
# ## single time test
# set.seed(1)
# beta.1rep <- array(data = NA, dim = c(length(nss), num.method, d))
# FULL.DATA <- generate.rare.data(N, case = case, beta0)
# X.full <- FULL.DATA$X
# Y.full <- FULL.DATA$Y
# print(mean(Y.full))
# print(sum(Y.full))
# beta0 <- FULL.DATA$beta0
# for(i in seq_along(nss)){
#   n1 <- nss[i]
#   beta.1rep[i, 1, ] <- logistic.coef.estimate(X=X.full, Y=Y.full)$beta #beta.full
#   beta.1rep[i, 2, ] <- RareLogistic(X.full, Y.full, n0, n1, method = 'UniW')$beta.est
#   beta.1rep[i, 3, ] <- RareLogistic(X.full, Y.full, n0, n1, method = 'Uni')$beta.est
#   beta.1rep[i, 4, ] <- RareLogistic(X.full, Y.full, n0, n1, method = 'Weighted')$beta.est
#   # beta.pilot <- RareLogistic(X.full, Y.full, n0, nss, method = 'Weighted')$beta.pilot
#   beta.1rep[i, 5, ] <- RareLogistic(X.full, Y.full, n0, n1, method = 'LogOddsCorrection')$beta.est
#   beta.1rep[i, 6, ] <- RareLogistic(X.full, Y.full, n0, n1, method = 'LCC')$beta.est
#   # beta.1rep[i, 7, ] <- RareLogistic(X.full, Y.full, n0, nss, method = 'Weighted', criteria = 'optA')$beta.est
#   # beta.1rep[i, 8, ] <- RareLogistic(X.full, Y.full, n0, nss, method = 'LogOddsCorrection', criteria = 'optA')$beta.est
# }


cl <- makeCluster(10) # how many CPU cores are called
t1 <- proc.time()
clusterExport(cl=cl,
              varlist=c('n0', 'case', 'nss', 'beta0', 'num.method', 'd', 'N'),
              envir=environment()) #import environment variables into the function 'subsampling_simu'
results <- parLapply(cl, 1:rpt, rare_event_simu)
t2 <- proc.time()
time.cost <- t2 - t1
for(i in 1:rpt){
  Betas[i,,,] <- results[[i]][[1]]
  beta.full[i,] <- results[[i]][[2]]
}

mse <- array(data = NA, dim = c(num.method, length(nss)))
for (m in 1:num.method){
  for (idn in seq_along(nss)) {
    mse[m, idn] <- mean(apply(sweep(Betas[, idn, m,], 2, beta0) ^ 2 , 1, sum))
    # fname <- paste0("rare_event_simu_output/csv/case",case,"N",N,"method",m,"n",n,".csv")
    # write.csv(Betas[,idn,,m],file=fname,row.names=FALSE)
  }
}
mse.full <- mean(apply(sweep(beta.full, 2, beta0) ^ 2 , 1, sum))
mse <- log(mse)
mse.full <- log(mse.full)
mse_plot <- data.frame(rou = nss/N,
                       uniW = mse[1,],
                       uniLik = mse[2,],
                       optW_A = mse[3,],
                       optLik_A = mse[4,],
                       LCC = mse[5,],
                       optW_L = mse[6,],
                       optLik_L = mse[7,]
                       )
mse_plot<-melt(mse_plot, id.vars = 'rou', variable.name="method",
               value.name="mse") # convert data structure to fit ggplot
mse_figure <- ggplot(mse_plot, aes(x = rou, y = mse)) +
  geom_hline(aes(yintercept = mse.full), linetype = "dashed") +
  geom_line(aes(colour = method, linetype=method)) + geom_point() +
  scale_x_continuous(breaks = nss / N) +
  theme_set(theme_bw()) +
  scale_color_manual(
    values = c(
      "uniW" = "blue",
      "uniLik" = "orange",
      "optW_A" = "green",
      "optLik_A" = "purple",
      "LCC" = "brown",
      "optW_L" = "green",
      "optLik_L" = "purple"
    )) +
  scale_linetype_manual(
    values = c(
      "uniW" = "solid",
      "uniLik" = "solid",
      "optW_A" = "solid",
      "optLik_A" = "solid",
      "LCC" = "solid",
      "optW_L" = "dashed",
      "optLik_L" = "dotted"
    )
  )
mse_figure
print(paste0('total time: ',time.cost[3][[1]],'s'))
# ggsave(paste('rare_event_simu_output/fig/mse_case', case, '.png'), width = 6, height = 3, dpi = 300)
