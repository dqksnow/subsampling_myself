proptional_ssp <- function(N, K, y){
  prob <- 1/((K+1) * as.data.frame(table(y))$Freq)
  return(prob[y+1])
}
###############################################################################
# psi_multi <- function(X, p.plt, P.plt, n.plt, K, d, y.matrix.plt) {
#   num <- min(K,2)
#   comb <- combn(1:K, num)
#   pi_s2_s <- p.plt / n.plt
#   r_pi <- 1 - p.plt
#
#   psii <- matrix(NA, nrow = d^2, ncol = K^2)
#   psii[, (comb[1,]-1)*K+comb[num,]] <-
#     psii[,(comb[num,]-1)*K+comb[1,]] <-
#     t(X) %*% ((as.matrix((y.matrix.plt-P.plt)[,comb[1,]] *
#                            (y.matrix.plt-P.plt)[,comb[num,]])[,rep(1:dim(comb)[2], each=d)]*
#                  (r_pi/(pi_s2_s^2)))*X[,rep(1:d, dim(comb)[2])])
#   psii[, seq(1,(K^2), by = (K+1))] <-
#     t(X) %*% ((((y.matrix.plt-P.plt)^2)[,rep(1:K, each=d)]*
#                  (r_pi/(pi_s2_s^2)))*X[,rep(1:d, K)])
#   psii <- t(matrix(psii, nrow = d))
#   psi <- matrix(psii, nrow = K*d)
#   psi <- t(psi[, c(rep(seq(1,K*d,K), K) + rep(0:(K-1), each = d))])
#   psi
# }
###############################################################################
random.index <- function (N, n, p=NULL) {
  ifelse(is.null(p),
         index <- sample(N, n, replace = TRUE),
         index <- sample(N, n, replace = TRUE, prob = p))
  return(as.vector(index))
}
softmax_coef_estimate <- function(X, y, weights){
  fit <- nnet::multinom(y ~ X - 1, weights = weights)
  t(coef(fit))
}
###############################################################################
pbeta_multi <- function(X, beta){
  exp(X %*% beta - matrixStats::rowLogSumExps(cbind(0, X %*% beta)))
}
###########################
MN_multi <- function(X, P, K, d, p, N, n){
  X_Kd <- X[, rep(1:d, K)]
  diagp <- P[, rep(1:K, each = d)]
  spi <- sqrt(p)
  diagpX_pd <- (diagp * X_Kd)/spi
  MN_p1 <- t(X_Kd/spi) %*% diagpX_pd
  MN_p2 <- t(diagpX_pd) %*% (diagpX_pd)
  MN <- MN_p1 * (diag(1, K) %x% matrix(1, nrow = d, ncol = d)) - MN_p2
  MN <- MN / (n * N)
  return(MN)
}
Psi_multi <- function(X, Y.matrix, P, K, d, p, N, n){
  X_Kd <- X[, rep(1:d, K)]
  psi <- Y.matrix[, rep(1:K, each = d)] - P[, rep(1:K, each = d)]
  spi <- p
  psiX <- psi * X_Kd / spi
  # psiX <- psi * X_Kd / p
  PsiX <- t(psiX) %*% (psiX) / (n^2 * N^2)
  return(PsiX)
}
###############################################################################
library(MASS)
rm(list = setdiff(ls(), lsf.str()))
d <- 3 # dim of covariates
K <- 5 # K + 1 classes
G <- rbind(rep(-1/(K+1), K), diag(K) - 1/(K+1)) %x% diag(d)
N <- 1e4
beta.true <- 0.2 * matrix(-1, d, K)
beta.true.sum <- cbind(rep(1, d), beta.true)
# beta.true <- matrix(c(rep(0,d), rep(-1, K*d)), d, K+1) # baseline constraint
# set.seed(1)
mu <- rep(0, d)
sigma <- matrix(0.5, nrow = d, ncol = d)
diag(sigma) <- rep(1, d)
X <- mvrnorm(N, mu, sigma)

# prob.K <- exp(X %*% beta.true)
# prob.K <- prob.K / (1 + rowSums(prob.K)) # can I use pbeta_multi?
# prob.Kplus1 <- cbind(1 - rowSums(prob.K), prob.K)
# y <- apply(prob.Kplus1, 1, function(row) sample(0:K, size = 1, prob = row))
# table(y)
# unique(y)
prob.Kplus2 <- exp( X %*% beta.true.sum)
prob.Kplus2 <- prob.Kplus2 / rowSums(prob.Kplus2)
y <- apply(prob.Kplus2, 1, function(row) sample(0:K, size = 1, prob = row))
table(y)
fit.full <- nnet::multinom(y ~ X - 1)
beta.full <- G %*% as.vector(t(coef(fit.full)))
beta.full <- matrix(beta.full, nrow = d)
mse.full <- sum((beta.full - beta.true.sum)^2)

y.matrix <- matrix(0, nrow = N, ncol = K)
y.matrix[cbind(seq_along(y), y)] <- 1
#
n.plt <- 300
n.ssp <- 3000
index.plt <- random.index(N, n.plt)
p.plt <- rep(1 / N, n.plt) # pilot sampling probability
x.plt <- X[index.plt,]
y.plt <- y[index.plt]
y.matrix.plt <- y.matrix[index.plt, ]
beta.plt <- softmax_coef_estimate(x.plt, y.plt, weights = 1 / p.plt) # weigths=1
mse.plt <- sum((matrix(G %*% as.vector(t(beta.plt)), nrow = d) - beta.true.sum)^2)
P.step1 <- pbeta_multi(X, beta.plt) # N*K, i-th row = Prob(y_ik=1) for i= 1 to K
P.plt <- P.step1[index.plt, ]
MN.plt <- MN_multi(X = x.plt, P = P.plt, K, d, p = p.plt, N, n.plt)

# psi.plt <- psi_multi(x.plt, p.plt,
#                   P.plt, n.plt, K, d, y.matrix.plt)
Psi.plt <- Psi_multi(X = x.plt, Y.matrix = y.matrix.plt, P = P.plt, K = K,
                     d = d, p = p.plt, N = N, n = n.plt)
sixi <- (y.matrix-P.step1)[,rep(seq(K), each = d)] * X[,rep(seq(d), K)]
G <- rbind(rep(-1/(K+1), K), diag(K) - 1/(K+1)) %x% diag(d)

constraint = "baseline" # summation
if (constraint == "baseline"){
  pi_num <- sqrt(colSums((solve(MN.plt, t(sixi))^2)))
} else if (constraint == "summation"){
  pi_num <- sqrt(colSums(((G %*% solve(MN.plt, t(sixi)))^2)))
}
ossp <- pi_num/sum(pi_num)
index.ssp <- random.index(N, n.ssp, p = ossp)
p.ssp <- ossp[index.ssp]
x.ssp <- X[index.ssp,]
y.ssp <- y[index.ssp]
y.matrix.ssp <- y.matrix[index.ssp, ]
beta.ssp <- softmax_coef_estimate(x.ssp, y.ssp, weights = 1 / p.ssp)
mse.ssp <- sum((matrix(G %*% as.vector(t(beta.ssp)), nrow = d) - beta.true.sum)^2)
P.ssp <- pbeta_multi(x.ssp, beta.ssp)
# MN.ssp <- MN_multi(x.ssp, P.ssp, K, d, pi = p.ssp, n.ssp)
MN.ssp <- MN_multi(X = x.ssp, P = P.ssp, K, d, p = p.ssp, N, n.ssp)
# psi.ssp <- psi_multi(x.ssp, p.ssp, P.ssp, n.ssp, K, d, y.matrix.ssp)
Psi.ssp <- Psi_multi(X = x.ssp, Y.matrix = y.matrix.ssp, P = P.ssp, K = K,
                     d = d, p = p.ssp, N = N, n = n.ssp)

MN.plt <- n.plt * MN.plt
MN.ssp <- n.ssp * MN.ssp
Psi.plt <- n.plt^2 * Psi.plt
Psi.ssp <- n.ssp^2 * Psi.ssp
solveMN <- solve(MN.plt + MN.ssp)
beta.cmb <- solveMN %*% MN.plt %*% c(beta.plt) + solveMN %*% MN.ssp %*% c(beta.ssp)
beta.cmb <- matrix(beta.cmb, ncol = K)
mse.cmb <- sum((matrix(G %*% as.vector(t(beta.cmb)), nrow = d) - beta.true.sum)^2)

var.cmb <- solveMN %*% (Psi.plt + Psi.ssp) %*% solveMN
mse.cmb.the <- sum(diag(var.cmb))
var.plt <- solve(MN.plt) %*% (Psi.plt) %*% solve(MN.plt)
mse.plt.the <- sum(diag(var.plt))
var.ssp <- solve(MN.ssp) %*% (Psi.ssp) %*% solve(MN.ssp)
mse.ssp.the <- sum(diag(var.ssp))
###############################################################################







#####
num <- min(K,2)
comb <- combn(1:K, num)
pi_s2_s <- p.plt/n.plt
r_pi <- 1 - p.plt

psii <- matrix(NA, nrow = d^2, ncol = K^2)
psii[, (comb[1,]-1)*K+comb[num,]] <-
  psii[,(comb[num,]-1)*K+comb[1,]] <-
  t(x.plt) %*% ((as.matrix((y.matrix.plt-P.plt)[,comb[1,]] *
                         (y.matrix.plt-P.plt)[,comb[num,]])[,rep(1:dim(comb)[2], each=d)] *
               (r_pi/(pi_s2_s^2))) *
                 x.plt[,rep(1:d, dim(comb)[2])])

psii[, seq(1,(K^2), by = (K+1))] <-
  t(x.plt) %*% ((((y.matrix.plt-P.plt)^2)[,rep(1:K, each=d)]*
               (r_pi/(pi_s2_s^2)))
               * x.plt[,rep(1:d, K)])
psii <- t(matrix(psii, nrow = d))
psi <- matrix(psii, nrow = K*d)
psi <- t(psi[, c(rep(seq(1,K*d,K), K) + rep(0:(K-1), each = d))])


#############
compute_phi_i <- function(Xi) {
  diag(Xi) - Xi %*% t(Xi)
}
compute_outer_i <- function(Xi) {
  Xi %*% t(Xi)
}


A <- array(apply(P.plt,
                 1,
                 function(row) {
                   diag(row) - row %*% t(row)
                 }),
           dim = c(n.plt, K, K)
           )
B <- array(apply(x.plt,
                 1,
                 function(row) {
                   row %*% t(row)
                 }),
           dim = c(n.plt, d, d)
           )

# C <- array(mapply(kronecker, A, B, SIMPLIFy = 'array'),
#            dim = c(n.plt, K*d, K*d))

C <- array(0, dim = c(n.plt, K*d, K*d))
for (i in 1:n.plt) {
  C[i,,] <- kronecker(A[i,,], B[i,,])
}

C <- matrix(0, nrow = K*d, ncol = K*d)
for (i in 1:n.plt) {
  C <- C + w[i] * (X[i, ] %*% t(X[i, ])) %x% (S[i, ] %*% t(S[i, ]))
}

A <- alply(A, 1, function(a) matrix(a, ncol=K, nrow=K))
B <- alply(B, 1, function(a) matrix(a, ncol=d, nrow=d))
mapply(function(a,b) a%*%b, A, B, SIMPLIFy=FALSE)

M.N <- Map(kronecker, A, B)




M.N <- Reduce("+", Map(kronecker, A_list, Q_list))

A_list <- lapply(split(P.plt, row(P.plt)), compute_phi_i)
Q_list <- lapply(split(x.plt, row(x.plt)), compute_outer_i)
B <- Reduce("+", Map(kronecker, A_list, Q_list))
# A_matrix <- do.call(rbind, A_list)




test1 <- c(0, 1, 1)
test2 <- c(1, 1)
matrixStats::logSumExp(test1)
matrixStats::logSumExp(test2) + log(1)
log(exp(1)+exp(1))
log(exp(0)+exp(1)+exp(1))

MN1 <- MN_multi(x.plt, pbeta_pilot[pilot_indx,],
                k, d, pi = pinv_pilot)
psi1 <- psi_multi(X[pilot_indx,], pilot_ssp[pilot_indx],
                  pbeta_pilot[pilot_indx,], r0, k, d, I[pilot_indx,])
#
# k <- length(unique(y))-1
# N <- nrow(X)
# d <- ncol(X)
# I <- matrix(0, nrow = N, ncol = k)
# I[cbind(seq_along(y), y)] <- 1
# pilot_ssp <- proptional_ssp(N, k, y)
# pilot_indx <- swr_indx(N, r0, pilot_ssp)


###############################################################################
pinv_pilot <- 1/pilot_ssp[pilot_indx]

beta_pilot <- softmax_coef_estimate(X, y, pinv_pilot, pilot_indx)
pbeta_pilot <- pbeta_multi(X, beta_pilot)
###########################
MN_multi <- function(X, p, k, d, pi){
  X_kd <- X[, rep(1:d, k)]
  diagp <- p[, rep(1:k, each = d)]
  spi <- sqrt(pi)
  diagpX_pd <- (diagp * X_kd)/spi
  MN_p1 <- t(X_kd/spi) %*% diagpX_pd

  MN_p2 <- t(diagpX_pd) %*% (diagpX_pd)
  MN <- MN_p1 * (diag(1, k) %x% matrix(1, nrow = d, ncol = d)) - MN_p2
  MNjnm
}
###########################
MN1 <- MN_multi(X[pilot_indx,], pbeta_pilot[pilot_indx,],
                k, d, pi = pinv_pilot)
psi1 <- psi_multi(X[pilot_indx,], pilot_ssp[pilot_indx],
                  pbeta_pilot[pilot_indx,], r0, k, d, I[pilot_indx,])

sixi <- (I-pbeta_pilot)[,rep(seq(1:k), each = ncol(X))] * # why N * kd?
  X[,rep(seq(ncol(X)), k)]
pi_num <- sqrt(colSums((solve(MN1, t(sixi))^2)))


ossp <- pi_num/sum(pi_num)
second_indx <- swr_indx(N, r, ossp)

pinv_s2 <- 1/ossp[second_indx]
invisible(capture.output(
  beta_s2 <- softmax_coef_estimate(X, y, pinv_s2, second_indx)))
pbeta_s2  <- pbeta_multi(X[second_indx,], beta_s2)

MN2 <- MN_multi(X[second_indx,], pbeta_pilot[second_indx,],
                k, d, pi = pinv_s2)
psi2 <- psi_multi(X[second_indx,], pilot_ssp[second_indx],
                  pbeta_pilot[second_indx,], r, k, d, I[second_indx,])

beta_cmb <- solve(MN1 + MN2, (MN1 %*% c(beta_pilot) + MN2 %*% c(beta_s2)))
beta_cmb <- matrix(beta_cmb, ncol = k)
var <- (solve(MN1 + MN2, psi1 + psi2) %*% solve(MN1 + MN2))
trvar <- diag(var)
