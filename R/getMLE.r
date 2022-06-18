getMLE <- function(x, y, w=1) {
    beta <- rep(0, d)
dig <- diag(0.00001, d)
    loop  <- 1
    Loop  <- 100
    msg <- "NA"
    while (loop <= Loop) {
        pr <- c(1 - 1 / (1 + exp(x %*% beta)))
        H <- t(x) %*% (pr * (1 - pr) * w * x)
        S <- colSums((y - pr) * w * x)
        tryCatch(
            {shs <- NA
             shs <- solve(H + dig, S) }, 
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
        if (loop == Loop) {
            warning("Maximum iteration reached")
            ## beta <- NA
        }
        loop  <- loop + 1
    }
    list(par=beta, message=msg, iter=loop)
}

## ## r1 <- r1+r; r2=0

## AlgTwoStp <- function(r1=r1, r2=r2) {
##         ## mVc
##         PI.mVc <- sqrt((Y - P.prop)^2 * rowSums(X^2))
##         PI.mVc <- PI.mVc / sum(PI.mVc)
##         idx.mVc <- sample(1:n, r2, T, PI.mVc)
##         x.mVc <- X[c(idx.mVc, idx.prop),]
##         y.mVc <- Y[c(idx.mVc, idx.prop)]
##         pinv.mVc <- c(1 / PI.mVc[idx.mVc], pinv.prop)
##         fit.mVc <- getMLE(x=x.mVc, y=y.mVc, w=pinv.mVc)
##         ## mMSE
##         W.prop <- solve(t(x.prop) %*% (x.prop * w.prop * pinv.prop))
##         PI.mMSE <- sqrt((Y - P.prop)^2 * rowSums((X%*%W.prop)^2))
##         PI.mMSE <- PI.mMSE / sum(PI.mMSE)
##         idx.mMSE <- sample(1:n, r2, T, PI.mMSE)
##         x.mMSE <- X[c(idx.mMSE, idx.prop),]
##         y.mMSE <- Y[c(idx.mMSE, idx.prop)]
##         pinv.mMSE <- c(1 / PI.mMSE[idx.mMSE], pinv.prop)
##         fit.mMSE <- getMLE(x=x.mMSE, y=y.mMSE, w=pinv.mMSE)

##         ## mVcUW
##         x.mVcuw <- X[c(idx.mVc),]
##         y.mVcuw <- Y[c(idx.mVc)]
##         fit.mVcuw <- getMLE(x=x.mVcuw, y=y.mVcuw, w=1)
## beta.mVcuw <- fit.mVcuw$par+beta.prop
## p.mVcuw  <- 1 - 1 / (1 + exp(c(x.mVcuw %*% beta.mVcuw)))
## w.mVcuw <- p.mVcuw * (1 - p.mVcuw)
## ldd.mVcuw <- t(x.mVcuw) %*% (x.mVcuw * w.mVcuw)
## beta.mVcuwcb <- solve(ldd.prop+ldd.mVcuw,
##       ldd.prop %*% beta.prop + ldd.mVcuw %*% beta.mVcuw)
        
##         ## mMSEUW
##         x.mMSEuw <- X[c(idx.mMSE),]
##         y.mMSEuw <- Y[c(idx.mMSE)]
##         fit.mMSEuw <- getMLE(x=x.mMSEuw, y=y.mMSEuw, w=1)
## beta.mMSEuw <- fit.mMSEuw$par+beta.prop
## p.mMSEuw  <- 1 - 1 / (1 + exp(c(x.mMSEuw %*% beta.mMSEuw)))
## w.mMSEuw <- p.mMSEuw * (1 - p.mMSEuw)
## ldd.mMSEuw <- t(x.mMSEuw) %*% (x.mMSEuw * w.mMSEuw)
## beta.mMSEuwcb <- solve(ldd.prop+ldd.mMSEuw,
##       ldd.prop %*% beta.prop + ldd.mMSEuw %*% beta.mMSEuw)

##         ## LCC
##         PI.LCC <- abs(Y - P.prop)
##         PI.LCC <- PI.LCC / sum(PI.LCC)
##         idx.LCC <- sample(1:n, r2, T, PI.LCC)
##         x.LCC <- X[c(idx.LCC),]
##         y.LCC <- Y[c(idx.LCC)]
##         pinv.LCC <- 1
##         fit.LCC <- getMLE(x=x.LCC, y=y.LCC, w=pinv.LCC)
## beta.LCC <- fit.LCC$par+beta.prop
## p.LCC  <- 1 - 1 / (1 + exp(c(x.LCC %*% beta.LCC)))
## w.LCC <- p.LCC * (1 - p.LCC)
## ldd.LCC <- t(x.LCC) %*% (x.LCC * w.LCC)
## beta.LCCcb <- solve(ldd.prop+ldd.LCC,
##       ldd.prop %*% beta.prop + ldd.LCC %*% beta.LCC)



##         opt <- cbind(fit.mMSE$par, fit.mVc$par, beta.LCCcb,
##                      beta.mMSEuwcb, beta.mVcuwcb)
##         msg <- c(fit.prop$message, fit.mMSE$message, fit.mVc$message,
##                  fit.LCC$message, fit.mMSEuw$message, fit.mVcuw$message)
##         return(list(opt=opt, msg=msg))
## }
