#' MScPack
#' Description: Testing SampleDfmLoads function
#' Author: Rafael Barcellos
#' Last updated: 18th June, 2014
#' R 3.1.0

library(MScPack)

# sourcing DFM(2, 3) ------------------------------------------------------

source("tests/RandGenDfm2Var3.R")

# calling functions -------------------------------------------------------

require(Rcpp)
src <- paste(readLines("src/SampleDynFactors.cpp"), collapse = "\n")
sourceCpp(code = src)
src <- paste(readLines("src/SampleVarParms.cpp"), collapse = "\n")
sourceCpp(code = src)
src <- paste(readLines("src/SampleDfmLoads.cpp"), collapse = "\n")
sourceCpp(code = src)

# Testing SampleDfmLoads alone --------------------------------------------

N <- 5e3
LambdaBar.sim <- array(NA, c(q, k*(s+1), N))
system.time(
  for (i in 1:N){
    LambdaBar.sim[,, i] <- SampleDfmLoads(y, factors, psi, s, 1e1)  
  }
)

# analyzing posterior draws -----------------------------------------------

# analytic solution
x <- NULL
TpH <- nrow(factors)
for (i in 0:s){
  x <- cbind(x, factors[(h-i+1):(TpH-i), ])
}

Xlambda <- kronecker(diag(1, q), x)

ystar <- as.vector(y)

inv.L1 <- t(Xlambda) %*% kronecker(diag(1/psi), diag(1, TT)) %*% Xlambda + 1/1e1
L1 <- solve(inv.L1)
l1 <- L1 %*% (t(Xlambda) %*% kronecker(diag(1/psi), diag(1, TT)) %*% ystar)

LambdaBar.exact.mean <- t(array(l1, c(k*(s+1), q)))
LambdaBar.exact.sd <- sqrt(t(array(diag(L1), c(k*(s+1), q))))

# comparison plot
par(mfrow = dim(LambdaBar), mar = c(2.1, 2.1, 1.0, 0.5))
for (i in 1:q){
  for (j in 1:(k*(s+1))){
    hist(LambdaBar.sim[i, j, ], breaks = 25, col = "darkblue", 
         main = "", freq = FALSE)#paste("LambdaBar[", i, ",", j, "]", sep = ""))
    abline(v = LambdaBar[i, j], col = "red")
    curve(dnorm(x, LambdaBar.exact.mean[i, j], LambdaBar.exact.sd[i, j]),
          add = T, col = "green")
  }
}

par(mfrow = dim(LambdaBar), mar = c(2.1, 2.1, 1.0, 0.5))
for (i in 1:q){
  for (j in 1:(k*(s+1))){
    dLambda <- density(LambdaBar.sim[i, j, ], n = 1000)
    plot(dLambda, main = "", ylab = "", xlab = "", col = "darkblue")
    abline(v = LambdaBar[i, j], col = "red")
    curve(dnorm(x, LambdaBar.exact.mean[i, j], LambdaBar.exact.sd[i, j]),
          add = T, col = "green")
  }
}

# how 5000 draws from a normal behave
par(mfrow = dim(LambdaBar), mar = c(2.1, 2.1, 1.0, 0.5))
for (i in 1:q){
  for (j in 1:(k*(s+1))){
    plot(density(rnorm(5e3)), main = "", ylab = "", xlab = "")
    curve(dnorm(x, 0, 1), add = T, col = "green")
  }
}

# comparing posterior draws with ols --------------------------------------

# posterior mean
LambdaBar.mean <- apply(LambdaBar.sim, c(1, 2), mean)

# OLS regression
reg <- lm(ystar ~ Xlambda-1)
lambda.est <- coef(reg)

# scatter plot between OLS and posterior estimates
par(mfrow = c(1, 1), mar = c(4.1, 2.1, 2.1, 0.5))
plot(lambda.est, as.vector(t(LambdaBar.mean)), pch = 20, 
     col = rgb(0, 0, 0.5, 0.3), bty = "l", 
     ylim = range(LambdaBar, LambdaBar.mean),
     xlab = "OLS estimates")
points(lambda.est, as.vector(t(LambdaBar)), pch = 15,
       col = rgb(0, 0.5, 0, 0.3))
abline(coef = c(0, 1), col = 2)
legend("topleft", legend = c("Posterior estimate", "True values", 
                             "identity between OLS and posterior"), 
       lty = c(NA, NA, 1), col = c(rgb(0, 0, 0.5, 0.3), rgb(0, 0.5, 0, 0.3), "red"),
       pch = c(20, 15, NA), pt.cex = 1.5, bty = "n")
grid()
