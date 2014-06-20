##' MScPack
##' Description: Tests of SampleVarParms function
##' Author: Rafael Barcellos
##' Last updated 19th June, 2014
##' R 3.1.0

# sourcing DFM(2, 3) ------------------------------------------------------

source("tests/RandGenDfm2Var3.R")

# calling function --------------------------------------------------------

require(Rcpp)
src <- paste(readLines("src/SampleVarParms.cpp"), collapse = "\n")
sourceCpp(code = src)

# sampling VAR parms ------------------------------------------------------

N = 5e3
PhiBar.sim <- array(NA, c(k, k*h, N))
system.time(
for (i in 1:N){
  PhiBar.sim[,, i] <- SampleVarParms(factors, h = h, c0 = 1e3)
}
)

# posterior histograms
par(mfrow = c(2, 6), mar = c(2.1, 2.1, 1.0, 0.5))
for (i in 1:k){
  for (j in 1:(h*k)){
    hist(PhiBar.sim[i, j, ], breaks = 25, col = "darkblue", 
         main = paste("PhiBar[", i, ",", j, "]", sep = ""),
         xlim = range(PhiBar.sim))
    abline(v = PhiBar[i, j], col = "red")
  }
}

# posterior histograms different scales
par(mfrow = c(2, 6), mar = c(2.1, 2.1, 1.0, 0.5))
for (i in 1:k){
  for (j in 1:(h*k)){
    hist(PhiBar.sim[i, j, ], breaks = 25, col = "darkblue", 
         main = paste("PhiBar[", i, ",", j, "]", sep = ""),
         xlim = range(PhiBar.sim[i, j, ], PhiBar[i, j]))
    abline(v = PhiBar[i, j], col = "red")
  }
}

# comparing posterior with OLS --------------------------------------------
library(vars)
PhiBar.mean <- apply(PhiBar.sim, c(1, 2), mean)
mod.var <- VAR(factors, p = 3, type = "none")
PhiBar.ols <- t(sapply(coefficients(mod.var), function(x) x[,1]))

PhiBar; PhiBar.mean; PhiBar.ols
PhiBar.mean/PhiBar.ols
plot(as.vector(PhiBar.mean), as.vector(PhiBar.ols))
abline(coef = c(0, 1))

# Gibbs sampler for factors and VAR parms ---------------------------------

# src <- paste(readLines("src/SampleDynFactors.cpp"), collapse = "\n")
# sourceCpp(code = src)

N <- 1e3
burn <- 2e3
PhiBar.sim <- array(0, c(k, k*h, N))
factors.sim <- array(0, c(k, TT+h, N))
for (j in 1:burn){
  factors.sim[,, 1] <- SampleDynFactors(y, LambdaBar, PhiBar.sim[,, 1], psi)
  PhiBar.sim[,, 1] <- SampleVarParms(t(factors.sim[,, 1]), h = h, c0 = 1e3)
}
system.time(
  for (j in 1:(N-1)){
    factors.sim[,, j+1] <- SampleDynFactors(y, LambdaBar, PhiBar.sim[,, j], psi)
    PhiBar.sim[,, j+1] <- SampleVarParms(t(factors.sim[,, j+1]), h = h, c0 = 1e3)
  }
)

# PhiBar graphs -----------------------------------------------------------

par(mfrow = c(2, 6), mar = c(2.1, 2.1, 1.0, 0.5))
for (i in 1:k){
  for (j in 1:(h*k)){
    hist(PhiBar.sim[i, j, ], breaks = 25, col = "darkblue", 
         main = paste("PhiBar[", i, ",", j, "]", sep = ""),
         xlim = range(PhiBar.sim[i, j, ], PhiBar[i, j]))
    abline(v = PhiBar[i, j], col = "red")
  }
}

# factors graphs ----------------------------------------------------------

factors.mean <- apply(factors.sim, c(1, 2), mean)
factors.upr <- apply(factors.sim, c(1, 2), quantile, probs = 0.95)
factors.lwr <- apply(factors.sim, c(1, 2), quantile, probs = 0.05)

par(mfcol = c(4, 2), mar = c(2.1, 2.1, 2.1, 0.5))
interval <- findInterval(1:(TT+h), c(1:4*(TT+h)/4))
for (i in 1:2){
  for (j in 0:3){
    xcoord <- (1:(TT+h))[interval == j]
    plot(xcoord, factors[interval == j, i], type = "l", lty = 1, 
         col = rgb(0, 0, 0.5, 0.3), lwd = 2, bty = "l",
         main = paste("Factor", i), ylim = range(factors.lwr[i, interval == j], 
                                                 factors.upr[i, interval == j]))
    points(xcoord, factors.mean[i, interval == j], type = "l", lty = 2, col = 2)
    ycoord <- c(factors.lwr[i, interval == j], 
                rev(factors.upr[i, interval == j]))
    xcoord <- c(xcoord, rev(xcoord))
    polygon(xcoord, ycoord, col = rgb(0,0.5,0,0.3), border = rgb(0,0.5,0,0.5))  
  }
}

par(mfcol = c(1, 2), mar = c(2.1, 2.1, 2.1, 0.5))
for (i in 1:2){
  plot(factors[, i] - factors.mean[i, ], type = "p", pch = 20, 
       col = rgb(0, 0, 0.5, 0.3), main = paste("Factor", i), 
       ylim = range(factors[, i] - factors.mean[i, ],
                    factors.lwr[i, ] - factors.mean[i, ], 
                    factors.upr[i, ] - factors.mean[i, ]))
  points(factors.lwr[i, ]-factors.mean[i, ], type = "l", lty = 2, col = 2)
  points(factors.upr[i, ]-factors.mean[i, ], type = "l", lty = 2, col = 2)
}

mean(factors > t(factors.upr))
mean(factors < t(factors.lwr))
plot(factors[,1], factors.mean[1,], pch = 20,col = rgb(0,0,0.5,0.3))
abline(coef = c(0, 1), col = "red")
plot(factors[,2], factors.mean[2,], pch = 20,col = rgb(0,0,0.5,0.3))
abline(coef = c(0, 1), col = "red")

summary(lm(factors[,1]~factors.mean[1,]))$r.sq
summary(lm(factors[,2]~factors.mean[2,]))$r.sq
