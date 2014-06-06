##' MScPack
##' Description: Testing if SampleDynFactors are estimating parms correctly
##' Author: Rafael Barcellos
##' Last updated 05 June 2014
##' R 3.0.2

# calling SampleDynFactors ------------------------------------------------

require(Rcpp)
src <- paste(readLines("src/SampleDynFactors.cpp"), collapse = "\n")
sourceCpp(code = src)

# sourcing DFM(0, 1) artificial data --------------------------------------

source("tests/RandGenDfm0Var1.R")

# generating draws from posterior -----------------------------------------

N <- 1e3 # number of draws
factors.sim <- array(NA, c(k, TT+h, N))
system.time(
  for (i in 1:N){
    factors.sim[,, i] <- SampleDynFactors(y, LambdaBar, PhiBar, psi)  
  }
)

# calculating mean and confidence limits ----------------------------------

factors.mean <- apply(factors.sim, c(1, 2), mean)
factors.upr <- apply(factors.sim, c(1, 2), quantile, probs = 0.95)
factors.lwr <- apply(factors.sim, c(1, 2), quantile, probs = 0.05)

par(mfcol = c(1, 2), mar = c(2.1, 2.1, 2.1, 0.5))
for (i in 1:2){
  plot(factors[, i], type = "l", col = rgb(0, 0, 0.5, 0.5), 
       main = paste("Factor", i), bty = "l",
       ylim = range(factors.lwr[i, -1], 
                    factors.upr[i, -1]))
  points(factors.mean[i, -1], type = "l", lty = 2, col = rgb(0.5, 0, 0, 0.5))
  points(factors.lwr[i, -1], type = "l", lty = 2, col = rgb(0.5, 0, 0, 0.5))
  points(factors.upr[i, -1], type = "l", lty = 2, col = rgb(0.5, 0, 0, 0.5))
}

par(mfcol = c(4, 2), mar = c(2.1, 2.1, 2.1, 0.5))
interval <- findInterval(time.id+1, c(1:4*(TT+1+0.01)/4))
for (i in 1:2){
  for (j in 0:3){
    xcoord <- time.id[interval == j]+1
    plot(xcoord, factors[interval == j, i], type = "l", lwd = 2, 
         col = rgb(0, 0, 0.5, 0.5), bty = "l",
         main = if (j == 0){ paste("Factor", i) }, 
         ylim = range(factors.lwr[i, interval == j], 
                      factors.upr[i, interval == j]))
    points(xcoord, factors.mean[i, interval == j], type = "l", lty = 2, col = 2)
    ycoord <- c(factors.lwr[i, interval == j], rev(factors.upr[i, interval == j]))
    xcoord <- c(xcoord, rev(xcoord))
    polygon(xcoord, ycoord, col = rgb(0, 0.5, 0, 0.5), border = rgb(0, 0.5, 0, 0.5))  
  }
}

factors.sd <- apply(factors.sim, c(1,2), sd)
factors.resid <- (t(factors) - factors.mean)/factors.sd
par(mfcol = c(1, 2), mar = c(2.1, 2.1, 2.1, 0.5))
for (i in 1:2){
  plot(factors.resid[i, ], x = time.id, 
       type = "p", pch = 20, col = rgb(0, 0, 0.5, 0.3), 
       main = paste("Factor", i), bty = "l",
       ylim = range(factors.resid[i, ],
                    (factors.lwr[i, ] - factors.mean[i, ])/factors.sd[i,], 
                    (factors.upr[i, ] - factors.mean[i, ])/factors.sd[i,]))
  points((factors.lwr[i, ] - factors.mean[i, ])/factors.sd[i, ], 
         type = "l", lty = 2, col = 2)
  points((factors.upr[i, ] - factors.mean[i, ])/factors.sd[i, ], 
         type = "l", lty = 2, col = 2)
  abline(h = qnorm(c(0.05, 0.95)), lty = 3)
}

mean(factors > t(factors.upr))
mean(factors < t(factors.lwr))

plot(factors[, 1], factors.mean[1, ], bty = "l",
     pch = 20, col = rgb(0, 0, 0.5, 0.3))
abline(coef = c(0, 1), col = "red")
plot(factors[, 2], factors.mean[2, ], bty = "l",
     pch = 20, col = rgb(0, 0, 0.5, 0.3))
abline(coef = c(0, 1), col = "red")

summary(lm(factors[,1]~factors.mean[1, ]))$r.sq
summary(lm(factors[,2]~factors.mean[2, ]))$r.sq

used.packs <- c("Rcpp", "RcppArmadillo", "truncnorm")
sapply(used.packs, function(x) toBibtex(citation(x)))



# sourcing DFM(2, 3) ------------------------------------------------------
rm(list = ls()[ls() != "SampleDynFactors"])
source("tests/RandGenDfm2Var3.R")

# generating draws from posterior -----------------------------------------

N <- 2e3
factors.sim <- array(NA, c(k, TT+h, N))
system.time(
  for (i in 1:N){
    factors.sim[,, i] <- SampleDynFactors(y, LambdaBar, PhiBar, psi)  
  }
)

# calculating mean and confidence limits ----------------------------------

factors.mean <- apply(factors.sim, c(1, 2), mean)
factors.upr <- apply(factors.sim, c(1, 2), quantile, probs = 0.95)
factors.lwr <- apply(factors.sim, c(1, 2), quantile, probs = 0.05)

par(mfcol = c(1, 2), mar = c(2.1, 2.1, 2.1, 0.5))
for (i in 1:2){
  plot(factors[, i], type = "l", col = rgb(0, 0, 0.5, 0.5), 
       main = paste("Factor", i), bty = "l",
       ylim = range(factors.lwr[i, -1], 
                    factors.upr[i, -1]))
  points(factors.mean[i, -1], type = "l", lty = 2, col = rgb(0.5, 0, 0, 0.5))
  points(factors.lwr[i, -1], type = "l", lty = 2, col = rgb(0.5, 0, 0, 0.5))
  points(factors.upr[i, -1], type = "l", lty = 2, col = rgb(0.5, 0, 0, 0.5))
}

par(mfcol = c(4, 2), mar = c(2.1, 2.1, 2.1, 0.5))
interval <- findInterval(time.id+1, c(1:4*(TT+1+0.01)/4))
for (i in 1:2){
  for (j in 0:3){
    xcoord <- time.id[interval == j]+1
    plot(xcoord, factors[interval == j, i], type = "l", lwd = 2, 
         col = rgb(0, 0, 0.5, 0.5), bty = "l",
         main = if (j == 0){ paste("Factor", i) }, 
         ylim = range(factors.lwr[i, interval == j], 
                      factors.upr[i, interval == j]))
    points(xcoord, factors.mean[i, interval == j], type = "l", lty = 2, col = 2)
    ycoord <- c(factors.lwr[i, interval == j], rev(factors.upr[i, interval == j]))
    xcoord <- c(xcoord, rev(xcoord))
    polygon(xcoord, ycoord, col = rgb(0, 0.5, 0, 0.5), border = rgb(0, 0.5, 0, 0.5))  
  }
}

factors.sd <- apply(factors.sim, c(1,2), sd)
factors.resid <- (t(factors) - factors.mean)/factors.sd
par(mfcol = c(1, 2), mar = c(2.1, 2.1, 2.1, 0.5))
for (i in 1:2){
  plot(factors.resid[i, ], x = time.id, 
       type = "p", pch = 20, col = rgb(0, 0, 0.5, 0.3), 
       main = paste("Factor", i), bty = "l",
       ylim = range(factors.resid[i, ],
                    (factors.lwr[i, ] - factors.mean[i, ])/factors.sd[i,], 
                    (factors.upr[i, ] - factors.mean[i, ])/factors.sd[i,]))
  points((factors.lwr[i, ] - factors.mean[i, ])/factors.sd[i, ], 
         type = "l", lty = 2, col = 2)
  points((factors.upr[i, ] - factors.mean[i, ])/factors.sd[i, ], 
         type = "l", lty = 2, col = 2)
  abline(h = qnorm(c(0.05, 0.95)), lty = 3)
}

mean(factors > t(factors.upr))
mean(factors < t(factors.lwr))

plot(factors[, 1], factors.mean[1, ], bty = "l",
     pch = 20, col = rgb(0, 0, 0.5, 0.3))
abline(coef = c(0, 1), col = "red")
plot(factors[, 2], factors.mean[2, ], bty = "l",
     pch = 20, col = rgb(0, 0, 0.5, 0.3))
abline(coef = c(0, 1), col = "red")

summary(lm(factors[, 1] ~ factors.mean[1, ]))$r.sq
summary(lm(factors[, 2] ~ factors.mean[2, ]))$r.sq
