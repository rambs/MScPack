
library(Rcpp)
library(RcppArmadillo)


src <- paste(readLines("src/SampleDynFactors.cpp"), collapse = "\n")
sourceCpp(code = src)


# dados artificiais ####
TT = 500
psi = c(0.02, 0.19, 0.36, 0.02, 0.02, 0.19, 0.19, 0.36, 0.36)
q = length(psi)
k = 2

PhiBar <- matrix(c(1, -0.1, 0, 0.7), 2, 2)
set.seed(7623)
factors <- array(rnorm(TT*2), c(TT, 2))
for (i in 1:TT){
  if (i == 1){
    next
  }
  factors[i, ] <- PhiBar %*% factors[i-1, ] + factors[i, ]
}
plot(factors)
par(mfrow = c(1, 2), mar = c(2.1, 2.1, 2.1, 0.2))
apply(factors, 2, plot, type = "l")

set.seed(4091)
LambdaBar  <- matrix(runif(2*q, -0.1, 0.1), q, 2)

set.seed(3708)
y <- t(array(rnorm(q*TT, 0, sqrt(psi)), c(q, TT)))
y <- factors %*% t(LambdaBar) + y
par(mfrow = c(3, 3))
apply(y, 2, plot, type = "l")


N <- 1e4
factors.sim <- array(NA, c(2, TT, N))
system.time(
for (i in 1:N){
  factors.sim[,, i] <- SampleDynFactors(y, LambdaBar, PhiBar, psi)  
}
)
factors.mean <- apply(factors.sim, c(1, 2), mean)
factors.upr <- apply(factors.sim, c(1, 2), quantile, probs = 0.95)
factors.lwr <- apply(factors.sim, c(1, 2), quantile, probs = 0.05)

par(mfcol = c(1, 2), mar = c(2.1, 2.1, 2.1, 0.5))
for (i in 1:2){
  plot(factors[, i], type = "p", pch = 20, col = rgb(0, 0, 0.5, 0.3), 
       main = paste("Factor", i), ylim = range(factors.lwr[i,], factors.upr[i,]))
  points(factors.mean[i, ], type = "l", lty = 2, col = 2)
  points(factors.lwr[i, ], type = "l", lty = 2, col = 2)
  points(factors.upr[i, ], type = "l", lty = 2, col = 2)
}


par(mfcol = c(1, 2), mar = c(2.1, 2.1, 2.1, 0.5))
for (i in 1:2){
  plot(factors[, i] - factors.mean[i, ], type = "p", pch = 20, col = rgb(0, 0, 0.5, 0.3), 
       main = paste("Factor", i), ylim = range(factors[, i] - factors.mean[i,],
                                               factors.lwr[i,]-factors.mean[i,], 
                                               factors.upr[i,]-factors.mean[i,]))
  points(factors.lwr[i, ]-factors.mean[i, ], type = "l", lty = 2, col = 2)
  points(factors.upr[i, ]-factors.mean[i, ], type = "l", lty = 2, col = 2)
}

mean(factors > t(factors.upr))
mean(factors < t(factors.lwr))
plot(factors[,1], factors.mean[1,])
abline(coef = c(0, 1), col = "red")
plot(factors[,2], factors.mean[2,])
abline(coef = c(0, 1), col = "red")

summary(lm(factors[,1]~factors.mean[1,]))
summary(lm(factors[,2]~factors.mean[2,]))


# dados artificiais VAR(3) - DFM(2)####
TT <- 2500
psi <- c(0.02, 0.19, 0.36, 0.02, 0.02, 0.19, 0.19, 0.36, 0.36)
q <- length(psi)
k <- 2
h <- 3
r <- h*k

Phi1 <- matrix(c(0, -0.1, 0, -0.7), 2, 2) + diag(rep(1, k))
Phi2 <- matrix(c(0.02, -0.08, 0.4, -0.2), 2, 2)
Phi3 <- matrix(c(-0.06, 0.07, -0.6, 0.35), 2, 2)
PhiBar <- cbind(Phi1, Phi2, Phi3)

set.seed(7623)
factors <- array(rnorm(TT*2), c(TT, 2))
phi = rep(0, r)
for (i in 1:TT){
  factors[i, ] <- PhiBar %*% phi + factors[i, ]
  phi <- c(factors[i, ], phi[1:(r-k)])
}

plot(factors)
par(mfrow = c(1, 2), mar = c(2.1, 2.1, 2.1, 0.2))
apply(factors, 2, plot, type = "l")

set.seed(4091)
s <- 2
LambdaBar  <- matrix(runif(k*q*(s+1), -0.1, 0.1), q, k*(s+1))

set.seed(3708)
y <- t(array(rnorm(q*TT, 0, sqrt(psi)), c(q, TT)))
f.star <- rep(0, k*(s+1))
for (i in 1:TT){
  f.star <- c(factors[i, ], f.star[1:(k*s)])
  y[i, ] <- LambdaBar %*% f.star + y[i, ]
}
par(mfrow = c(3, 3))
apply(y, 2, plot, type = "l")

# calling function ####
src <- paste(readLines("src/SampleDynFactors.cpp"), collapse = "\n")
sourceCpp(code = src)
SampleDynFactors(y, LambdaBar, PhiBar, psi)  

N <- 1e3
factors.sim <- array(NA, c(2, TT, N))
system.time(
  for (i in 1:N){
    factors.sim[,, i] <- SampleDynFactors(y, LambdaBar, PhiBar, psi)  
  }
)
factors.mean <- apply(factors.sim, c(1, 2), mean)
factors.upr <- apply(factors.sim, c(1, 2), quantile, probs = 0.95)
factors.lwr <- apply(factors.sim, c(1, 2), quantile, probs = 0.05)

par(mfcol = c(4, 2), mar = c(2.1, 2.1, 2.1, 0.5))
interval <- findInterval(1:TT, c(1:4*(2500+1)/4))
for (i in 1:2){
  for (j in 0:3){
    xcoord <- (1:TT)[interval == j]
    plot(xcoord, factors[interval == j, i], type = "l", lty = 1, col = rgb(0, 0, 0.5, 0.3), 
         main = paste("Factor", i), ylim = range(factors.lwr[i, interval == j], 
                                                 factors.upr[i, interval == j]))
    points(xcoord, factors.mean[i, interval == j], type = "l", lty = 2, col = 2)
    ycoord <- c(factors.lwr[i, interval == j], rev(factors.upr[i, interval == j]))
    xcoord <- c(xcoord, rev(xcoord))
#     polygon(xcoord, ycoord, col = rgb(0,0.5,0,0.3), border = rgb(0,0.5,0,0.5))  
  }
  
#   points(factors.lwr[i, ], type = "l", lty = 2, col = 3)
#   points(factors.upr[i, ], type = "l", lty = 2, col = 3)
}


par(mfcol = c(1, 2), mar = c(2.1, 2.1, 2.1, 0.5))
for (i in 1:2){
  plot(factors[, i] - factors.mean[i, ], type = "p", pch = 20, col = rgb(0, 0, 0.5, 0.3), 
       main = paste("Factor", i), ylim = range(factors[, i] - factors.mean[i,],
                                               factors.lwr[i,]-factors.mean[i,], 
                                               factors.upr[i,]-factors.mean[i,]))
  points(factors.lwr[i, ]-factors.mean[i, ], type = "l", lty = 2, col = 2)
  points(factors.upr[i, ]-factors.mean[i, ], type = "l", lty = 2, col = 2)
}

mean(factors > t(factors.upr))
mean(factors < t(factors.lwr))
plot(factors[,1], factors.mean[1,], pch = 20,col = rgb(0,0,0.5,0.3))
abline(coef = c(0, 1), col = "red")
plot(factors[,2], factors.mean[2,], pch = 20,col = rgb(0,0,0.5,0.3))
abline(coef = c(0, 1), col = "red")

summary(lm(factors[,1]~factors.mean[1,]))
summary(lm(factors[,2]~factors.mean[2,]))
