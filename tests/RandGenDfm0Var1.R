##' MScPack
##' Description: Random generation of DFM with s = 0 and h = 1
##' Author: Rafael Barcellos
##' Last updated 05 June 2014
##' R 3.0.2


# defining parms ----------------------------------------------------------

TT <- 500 # time span
psi <- c(0.02, 0.19, 0.36, 0.02, 0.02, 0.19, 
         0.19, 0.36, 0.36) # idiosyncratic variances
q <- length(psi) # number of variables
k <- 2 # number of factors
h <- 1

time.id <- 0:TT
PhiBar <- matrix(c(1, -0.1, 0, 0.7), 2, 2) # VAR parms
set.seed(4091)
LambdaBar  <- matrix(runif(2*q, -0.1, 0.1), q, 2) # loadings

# generating factors ------------------------------------------------------

set.seed(7623)
factors <- array(rnorm((TT+h)*2), c(TT+h, 2))
for (i in time.id+1){
  if (i == 1){
    next
  }
  factors[i, ] <- PhiBar %*% factors[i-1, ] + factors[i, ]
}

par(mfrow = c(1, 1), mar = c(4.1, 4.2, 2.1, 1.0))
plot(factors, pch = 20, col = rgb(0, 0, 0.5, 0.5),
     xlab = "factor 1", ylab = "factor 2")

par(mfrow = c(1, 2), mar = c(2.1, 2.1, 2.1, 0.2))
apply(factors, 2, plot, x = time.id, type = "l", 
      col = rgb(0,0,0.5,0.5), bty = "l")

# generating data ---------------------------------------------------------

set.seed(3708)
y <- t(array(rnorm(q*TT, 0, sqrt(psi)), c(q, TT)))
y <- factors[-(1:h), ] %*% t(LambdaBar) + y

par(mfrow = c(3, 3), mar = c(2.1, 2.1, 0.5, 0.5))
apply(y, 2, plot, x = time.id[-(1:h)], type = "l", 
      col = rgb(0, 0.5, 0, 0.5), bty = "l")
