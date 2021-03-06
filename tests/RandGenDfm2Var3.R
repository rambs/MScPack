##' MScPack
##' Description: Random generation of DFM with s = 2 and h = 3
##' Author: Rafael Barcellos
##' Last updated 21st June 2014
##' R 3.1.0

# defining parms ----------------------------------------------------------

TT <- 500 # time span
psi <- c(0.02, 0.19, 0.36, 0.02, 0.02, 
         0.19, 0.19, 0.36, 0.36, 
         0.02, 0.19, 0.36, 0.02, 0.02, 
         0.19, 0.19, 0.36, 0.36) # idiosyncratic variances
q <- length(psi) # number of variables
k <- 2 # number of factors
h <- 3 # VAR order
r <- h*k # number of state parms in FFBS

time.id <- (1-h):TT
Phi1 <- matrix(c(0, -0.1, 0, -0.7), 2, 2) + diag(rep(1, k)) # VAR parms
Phi2 <- matrix(c(0.02, -0.08, 0.4, -0.2), 2, 2)
Phi3 <- matrix(c(-0.06, 0.07, -0.6, 0.35), 2, 2)
PhiBar <- cbind(Phi1, Phi2, Phi3)

s <- scan(n = 1)
set.seed(4091) # generating loadings
LambdaBar  <- matrix(runif(k*q*(s+1), -0.1, 0.1), q, k*(s+1))

# generating factors ------------------------------------------------------

set.seed(7623) 
factors <- array(rnorm((TT+h)*k), c(TT+h, k))
phi = rep(0, r)
for (i in time.id+h){
  factors[i, ] <- PhiBar %*% phi + factors[i, ]
  phi <- c(factors[i, ], phi[1:(r-k)])
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
f.star <- as.vector(t(factors[h:(h-s),])) # auxiliar vector to handle lagged factors
for (i in 1:TT){
  f.star <- c(factors[i+h, ], f.star[-((k*s+1):(k*s+k))])
  y[i, ] <- LambdaBar %*% f.star + y[i, ]
}

par(mfrow = c(3, 3), mar = c(2.1, 2.1, 0.5, 0.5))
apply(y, 2, plot, x = time.id[-(1:h)], type = "l", 
      col = rgb(0, 0.5, 0, 0.5), bty = "l")
