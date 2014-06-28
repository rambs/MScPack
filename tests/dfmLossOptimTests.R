# MScPack
# Description: Testing the ex-post loss optimization
# Author: Rafael Barcellos
# Last updated 21st June, 2014
# R 3.0.1

# loading packages --------------------------------------------------------

library(MScPack)

# generating true parameters ----------------------------------------------

q <- 9 # number of variables
k <- 4 # number of factors
s <- 2 # number of lagged factors in observational equation
h <- 3 # number of factors' VAR process

set.seed(4091) # generating loadings
LambdaBar  <- matrix(runif(k*q*(s+1), -0.1, 0.1), q, k*(s+1))
Lambda <- as.vector(LambdaBar)
dim(Lambda) <- c(q, k, (s+1))
Lambda <- Reduce("rbind", lapply(1:(s+1), function(x) Lambda[,, x]))

set.seed(8092)
PhiBar <- array(rnorm(k*k*h, 0, 0.1), c(k, k*h))

# rotated parms -----------------------------------------------------------

set.seed(3478)
rotation.parms <- runif(k*(k-1)/2, -pi, pi)
D <- BuildOrthMat(rotation.parms)
LambdaStar <- Lambda %*% D + 0.01
PhiStar <- RotatePhi(PhiBar, D)

# optimizing quadratic loss -----------------------------------------------

rtd.loss <- BuildExPostRtdLoss(Lambda, LambdaStar, PhiBar, PhiStar)
n.values <- 200
set.seed(9046)
init.values <- array(runif(k*(k-1)/2*n.values), c(n.values, k*(k-1)/2))
init.values <- rbind(rep(0, k*(k-1)/2), init.values)
opt.loss <- apply(init.values, 1, function(x){
  optim(x, rtd.loss, reflexion = 1, method = "L-BFGS-B", lower = -pi+1e-14, 
        upper = pi-1e-14)
})

plot(sapply(opt.loss, function(x) x$val))
opt.posit <- which.min(sapply(opt.loss, function(x) x$val))
opt.loss[[opt.posit]]$val
opt.loss[[1]]$val
opt.loss[[opt.posit]]$par
opt.loss[[1]]$par
rotation.parms
BuildOrthMat(opt.loss[[opt.posit]]$par)/D
BuildOrthMat(opt.loss[[1]]$par)/D

BuildOrthMat(opt.loss[[opt.posit]]$par)/BuildOrthMat(opt.loss[[1]]$par)

# optimizing quadratic loss after WOP -------------------------------------

svd.S <- svd(t(Lambda) %*% LambdaStar)
D.wop <- svd.S$u %*% t(svd.S$v)
L1 <- sum((Lambda %*% D.wop- LambdaStar)^2)
sum((Lambda %*% D - LambdaStar)^2)
L2 <- sum((RotatePhi(PhiBar, D.wop) - PhiStar)^2)
L1+L2
rtd.loss <- BuildExPostRtdLoss(Lambda %*% D.wop, LambdaStar, 
                               RotatePhi(PhiBar, D.wop), PhiStar)
n.values <- 100
init.values <- array(runif(k*(k-1)/2*n.values), c(n.values, k*(k-1)/2))
init.values <- rbind(rep(0, k*(k-1)/2), init.values)
opt.loss <- apply(init.values, 1, function(x){
  optim(x, rtd.loss, reflexion = 1, method = "L-BFGS-B", lower = -pi+1e-14, 
        upper = pi-1e-14)
})

plot(sapply(opt.loss, function(x) x$val))
opt.posit <- which.min(sapply(opt.loss, function(x) x$val))
opt.loss[[opt.posit]]$val
opt.loss[[1]]$val
opt.loss[[opt.posit]]$par
opt.loss[[1]]$par
rotation.parms
D.wop %*% BuildOrthMat(opt.loss[[opt.posit]]$par)/D

BuildOrthMat(opt.loss[[opt.posit]]$par)/BuildOrthMat(opt.loss[[1]]$par)
