

# dados artificiais VAR(3) - DFM(2)####
TT <- 1000
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
factors <- array(rnorm((TT+h)*k), c(TT+h, k))
phi = rep(0, r)
for (i in 1:(TT+h)){
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
f.star <- as.vector(t(factors[1:h,]))
for (i in 1:TT){
  f.star <- c(factors[i+h, ], f.star[1:(k*s)])
  y[i, ] <- LambdaBar %*% f.star + y[i, ]
}
par(mfrow = c(3, 3))
apply(y, 2, plot, type = "l")

# calling function ####
require(Rcpp)
src <- paste(readLines("src/SampleDfmLoads.cpp"), collapse = "\n")
sourceCpp(code = src)

N = 1e3
LambdaBar.sim <- array(NA, c(q, k*(s+1), N))
for (i in 1:N){
  LambdaBar.sim[,, i] <- SampleDfmLoads(y, factors, psi, s, 1e1)  
}

par(mfrow = dim(LambdaBar), mar = c(2.1, 2.1, 1.0, 0.5))
for (i in 1:q){
  for (j in 1:(k*(s+1))){
    hist(LambdaBar.sim[i, j, ], breaks = 25, col = "darkblue", 
         main = "")#paste("LambdaBar[", i, ",", j, "]", sep = ""))
    abline(v = LambdaBar[i, j], col = "red")
  }
}


# sampling factors and VAR parms in a Gibbs sampler ####
src <- paste(readLines("src/SampleDynFactors.cpp"), collapse = "\n")
sourceCpp(code = src)
src <- paste(readLines("src/SampleVarParms.cpp"), collapse = "\n")
sourceCpp(code = src)

N <- 1e3
burn <- 1e4
LambdaBar.sim <- array(0, c(q, k*(s+1), N))
PhiBar.sim <- array(0, c(k, k*h, N))
factors.sim <- array(0, c(k, TT+h, N))
for (j in 1:burn){
  LambdaBar.sim[,, 1] <- SampleDfmLoads(y, t(factors.sim[,, 1]), psi, s, 1e1)  
  factors.sim[,, 1] <- SampleDynFactors(y, LambdaBar.sim[,, 1], PhiBar.sim[,, 1], psi)
  PhiBar.sim[,, 1] <- SampleVarParms(t(factors.sim[,, 1]), h = h, c0 = 1e3)
}
system.time(
  for (j in 1:(N-1)){
    LambdaBar.sim[,, j+1] <- SampleDfmLoads(y, t(factors.sim[,, j]), psi, s, 1e1)  
    factors.sim[,, j+1] <- SampleDynFactors(y, LambdaBar.sim[,, j+1], 
                                            PhiBar.sim[,, j], psi)
    PhiBar.sim[,, j+1] <- SampleVarParms(t(factors.sim[,, j+1]), h = h, c0 = 1e3)
  }
)

par(mfrow = dim(LambdaBar), mar = c(2.1, 2.1, 1.0, 0.5))
for (i in 1:q){
  for (j in 1:(k*(s+1))){
    hist(LambdaBar.sim[i, j, ], breaks = 25, col = "darkblue", 
         main = "")#paste("LambdaBar[", i, ",", j, "]", sep = ""))
    abline(v = LambdaBar[i, j], col = "red")
  }
}

for (i in 1:q){
  for (j in 1:(k*(s+1))){
    plot(LambdaBar.sim[i, j, ], type = "l", col = "darkblue", 
         main = "")#paste("LambdaBar[", i, ",", j, "]", sep = ""))
    abline(h = LambdaBar[i, j], col = "red")
  }
}

library(MScPack)
qr.Lambda0 <- qr(t(LambdaBar[, 1:2]))
Lambda0.plt <- t(qr.R(qr.Lambda0))
Q.plt <- t(qr.Q(qr.Lambda0))
Lambda1.plt <- LambdaBar[, 3:4] %*% Q.plt
Lambda2.plt <- LambdaBar[, 5:6] %*% Q.plt
LambdaBar.plt <- cbind(Lambda0.plt, Lambda1.plt, Lambda2.plt)

WOP <- wop.fdlm(LambdaBar.sim[, 1:2,])

LambdaBarWOP <- array(NA, dim(LambdaBar.sim))
for (i in 1:N){
  LambdaBarWOP[, 1:2, i] <- WOP$Lambda[,, i]
  LambdaBarWOP[, 3:4, i] <- LambdaBar.sim[, 3:4, i] %*% WOP$D[,, i]
  LambdaBarWOP[, 5:6, i] <- LambdaBar.sim[, 5:6, i] %*% WOP$D[,, i]
}

PLT <- plt.fdlm(LambdaBarWOP[,1:2,])
str(PLT)
LambdaBarPLT <- array(NA, dim(LambdaBar.sim))
for (i in 1:N){
  LambdaBarPLT[, 1:2, i] <- PLT$Lambda[,, i]
  LambdaBarPLT[, 3:4, i] <- LambdaBarWOP[, 3:4, i] %*% PLT$Q[,, i]
  LambdaBarPLT[, 5:6, i] <- LambdaBarWOP[, 5:6, i] %*% PLT$Q[,, i]
}

for (i in 1:q){
  for (j in 1:(k*(s+1))){
    plot(LambdaBarPLT[i, j, ], type = "l", col = "darkblue", 
         main = "")#paste("LambdaBar[", i, ",", j, "]", sep = ""))
    abline(h = LambdaBar.plt[i, j], col = "red")
  }
}


par(mfrow = c(2, 6), mar = c(2.1, 2.1, 1.0, 0.5))
for (i in 1:k){
  for (j in 1:(h*k)){
    hist(PhiBar.sim[i, j, ], breaks = 25, col = "darkblue", 
         main = paste("PhiBar[", i, ",", j, "]", sep = ""),
         xlim = range(PhiBar.sim))
    abline(v = PhiBar[i, j], col = "red")
  }
}
