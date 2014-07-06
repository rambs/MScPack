

# simulacao MVE ####
mu <- c(-11, -10, -9)
phi <- c(0.98, 0.6, 0.9)
W <- matrix(c(5, 1, -1, 1, 2, 0, -1, 0, 2), 3, 3)
U <- W - diag(phi) %*% W %*% diag(phi)

TT <- 1000
set.seed(23907)
eta <- matrix(rnorm(TT*length(mu)), length(mu), TT)
At <- t(chol(U))
for (i in 1:TT) {
  if (i <= 1){
    eta[, i] <- mu + t(chol(W)) %*% eta[, i]
    next
  }
  eta[, i] <- mu + phi * (eta[, i-1] - mu) + At %*% eta[, i]
}

factors <- t(apply(exp(eta/2), 2, rnorm, n = 3, mean = 0))

q <- 9
k <- 3
set.seed(6730)
Lambda <- array(runif(q*k, -1, 1), c(q, k))
Lambda[upper.tri(Lambda)] <- 0
diag(Lambda) <- 1

psi <- c(0.02, 0.19, 0.36, 0.02, 0.02, 
         0.19, 0.19, 0.36, 0.36)/1e3 # idiosyncratic variances
set.seed(8480)
y <- t(array(rnorm(q*TT, 0, sqrt(psi)), c(q, TT)))
y <- factors %*% t(Lambda) + y
par(mfrow = c(4, 3), mar = c(2.1, 2.1, 2.1, 0.1))
invisible(apply(y, 2, plot, type = "l", ylab = "", xlab = "Tempo", main = "T = 1000"))
invisible(apply(exp(eta/2), 1, plot, type = "l", ylab = "", xlab = "Tempo", 
                main = expression(exp(eta[t]/2))))


# calling functions -------------------------------------------------------

library(MScPack)
library(MCMCpack)
library(truncnorm)
sourceCpp("src/SampleSv.cpp")
sourceCpp("src/SampleSvPhi.cpp")
sourceCpp("src/SampleSvMean.cpp")
sourceCpp("src/SampleSvVariance.cpp")
sourceCpp("src/SampleFsvFactors.cpp")
sourceCpp("src/SampleFsvLoads.cpp")

# Gibbs algorithm ---------------------------------------------------------

N <- 1e4
factors.sim <- array(NA, c(k, TT, N))
Lambda.sim <- array(NA, c(q, k, N))
eta.sim <- array(NA, c(dim(eta), N))
phi.sim <- array(NA, c(N, 3))
mu.sim <- array(NA, c(N, 3))
U.sim <- array(NA, c(3, 3, N))
system.time(
  for (i in 1:N){
    if(i <= 1){
      factors.sim[,, i] <- SampleFsvFactors(y, Lambda, eta, psi)
      Lambda.sim[,, i] <- SampleFsvLoads(y, t(factors.sim[,, i]), psi, 0, 1e3)
      eta.sim[,, i] <- SampleSv(t(factors.sim[,,i]), eta, mu, phi, U)
      mu.sim[i, ] <- SampleSvMean(eta.sim[,, i], phi, U, 0, 1e6)
      phi.sim[i, ] <- SampleSvPhi(phi, eta.sim[,, i], mu.sim[i, ], U)      
      U.sim[,, i] <- SampleSvVariance(U, eta.sim[,, i], phi.sim[i, ], 
                                      mu.sim[i, ], 0.1, diag(0.1, 3))
      next
    }
    factors.sim[,, i] <- SampleFsvFactors(y, Lambda.sim[,, i-1], 
                                          eta.sim[,, i-1], psi)
    Lambda.sim[,, i] <- SampleFsvLoads(y, t(factors.sim[,, i]), psi, 0, 1e3)
    eta.sim[,, i] <- SampleSv(t(factors.sim[,,i]), eta.sim[,, i-1], mu.sim[i-1, ], 
                              phi.sim[i-1, ], U.sim[,, i-1]) 
    mu.sim[i, ] <- SampleSvMean(eta.sim[,, i], phi.sim[i-1, ], U.sim[,, i-1], 
                                0, 1e6)
    phi.sim[i, ] <- SampleSvPhi(phi.sim[i-1, ], eta.sim[,, i], 
                                mu.sim[i, ], U.sim[,, i-1])
    U.sim[,, i] <- SampleSvVariance(U.sim[,, i-1], eta.sim[,, i], phi.sim[i, ], 
                                    mu.sim[i, ], 0.1, diag(0.1, 3))
    if(! i %% (N/10)){
      message(paste("Iteracao", i, "concluida"))
    }
  }
)

factors.mean <- apply(factors.sim, c(1, 2), mean)
factors.lwr <- apply(factors.sim, c(1, 2), quantile, probs = 0.05)
factors.upr <- apply(factors.sim, c(1, 2), quantile, probs = 0.95)

par(mfcol = c(2, 3), mar = c(2.1, 2.1, 2.1, 0.5))
for (i in 1:3){
  plot(factors.mean[i, ], type = "l", 
       ylim = range(factors.lwr[i, ], factors.upr[i, ], factors[, i]))
  points(factors[, i], type = "l", lty = 2, col = 2)
  points(factors.lwr[i, ], type = "l", col = 3)
  points(factors.upr[i, ], type = "l", col = 3)
  plot(factors.mean[i, ], factors[, i], pch = 20, col = rgb(0, 0, 0.5, 0.3))
  abline(coef = c(0, 1), col = "red")
  abline(lm(factors[, i]~factors.mean[i, ]))
}

mean(t(factors) < factors.lwr)
mean(t(factors) > factors.upr)

eta.mean <- apply(eta.sim, c(1, 2), mean)
eta.lwr <- apply(eta.sim, c(1, 2), quantile, probs = 0.05)
eta.upr <- apply(eta.sim, c(1, 2), quantile, probs = 0.95)

par(mfrow = c(3, 1), mar = c(2.1, 2.1, 2.1, 0.5))
for (i in 1:3){
  plot(eta.mean[i, ], type = "l", ylim = range(eta[i, ], eta.lwr[i, ], 
                                               eta.upr[i, ]))
  points(eta.lwr[i, ], type = "l", col = 3)
  points(eta.upr[i, ], type = "l", col = 3)
  points(eta[i, ], type = "l", col = 2)
}

mean(eta < eta.lwr)
mean(eta > eta.upr)

# analyzing Lambda output -------------------------------------------------
par(mfcol = c(q, k), mar = c(2.1, 2.1, 0.5, 0.5))
for (j in 1:k){
  for (i in 1:q){
    if (j >= i){
      plot(c(0, 1), type = "n", axes = F, ylab = "", xlab = "", main = "")
    } else {
      lsim <- Lambda.sim[i, j, ]
      plot(lsim, type = "l", ylim = range(lsim, Lambda[i, j]))
      abline(h = Lambda[i, j], col = 2)  
    }
  }
}

par(mfcol = c(q, k), mar = c(2.1, 2.1, 0.5, 0.5))
for (j in 1:k){
  for (i in 1:q){
    if (j >= i){
      plot(c(0, 1), type = "n", axes = F, ylab = "", xlab = "", main = "")
    } else {
      lsim <- Lambda.sim[i, j, ]
      hist(lsim, col = rgb(0, 0, 0.5), xlim = range(lsim, Lambda[i, j]),
           breaks = 50)
      abline(v = Lambda[i, j], col = 2)  
    }
  }
}


# analyzing mcmc objects --------------------------------------------------

lmcmc <- t(apply(Lambda.sim, 3, function(x) x[lower.tri(x)]))
plot(as.mcmc(lmcmc))
plot(as.mcmc(mu.sim))
plot(as.mcmc(phi.sim))
U.mcmc <- t(apply(U.sim, 3, function(x) x[lower.tri(x, TRUE)]))
colnames(U.mcmc) <- c("U1_1", "U1_2", "U1_3", "U2_2", "U2_3", "U3_3")
plot(as.mcmc(U.mcmc))

size <- 0
size <- size + object.size(factors.sim)
size <- size + object.size(eta.sim)
size <- size + object.size(mu.sim)
size <- size + object.size(phi.sim)
size <- size + object.size(U.sim)
tmp <- 2.5*3*size
print(tmp, units= "Gb")
