
library(MScPack)
library(MCMCpack)
library(truncnorm)
sourceCpp("src/SampleSv.cpp")
sourceCpp("src/SampleSvPhi.cpp")
sourceCpp("src/SampleSvMean.cpp")
sourceCpp("src/SampleSvVariance.cpp")

# simulacao MVE ####
mu = c(-11, -10, -9)
phi = c(0.98, 0.6, 0.9)
W = matrix(c(5, 1, -1, 1, 2, 0, -1, 0, 2), 3, 3)
U = W - diag(phi) %*% W %*% diag(phi)

TT = 1000
set.seed(23907)
eta = matrix(rnorm(TT*length(mu)), length(mu), TT)
At = t(chol(U))
for (i in 1:TT) {
  if (i <= 1){
    eta[, i] = mu + t(chol(W)) %*% eta[, i]
    next
  }
  eta[, i] = mu + phi * (eta[, i-1] - mu) + At %*% eta[, i]
}

y = t(apply(exp(eta/2), 2, rnorm, n = 3, mean = 0))
y1 = y[1:100, ]
y2 = y[1:500, ]

par(mfrow = c(4, 3), mar = c(2.1, 2.1, 2.1, 0.1))
invisible(apply(y1, 2, plot, type = "l", ylab = "", xlab = "Tempo", main = "T = 100"))
invisible(apply(y2, 2, plot, type = "l", ylab = "", xlab = "Tempo", main = "T = 500"))
invisible(apply(y, 2, plot, type = "l", ylab = "", xlab = "Tempo", main = "T = 1000"))
invisible(apply(exp(eta/2), 1, plot, type = "l", ylab = "", xlab = "Tempo", 
                main = expression(exp(eta[t]/2))))


N <- 1e2
eta.sim <- array(NA, c(dim(eta), N))
system.time(
for (i in 1:N){
  if(i <= 1){
    eta.sim[,, i] <- SampleSv(y, eta, mu, as.matrix(phi), U)  
    next
  }
  eta.sim[,, i] <- SampleSv(y, eta.sim[,, i-1], mu, as.matrix(phi), U)  
  if(! i %% (N/10)){
    message(paste("Iteracao", i, "concluida"))
  }
}
)

sourceCpp("src/SampleSvItMn.cpp")
N <- 1e2
eta.sim <- array(NA, c(dim(eta), N))
system.time(
  for (i in 1:N){
    if(i <= 1){
      eta.sim[,, i] <- SampleSvItMn(y, eta, mu, as.matrix(phi), U)  
      next
    }
    eta.sim[,, i] <- SampleSvItMn(y, eta.sim[,, i-1], mu, as.matrix(phi), U)  
    if(! i %% (N/10)){
      message(paste("Iteracao", i, "concluida"))
    }
  }
)

eta.mean <- apply(eta.sim, c(1, 2), mean)
eta.lwr <- apply(eta.sim, c(1, 2), quantile, probs = 0.05)
eta.upr <- apply(eta.sim, c(1, 2), quantile, probs = 0.95)

par(mfrow = c(3, 1))
for (i in 1:3){
  plot(eta.mean[i, ], type = "l", ylim = range(eta[i, ], eta.lwr[i, ], 
                                               eta.upr[i, ]))
  points(eta.lwr[i, ], type = "l", col = 3)
  points(eta.upr[i, ], type = "l", col = 3)
  points(eta[i, ], type = "l", col = 2)
}

mean(eta < eta.lwr)
mean(eta > eta.upr)


# Gibbs for eta, phi and mu -----------------------------------------------
sourceCpp("src/SampleSvVariance.cpp")

sourceCpp("src/tmp_vePhiSim.cpp")

sourceCpp("src/SampleSv.cpp")
sourceCpp("src/tmp_veFFBSv03.cpp")
SampleSv(y, eta, mu, phi, U)/veFFBSv03(y, eta, mu, phi, U)

tmp1 <- array(NA, c(3, 3, 100))
tmp2 <- array(NA, c(3, 3, 100))
for (i in 1:100){
  tmp1[,, i] <- SampleSv(y, eta, mu, phi, U)
  tmp2[,, i] <- veFFBSv03(y, eta, mu, phi, U)
}
for (i in 1:3){
  for (j in 1:3){
    hist(tmp1[i, j, ], xlim = range(tmp1[i, j, ], tmp2[i, j, ]))
    par(new = T)
    hist(tmp2[i, j, ])
  }
}

N <- 1e4
eta.sim <- array(NA, c(dim(eta), N))
phi.sim <- array(NA, c(N, 3))
mu.sim <- array(NA, c(N, 3))
U.sim <- array(NA, c(3, 3, N))
system.time(
  for (i in 1:N){
    if(i <= 1){
      eta.sim[,, i] <- SampleSv(y, eta, mu, phi, U)
      mu.sim[i, ] <- SampleSvMean(eta.sim[,, i], phi, U, 0, 1e6)
      phi.sim[i, ] <- SampleSvPhi(phi, eta.sim[,, i], mu.sim[i, ], U)      
      U.sim[,, i] <- SampleSvVariance(U, eta.sim[,, i], phi.sim[i, ], 
                                      mu.sim[i, ], 0.1, diag(0.1, 3))
      next
    }
    eta.sim[,, i] <- SampleSv(y, eta.sim[,, i-1], mu.sim[i-1, ], 
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

# resultados para comparar com relatorio na pasta 
# /users/u3rb/desktop/modelos dinamicos fatoriais copia/Rpres/
plot(as.mcmc(mu.sim))
plot(as.mcmc(phi.sim))
U.mcmc <- t(apply(U.sim, 3, function(x) x[lower.tri(x, TRUE)]))
colnames(U.mcmc) <- c("U1_1", "U1_2", "U1_3", "U2_2", "U2_3", "U3_3")
plot(as.mcmc(U.mcmc))


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

par(mfrow = c(3, 1))
for (i in 1:3){
  plot(phi.sim[, i], type = "l")
  abline(h = phi[i], col = 2)
}

for (i in 1:3){
  hist(phi.sim[, i], col = "darkblue", breaks = 25)
  abline(v = phi[i], col = 2)
}

par(mfrow = c(3, 1))
for (i in 1:3){
  plot(mu.sim[, i], type = "l")
  abline(h = mu[i], col = 2)
}

for (i in 1:3){
  hist(mu.sim[, i], col = "darkblue", breaks = 50, xlim = range(mu.sim))
  abline(v = mu[i], col = 2)
}

par(mfrow = c(3, 3))
for (i in 1:3){
  for (j in 1:3){
    if (i >= j){
      hist(U.sim[i, j, ], col = "darkblue", breaks = 100, xlim = range(U.sim),
           border = "darkblue")
      abline(v = U[i, j], col = 2)  
      next
    } 
    plot(U.sim[i, j, ], type = "l")
    abline(h = U[i, j], col = 2)  
  }  
}

Ucov.sim <- apply(U.sim, 3, cov2cor)
dim(Ucov.sim) <- c(3, 3, N)
par(mfrow = c(3, 3))
for (i in 1:3){
  for (j in 1:3){
    if (i > j){
      hist(Ucov.sim[i, j, ], col = "darkblue", breaks = 100, xlim = range(Ucov.sim),
           border = "darkblue")
      abline(v = cov2cor(U)[i, j], col = 2)  
      next
    } 
    plot(U.sim[i, j, ], type = "l")
    abline(h = U[i, j], col = 2)  
  }  
}

par(mfrow = c(3, 3))
for (i in 1:3){
  for (j in 1:3){
    plot(U.sim[i, j, ], type = "l")
    abline(h = U[i, j], col = 2)  
  }  
}