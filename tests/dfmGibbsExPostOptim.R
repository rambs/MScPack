#' MScPack
#' Description: DFM's Gibbs sampler and ex-post rotation
#' Author: Rafael Barcellos
#' Last updated: 19th June, 2014
#' R 3.1.0

# loading required packages -----------------------------------------------

library(MScPack)

# sourcing DFM(2, 3) ------------------------------------------------------

source("tests/RandGenDfm2Var3.R")

# sampling loads, factors, and VAR parms ----------------------------------

N <- 1e3
burn <- 5e3
LambdaBar.sim <- array(1, c(q, k*(s+1), N))
PhiBar.sim <- array(1, c(k, k*h, N))
factors.sim <- array(0, c(k, TT+h, N))
for (j in 1:burn){
  factors.sim[,, 1] <- SampleDynFactors(y, LambdaBar.sim[,, 1], PhiBar.sim[,, 1], psi)
  LambdaBar.sim[,, 1] <- SampleDfmLoads(y, t(factors.sim[,, 1]), psi, s, 1e1)  
  PhiBar.sim[,, 1] <- SampleVarParms(t(factors.sim[,, 1]), h = h, c0 = 1e3)
}
system.time(
  for (j in 1:(N-1)){
    factors.sim[,, j+1] <- SampleDynFactors(y, LambdaBar.sim[,, j], 
                                            PhiBar.sim[,, j], psi)
    LambdaBar.sim[,, j+1] <- SampleDfmLoads(y, t(factors.sim[,, j+1]), psi, s, 1e1)  
    PhiBar.sim[,, j+1] <- SampleVarParms(t(factors.sim[,, j+1]), h = h, c0 = 1e3)
  }
)

# LambdaBar posterior draws -----------------------------------------------

par(mfrow = dim(LambdaBar), mar = c(2.1, 2.1, 1.0, 0.5))
for (i in 1:q){
  for (j in 1:(k*(s+1))){
    hist(LambdaBar.sim[i, j, ], breaks = 25, col = "darkblue", 
         main = "", xlim = range(LambdaBar.sim[i, j, ], LambdaBar[i, j]))
    abline(v = LambdaBar[i, j], col = "red")
  }
}

for (i in 1:q){
  for (j in 1:(k*(s+1))){
    plot(LambdaBar.sim[i, j, ], type = "l", col = "darkblue", 
         main = "", ylim = range(LambdaBar.sim[i, j, ], 
                                 LambdaBar[i, j]))
    abline(h = LambdaBar[i, j], col = "red")
  }
}

# PhiBar posterior draws --------------------------------------------------

par(mfrow = c(2, 6), mar = c(2.1, 2.1, 1.0, 0.5))
for (i in 1:k){
  for (j in 1:(h*k)){
    hist(PhiBar.sim[i, j, ], breaks = 25, col = "darkblue", 
         main = paste("PhiBar[", i, ",", j, "]", sep = ""),
         xlim = range(PhiBar.sim))
    abline(v = PhiBar[i, j], col = "red")
  }
}

par(mfrow = c(2, 6), mar = c(2.1, 2.1, 1.0, 0.5))
for (i in 1:k){
  for (j in 1:(h*k)){
    plot(PhiBar.sim[i, j, ], col = "darkblue", type = "l",
         main = paste("PhiBar[", i, ",", j, "]", sep = ""),
         ylim = range(PhiBar.sim))
    abline(h = PhiBar[i, j], col = "red")
  }
}

# Factors posterior summaries ---------------------------------------------

factors.mean <- apply(factors.sim, c(1, 2), mean)
factors.lwr <- apply(factors.sim, c(1, 2), quantile, probs = 0.05)
factors.upr <- apply(factors.sim, c(1, 2), quantile, probs = 0.95)

par(mfrow = c(2, 1), mar = c(3.1, 3.1, 0.1, 3.1))
for (i in 1:2){
  plot(factors.mean[i, ], type = "n", ylab = "", xlab = "", las = 1,
       ylim = range(factors.lwr[i, ], factors.upr[i, ]), col = "gold")
  xcoord <- c(1:TpH, TpH:1)
  ycoord <- c(factors.lwr[i, ], rev(factors.upr[i, ]))
  polygon(xcoord, ycoord, border = "darkblue", col = rgb(0, 0, 0.5, 0.5))
  points(factors.mean[i, ], type = "l", col = "gold")
  par(new = T)
  plot(factors[, i], axes = FALSE, ylab = "", xlab = "", type = "l", lty = 2, 
       lwd = 2, col= "red")
  axis(4, las = 1)
}

# The behaviour of the paths of the factors' posterior are similar to those of the true
#  factors. We need to find the arguments that minimize the ex-post loss function.

library(MScPack)
# initial values
Lstar <- LambdaBar.sim[,, N]
dim(Lstar) <- c(q, k, s+1)
Lstar <- do.call("rbind", lapply(1:(s+1), function(i) Lstar[,,i]))
Pstar <- PhiBar.sim[,, N]

LambdaD <- array(NA, c(dim(Lstar), N))
PhiD <- array(NA, dim(PhiBar.sim))
D <- array(NA, c(ncol(Lstar), ncol(Lstar), N))

it <-0
tol <- 1e-10
max.iter <- 10
eps <- 1

while (eps > tol){
  if (it >= max.iter){
    break
  }
  Lstar0 <- Lstar
  Pstar0 <- Pstar
  for (j in 1:N){
    Ls <- LambdaBar.sim[,, j]
    dim(Ls) <- c(q, k, s+1)
    Ls <- do.call("rbind", lapply(1:(s+1), function(i) Ls[,, i]))
    Ps <- PhiBar.sim[,, j]
    out.opt <- ExPostLossOptim(Ls, Lstar, Ps, Pstar)
    LambdaD[,, j]  <- out.opt$Lambda.opt
    D[,, j] <- out.opt$D.opt
    PhiD[,, j] <- out.opt$Phi.opt
  }
  Lstar <- apply(LambdaD, c(1, 2), mean)
  Pstar <- apply(PhiD, c(1, 2), mean)
  eps <- sum((Lstar-Lstar0)^2)
  it <- it + 1
  cat(paste("\nIteracao", it, "concluida.\n"))
}
par(mfrow = dim(LambdaBar), mar = rep(0.1, 4))
apply(LambdaD, c(1, 2), plot, type = "l")
apply(PhiD, c(1, 2), plot, type = "l")


# applying PLT to original loads ------------------------------------------

LambdaS <- LambdaBar
dim(LambdaS) <- c(q, k, s+1)
LambdaS <- do.call("rbind", lapply(1:(s+1), function(i) LambdaS[,,i]))

qr.LambdaS <- qr(t(LambdaS))
LambdaS.plt <- t(qr.R(qr.LambdaS))
D.plt <- t(qr.Q(qr.LambdaS))

LambdaPLT <- plt.fdlm(LambdaD)$Lambda
Q.sim.plt <- plt.fdlm(LambdaD)$Q

dim(LambdaPLT)
par(mfrow = c(q*(s+1), k))
for (i in 1:(q*(s+1))){
  for (j in 1:k){
    plot(LambdaPLT[i, j, ], type = "l", col = "darkblue", 
         main = "")#paste("LambdaBar[", i, ",", j, "]", sep = ""))
    abline(h = LambdaS.plt[i, j], col = "red")
  }
}



par(mfrow = c(2, 6), mar = c(2.1, 2.1, 1.0, 0.5))
for (i in 1:k){
  for (j in 1:(h*k)){
    plot(PhiD[i, j, ], col = "darkblue", type = "l",
         main = paste("PhiBar[", i, ",", j, "]", sep = ""),
         ylim = range(PhiBar.sim))
    abline(h = PhiBar[i, j], col = "red")
  }
}

Phi.plt <- RotatePhi(PhiBar, D.plt)
Phi.sim.plt <- array(NA, c(k, k*h, N))
for (i in 1:N){
  Phi.sim.plt[,, i] <- RotatePhi(PhiD[,,i], D)
}
apply(PhiBar.sim, 3, RotatePhi, D = D.plt)
dim(Phi.sim.plt) <- c(k, k*h, N)


par(mfrow = c(2, 6), mar = c(2.1, 2.1, 1.0, 0.5))
for (i in 1:k){
  for (j in 1:(h*k)){
    plot(PhiD[i, j, ], col = "darkblue", type = "l",
         main = paste("PhiBar[", i, ",", j, "]", sep = ""),
         ylim = range(PhiBar.sim))
    abline(h = PhiBar[i, j], col = "red")
  }
}