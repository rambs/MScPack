# MScPack
# Description: DFM's Gibbs sampler and ex-post rotation
# Author: Rafael Barcellos
# Last updated: 26th June, 2014
# R 3.1.0

# loading required packages -----------------------------------------------

library(MScPack)

# sourcing DFM(2, 3) ------------------------------------------------------

source("tests/RandGenDfmSVarH.R")
500
1

# sampling loads, factors, and VAR parms ----------------------------------

N <- 1e3
thin <- 2
burn <- 5e3

LambdaBar.sim <- array(NA, c(q, k*(s+1), N))
PhiBar.sim <- array(NA, c(k, k*h, N))
factors.sim <- array(NA, c(k, TT+h, N))

# working objects with initial values
Ls <- array(0, c(q, k*(s+1)))
diag(Ls) <- 1
Ps <- array(0, c(k, k*h))
fs <- array(0, c(k, TT+h))

# progress bar
pb <- txtProgressBar(style = 3)
it <- 0
n.it <- burn + thin*N
start.time <- Sys.time()

# gibbs sampler
for (j in 1:N){
  if(j == 1){
    max.count <- burn + thin
  } else {
    max.count <- thin
=======

# sampling loads, factors, and VAR parms ----------------------------------

N <- 1e4
thin <- 2
burn <- 1e4

LambdaBar.sim <- array(NA, c(q, k*(s+1), N))
PhiBar.sim <- array(NA, c(k, k*h, N))
factors.sim <- array(NA, c(k, TT+h, N))

# working objects with initial values
Ls <- array(0, c(q, k*(s+1)))
diag(Ls) <- 1
Ps <- array(0, c(k, k*h))
fs <- array(0, c(k, TT+h))

# progress bar
pb <- txtProgressBar(style = 3)
it <- 0
n.it <- burn + thin*N
start.time <- Sys.time()

# gibbs sampler
for (j in 1:N){
  if(j == 1){
    max.count <- burn + thin
  } else {
    max.count <- thin
  }
  for (i in 1:max.count){
    fs <- SampleDynFactors(y, Ls, Ps, psi)
    Ls <- SampleDfmLoads(y, t(fs), psi, s, 1e1)
    Ps <- SampleVarParms(t(fs), h = h, c0 = 1e3)
    it <- it + 1
    setTxtProgressBar(pb, it/n.it)
  }
  factors.sim[,, j] <- fs
  LambdaBar.sim[,, j] <- Ls
  PhiBar.sim[,, j] <- Ps
}
end.time <- Sys.time()
duration <- end.time - start.time
duration

# LambdaBar posterior draws -----------------------------------------------

par(mfrow = dim(LambdaBar), mar = c(0.1, 0.1, 0.1, 0.1))
for (i in 1:q){
  for (j in 1:(k*(s+1))){
    plot(LambdaBar.sim[i, j, ], type = "l", col = "darkblue", 
         main = "", ylim = range(LambdaBar.sim[i, j, ], 
                                 LambdaBar[i, j]))
    abline(h = LambdaBar[i, j], col = "red")
>>>>>>> 78330ae168a348465f3583da84526799b54a60b7
  }
  for (i in 1:max.count){
    fs <- SampleDynFactors(y, Ls, Ps, psi)
    Ls <- SampleDfmLoads(y, t(fs), psi, s, 1e1)
    Ps <- SampleVarParms(t(fs), h = h, c0 = 1e3)
    it <- it + 1
    setTxtProgressBar(pb, it/n.it)
  }
  factors.sim[,, j] <- fs
  LambdaBar.sim[,, j] <- Ls
  PhiBar.sim[,, j] <- Ps
}
end.time <- Sys.time()
duration <- end.time - start.time
duration

# LambdaBar posterior draws -----------------------------------------------

par(mfrow = dim(LambdaBar), mar = c(0.1, 0.1, 0.1, 0.1))
for (i in 1:q){
  for (j in 1:(k*(s+1))){
    hist(LambdaBar.sim[i, j, ], breaks = 25, col = "darkblue", 
         main = "", xlim = range(LambdaBar.sim, LambdaBar[i, j]))
    abline(v = LambdaBar[i, j], col = "red")
  }
}

par(mfrow = dim(LambdaBar), mar = c(0.1, 0.1, 0.1, 0.1))
for (i in 1:q){
  for (j in 1:(k*(s+1))){
    hist(LambdaBar.sim[i, j, ], breaks = 25, col = "darkblue", 
         main = "", xlim = range(LambdaBar.sim, LambdaBar[i, j]))
    abline(v = LambdaBar[i, j], col = "red")
  }
}

# PhiBar posterior draws --------------------------------------------------

par(mfrow = c(2, 6), mar = c(2.1, 2.1, 0.1, 0.1))
for (i in 1:k){
  for (j in 1:(h*k)){
    plot(PhiBar.sim[i, j, ], col = "darkblue", type = "l",
         main = "", #paste("PhiBar[", i, ",", j, "]", sep = ""),
         ylim = range(PhiBar.sim))
    abline(h = PhiBar[i, j], col = "red")
  }
}

par(mfrow = c(2, 6), mar = c(2.1, 2.1, 1.0, 0.5))
for (i in 1:k){
  for (j in 1:(h*k)){
    hist(PhiBar.sim[i, j, ], breaks = 25, col = "darkblue", 
         main = "", #paste("PhiBar[", i, ",", j, "]", sep = ""),
         xlim = range(PhiBar.sim))
    abline(v = PhiBar[i, j], col = "red")
  }
}

# Factors posterior summaries ---------------------------------------------

factors.mean <- apply(factors.sim, c(1, 2), mean)
factors.lwr <- apply(factors.sim, c(1, 2), quantile, probs = 0.05)
factors.upr <- apply(factors.sim, c(1, 2), quantile, probs = 0.95)

TpH <- ncol(factors.mean) # factors' time interval = T + h

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
       lwd = 1.5, col= "red")
  axis(4, las = 1)
}

# The factors' posterior paths are similar to those of the true
#   factors. 
# We need to find the arguments that minimize the ex-post loss 
#   function.

# estimated levels of the series ------------------------------------------

level.sim <- array(NA, c(dim(y), N))
for (j in 1:N){
  x <- NULL
  for (i in 1:(s+1)){
    x <- cbind(x, t(factors.sim[,(h+2-i):(TpH-i+1), j]))
  }
  level.sim[,, j] <- x %*% t(LambdaBar.sim[,, j])  
}

level.mean <- apply(level.sim, c(1, 2), mean)
level.lwr <- apply(level.sim, c(1, 2), quantile, probs = 0.05)
level.upr <- apply(level.sim, c(1, 2), quantile, probs = 0.95)

x <- NULL
for (i in 1:(s+1)){
  x <- cbind(x, factors[(h+2-i):(TpH-i+1), ])
}
true.level <- x %*% t(LambdaBar)

par(mfrow = c(3, 3))
for (i in 1:q){
  plot(true.level[, i], type = "l")
  points(level.mean[, i], type = "l", col = 2)
  xcoord <- c(1:TT, rev(1:TT))
  ycoord <- c(level.lwr[, i], rev(level.upr[, i]))
  #polygon(xcoord, ycoord, border = "darkblue", col = rgb(0, 0, 0.5, 0.001))
}

level.mcmc <- level.sim
dim(level.mcmc) <- c(TT*q, N)
level.mcmc <- t(level.mcmc)
dim(level.mcmc)

dim(LambdaBar.sim)
tmp <- LambdaBar.sim[, 1:k,]
tmp <- apply(tmp, 3, tcrossprod)
tmp <- t(tmp)

LLt <- tcrossprod(LambdaBar[, 1:k])
par(mfrow = c(q, q), mar = rep(0.1, 4))
for (i in 1:(q^2)){
  plot(tmp[, i], type = "l")
  abline(h = as.vector(LLt)[i], col = 2)
}

library(MCMCpack)
gd <- geweke.diag(tmp)
rd <- raftery.diag(tmp)

par(mfrow = c(1, 1), mar = c(5.1, 4.1, 4.1, 2.1))
plot(abs(gd$z)) # no case grater than 1 in absolute values
abline(h = 1)

# quadratic loss optimization ---------------------------------------------

# initial values
Lstar <- LambdaBar.sim[,, N]
dim(Lstar) <- c(q, k, s+1)
Lstar <- do.call("rbind", lapply(1:(s+1), function(i) Lstar[,,i]))
Pstar <- PhiBar.sim[,, N]

LambdaD <- array(NA, c(dim(Lstar), N))
Phi.rtd <- array(NA, dim(PhiBar.sim))
D <- array(NA, c(ncol(Lstar), ncol(Lstar), N))

it <-0
tol <- 1e-10
max.iter <- 10
eps <- 1

pb <- txtProgressBar(style = 3)

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
    out.opt <- tmpExPostLossOptim(Ls, Lstar, Ps, Pstar)#, n.values = 0)
    LambdaD[,, j]  <- out.opt$Lambda.opt
    D[,, j] <- out.opt$D.opt
    Phi.rtd[,, j] <- out.opt$Phi.opt
    setTxtProgressBar(pb, j/N)
  }
  Lstar <- apply(LambdaD, c(1, 2), mean)
  Pstar <- apply(Phi.rtd, c(1, 2), mean)
  eps <- sum((Lstar-Lstar0)^2) + sum((Pstar-Pstar0)^2)
  it <- it + 1
  message(paste("Iteracao", it, "concluida."))
}

# testing optimisation function -------------------------------------------

tmp <- RunDfmExPostOptAlg(LambdaBar.sim, PhiBar.sim)
all.equal(tmp$Lambda, LambdaD)
str(tmp)

attach(tmp)

# analyzing optimisation results ------------------------------------------

LambdaS <- LambdaBar
dim(LambdaS) <- c(q, k, s+1)
LambdaS <- do.call("rbind", lapply(1:(s+1), function(i) LambdaS[,,i]))

true.opt <- ExPostLossOptim(LambdaS, tmp$Lstar, PhiBar, tmp$Pstar)
D.wop <- true.opt$D.opt
LambdaS.opt <- true.opt$Lambda.opt
PhiBar.opt <- true.opt$Phi.opt

<<<<<<< HEAD
par(mfrow = dim(LambdaBar), mar = rep(0.1, 4))
for (i in 1:nrow(LambdaS)){
  for (j in 1:ncol(LambdaS)){
    plot(LambdaD[i, j, ], type = "l", col = rgb(0, 0, 0.5))
    abline(h = LambdaS.opt[i, j], col = "red")
  }
}

par(mfrow = dim(LambdaBar), mar = rep(0.1, 4))
for (i in 1:nrow(LambdaS)){
  for (j in 1:ncol(LambdaS)){
     hist(LambdaD[i, j, ], col = rgb(0, 0, 0.5, 0.3), breaks = 30, 
          xlim = range(LambdaD), main = "", border = rgb(0,0,0.5,0.3))
     par(new = T)
    hist(tmp$Lambda[i, j, ], col = rgb(0, 0, 0.5, 0), breaks = 30, 
         xlim = range(tmp$Lambda), main = "", border = rgb(0,0,0,0.3))
    abline(v = LambdaS.opt[i, j], col = "red")
  }
}

=======
par(mfrow = dim(LambdaBar), mar = rep(0.1, 4))
for (i in 1:nrow(LambdaS)){
  for (j in 1:ncol(LambdaS)){
    plot(LambdaD[i, j, ], type = "l", col = rgb(0, 0, 0.5))
    abline(h = LambdaS.opt[i, j], col = "red")
  }
}

par(mfrow = dim(LambdaBar), mar = rep(0.1, 4))
for (i in 1:nrow(LambdaS)){
  for (j in 1:ncol(LambdaS)){
     hist(LambdaD[i, j, ], col = rgb(0, 0, 0.5, 0.3), breaks = 30, 
          xlim = range(LambdaD), main = "", border = rgb(0,0,0.5,0.3))
     par(new = T)
    hist(tmp$Lambda[i, j, ], col = rgb(0, 0, 0.5, 0), breaks = 30, 
         xlim = range(tmp$Lambda), main = "", border = rgb(0,0,0,0.3))
    abline(v = LambdaS.opt[i, j], col = "red")
  }
}

>>>>>>> 78330ae168a348465f3583da84526799b54a60b7
plot(LambdaD[1, 1,], LambdaD[2, 1, ])
plot(LambdaD[3, 2,], LambdaD[2, 2, ])
plot(LambdaD[3, 1,], LambdaD[4, 1, ])

par(mfrow = c(2, 6), mar = c(2.1, 2.1, 1.0, 0.5))
for (i in 1:k){
  for (j in 1:(h*k)){
    plot(Phi.rtd[i, j, ], col = "darkblue", type = "l",
         main = "", 
         ylim = range(Phi.rtd))
    abline(h = PhiBar.opt[i, j], col = "red")
  }
}

par(mfrow = c(2, 6), mar = c(2.1, 2.1, 1.0, 0.5))
for (i in 1:k){
  for (j in 1:(h*k)){
    hist(Phi.rtd[i, j, ], col = "darkblue", breaks = 30,
         main = "")#, xlim = range(Phi.rtd))
    abline(v = PhiBar.opt[i, j], col = "red")
  }
}
<<<<<<< HEAD

plot(Phi.rtd[2, 1, ], Phi.rtd[2, 3, ])
points(PhiBar.opt[2, 1], PhiBar.opt[2, 3], pch = 16, col = 2)

plot(Phi.rtd[2, 1, ], Phi.rtd[2, 6, ])
points(PhiBar.opt[2, 1], PhiBar.opt[2, 6], pch = 16, col = 2)
=======
>>>>>>> 78330ae168a348465f3583da84526799b54a60b7

# applying PLT to original loads ------------------------------------------

LambdaS <- LambdaBar
dim(LambdaS) <- c(q, k, s+1)
LambdaS <- do.call("rbind", lapply(1:(s+1), function(i) LambdaS[,,i]))

qr.LambdaS <- qr(t(LambdaS))
LambdaS.plt <- t(qr.R(qr.LambdaS))
reflexion <- diag(sign(diag(LambdaS.plt)))
LambdaS.plt <- LambdaS.plt %*% reflexion
D.plt <- qr.Q(qr.LambdaS) %*% reflexion # %*% reflexion

LambdaPLT <- plt.fdlm(LambdaD)$Lambda
<<<<<<< HEAD
D.sim.plt <- plt.fdlm(LambdaD)$D
=======
D.sim.plt <- plt.fdlm(LambdaD)$Q
>>>>>>> 78330ae168a348465f3583da84526799b54a60b7

par(mfrow = c(q, k*(s+1)), mar = rep(0.1, 4))
for (i in 1:(q*(s+1))){
  for (j in 1:k){
    plot(LambdaPLT[i, j, ], type = "l", col = "darkblue", 
         main = "")#paste("LambdaBar[", i, ",", j, "]", sep = ""))
    abline(h = LambdaS.plt[i, j], col = "red")
  }
}

par(mfrow = c(q, k*(s+1)), mar = c(2.1, 2.1, 0.1, 0.1))
lmat <- matrix(1:(k*(s+1)*q), ncol = k, byrow = T)
lmat <- t(lmat)
dim(lmat) <- c(k, q, (s+1))
lmat <- t(Reduce("rbind", lapply(1:(s+1), function(i) lmat[,, i])))
layout(lmat)
for (i in 1:(q*(s+1))){
  for (j in 1:k){
    hist(LambdaPLT[i, j, ], col = "darkblue", xlim = range(LambdaPLT),
         main = "", breaks = 100)
    abline(v = LambdaS.plt[i, j], col = "red")
  }
}

# analyzing bimodal results -----------------------------------------------

det.D <- apply(D, 3, det)
plot(det.D)
# the determinants are all equal to one. It must be studied if some of them
# should be negative to avoid the bimodal behaviour of the posterior.

Phi.plt <- RotatePhi(PhiBar, D.plt)
Phi.sim.plt <- array(NA, c(k, k*h, N))
for (i in 1:N){
  Phi.sim.plt[,, i] <- RotatePhi(Phi.rtd[,,i], D.sim.plt[,, i])
}

par(mfrow = c(2, 6), mar = c(2.1, 2.1, 1.0, 0.5))
for (i in 1:k){
  for (j in 1:(h*k)){
    plot(Phi.sim.plt[i, j, ], col = "darkblue", type = "p",
         main = "", #paste("PhiBar[", i, ",", j, "]", sep = ""),
         ylim = range(Phi.sim.plt))
    abline(h = Phi.plt[i, j], col = "red")
  }
}

par(mfrow = c(2, 6), mar = c(2.1, 2.1, 1.0, 0.5))
for (i in 1:k){
  for (j in 1:(h*k)){
<<<<<<< HEAD
    hist(Phi.sim.plt[i, j, ], col = "darkblue", breaks = 50, 
         border = rgb(0,0,0.5),
         main = "")#, #paste("PhiBar[", i, ",", j, "]", sep = ""),
         #xlim = range(Phi.sim.plt))
=======
    hist(Phi.sim.plt[i, j, ], col = "darkblue",
         main = "", #paste("PhiBar[", i, ",", j, "]", sep = ""),
         xlim = range(Phi.sim.plt))
>>>>>>> 78330ae168a348465f3583da84526799b54a60b7
    abline(v = Phi.plt[i, j], col = "red")
  }
}

x <- Phi.sim.plt
dim(x) <- c(k*k*h, N)
x <- t(x)
par(mfrow = c(11, 6),mar = rep(1.0, 4))
for (i in 1:11){
  for (j in (i+1):12){
    plot(x[, i], x[, j], pch = 16, 
         col = rgb(0, 0, 0.5, 0.2))
  }
}

# factors after opt -------------------------------------------------------

factors.sim.opt <- factors.sim
for (j in 1:N){
  factors.sim.opt[,, j] <- t(D[,, j]) %*% factors.sim[,, j]
}

dim(factors.sim.opt)
par(mfrow = c(1, 1))
plot(factors.sim.opt[1, 1, ], type = "l")
plot(factors.sim.opt[1, 503, ], type = "l")
hist(factors.sim.opt[1, 503, ], breaks = 100)
hist(factors.sim.opt[1, 502, ], breaks = 100, main = "")
hist(factors.sim.opt[2, 502, ], breaks = 100, main = "")
plot(factors.sim.opt[1, 502, ], factors.sim.opt[2, 502, ])
plot(factors.sim.opt[1, 503, ], factors.sim.opt[2, 503, ])
plot(factors.sim.opt[1, 50, ], factors.sim.opt[2, 50, ])
plot(factors.sim.opt[1, 150, ], factors.sim.opt[2, 150, ])
plot(factors.sim.opt[1, 250, ], factors.sim.opt[2, 250, ])

factors.opt.mean <- apply(factors.sim.opt, c(1, 2), mean)
factors.opt.lwr <- apply(factors.sim.opt, c(1, 2), quantile, probs = 0.05)
factors.opt.upr <- apply(factors.sim.opt, c(1, 2), quantile, probs = 0.95)

Ls <- LambdaBar
dim(Ls) <- c(q, k, s+1)
Ls <- do.call("rbind", lapply(1:(s+1), function(i) Ls[,, i]))
Ps <- PhiBar
opt.parms <- ExPostLossOptim(Ls, Lstar, Ps, Pstar, n.values = 10)
D.opt <- opt.parms$D.opt
factors.opt <- factors %*% D.opt

par(mfrow = c(2, 1), mar = c(3.1, 3.1, 0.1, 3.1))
for (i in 1:2){
  plot(factors.opt.mean[i, ], type = "n", ylab = "", xlab = "", las = 1,
       ylim = range(factors.opt.lwr[i, ], factors.opt.upr[i, ]), col = "gold")
  xcoord <- c(1:TpH, TpH:1)
  ycoord <- c(factors.opt.lwr[i, ], rev(factors.opt.upr[i, ]))
  polygon(xcoord, ycoord, border = "darkblue", col = rgb(0, 0, 0.5, 0.5))
  points(factors.opt.mean[i, ], type = "l", col = "gold")
  #par(new = T)
  points(factors.opt[, i], type = "l", lty = 2, 
         lwd = 1.5, col= "red")
  axis(4, las = 1)
}
mean(factors.opt < t(factors.opt.lwr))
mean(factors.opt > t(factors.opt.upr))

# factors after PLT -------------------------------------------------------

factors.sim.plt <- factors.sim
for (j in 1:N){
  factors.sim.plt[,, j] <- t(D.sim.plt[,, j]) %*% factors.sim[,, j]
}

dim(factors.sim.plt)
par(mfrow = c(1, 1))
plot(factors.sim.plt[1, 1, ], type = "l")
plot(factors.sim.plt[1, 503, ], type = "l")
hist(factors.sim.plt[1, 503, ], breaks = 100)
hist(factors.sim.plt[1, 502, ], breaks = 100, main = "")
hist(factors.sim.plt[2, 502, ], breaks = 100, main = "")
plot(factors.sim.plt[1, 502, ], factors.sim.plt[2, 502, ])
plot(factors.sim.plt[1, 503, ], factors.sim.plt[2, 503, ])
plot(factors.sim.plt[1, 50, ], factors.sim.plt[2, 50, ])
plot(factors.sim.plt[1, 150, ], factors.sim.plt[2, 150, ])
plot(factors.sim.plt[1, 250, ], factors.sim.plt[2, 250, ])
plot(factors.sim.plt[1, 350, ], factors.sim.plt[2, 350, ])
plot(factors.sim.plt[1, 450, ], factors.sim.plt[2, 450, ])

factors.plt.mean <- apply(factors.sim.plt, c(1, 2), mean)
factors.plt.lwr <- apply(factors.sim.plt, c(1, 2), quantile, probs = 0.05)
factors.plt.upr <- apply(factors.sim.plt, c(1, 2), quantile, probs = 0.95)

factors.plt <- factors %*% D.plt

par(mfrow = c(2, 1), mar = c(3.1, 3.1, 0.1, 3.1))
for (i in 1:2){
  plot(factors.plt.mean[i, ], type = "n", ylab = "", xlab = "", las = 1,
       ylim = range(factors.plt.lwr[i, ], factors.plt.upr[i, ]), col = "gold")
  xcoord <- c(1:TpH, TpH:1)
  ycoord <- c(factors.plt.lwr[i, ], rev(factors.plt.upr[i, ]))
  polygon(xcoord, ycoord, border = "darkblue", col = rgb(0, 0, 0.5, 0.5))
  points(factors.plt.mean[i, ], type = "l", col = "gold")
  #par(new = T)
  points(factors.plt[, i], type = "l", lty = 2, 
         lwd = 1.5, col= "red")
  axis(4, las = 1)
}


