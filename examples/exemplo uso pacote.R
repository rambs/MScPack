library(MScPack)

TT = 500
set.seed(293874)
xDyn = 3 + cumsum(rnorm(TT, 0, 0.05))
plot(xDyn, type = "l")

xFix = rep(1, TT)
parmsFix = -c(0.8, 0.92, 1.01, 1.1, 0.79, 0.98, 0.94, 1.07, 0.77)
q = length(parmsFix)
m0 = rep(1, q)
C0 = diag(0.01, q)
k = 3
Lambda.lim = 0.1
set.seed(1928)
Lambda = array(runif(q*k, 0, Lambda.lim), c(q, k))
Lambda[upper.tri(Lambda)] = 0
diag(Lambda) = c(0.99, 0.95, 0.9)*Lambda.lim
psi = c(0.02, 0.19, 0.36, 0.02, 0.02, 0.19, 0.19, 0.36, 0.36)*Lambda.lim/10
comunal = diag(tcrossprod(Lambda))
comunal/(comunal+psi)
mdfSim = mdfDiscW.sim(xFix, array(parmsFix, c(1, q)), xDyn, m0, C0, 0.95, Lambda, psi, 90165)
par(mfrow = c(3, 3), mar = c(2.1, 2.1, 0.1, 0.1))
apply(mdfSim$y, 2, plot, type = "l", bty = "l")
apply(exp(mdfSim$y), 2, plot, type = "l", bty = "l")



modelo = list(y = mdfSim$y, xFixReg = mdfSim$mod$xFixReg, xDynReg = mdfSim$mod$xDynReg,
              nFactors = ncol(mdfSim$mod$Lambda), 
              L0 = array(0, dim(mdfSim$mod$Lambda)), H0 = diag(1e2, ncol(mdfSim$mod$Lambda)),
              m0 = as.matrix(mdfSim$mod$m0), C0 = mdfSim$mod$C0,
              b0 = array(0, dim(mdfSim$mod$parmsFixReg)), B0 = diag(1e3, nrow(mdfSim$mod$parmsFixReg)),
              n0 = 1, s0sq = array(mdfSim$mod$psi, c(ncol(mdfSim$y), 1)),
              discW = mdfSim$mod$discW)
str(mdfSim)
str(modelo)
init.val = list(Lambda = mdfSim$mod$Lambda, psi = as.matrix(mdfSim$mod$psi), Beta = mdfSim$mod$parmsFixReg)
str(init.val)
MScPack:::.drmFFBSdiscW(modelo$y, modelo$xDynReg, diag(as.vector(modelo$s0sq)), 0.95,
                        modelo$m0, sqrt(modelo$C0))
MScPack:::.LambdaSimMV(modelo$y, array(rnorm(TT*k), c(TT, k)), modelo$s0sq, modelo$L0,
                       solve(modelo$H0))
MScPack:::.FactorSim(modelo$y, modelo$L0, modelo$s0sq)
MScPack:::.psiSim(modelo$y, modelo$b0, modelo$b0, solve(modelo$B0), modelo$L0, modelo$L0,
                  solve(modelo$H0), modelo$n0, modelo$s0sq)
                       
MScPack:::.parmsFixRegSim(modelo$y, modelo$xFixReg, modelo$B0, sqrt(modelo$B0), 
                          modelo$s0sq, solve(modelo$B0)%*%modelo$b0)

mdfGibbs = fdlmGibbs(1, 10000, 1, modelo, init.val)
str(mdfGibbs)
mdfGibbsOnly1 = fdlmGibbsOnly1(1, 10000, 1, modelo, init.val)
mdfGibbs$gibbs$dur
mdfGibbsOnly1$gibbs$dur
plot(mdfGibbs$values$theta[1,,1], type = "l", ylim = range(mdfGibbs$values$theta))
apply(mdfGibbs$values$theta[,,1], 1, points, type = "l")
apply(mdfGibbsOnly1$values$theta[,,1], 1, points, type = "l", col = 2)

par(mfrow = c(3,3), mar = c(2.1, 2.1, 0.1, 0.1))
for(i in 1:9){
  plot(mdfGibbs$values$Mu[,i,1], type = "l")
  points(mdfGibbsOnly1$values$Mu[,i,1], type = "l", col = 2)
}

# modelo com missings ####
modeloNA = modelo
modeloNA$y[c(50, 100, 150, 200, 250), 1:3] = NA
modeloNA$y[c(80, 180, 280, 380, 480), 6:9] = NA

mdfGibbsNA = fdlmGibbsNA(1000, 1000, 1, modeloNA, init.val, T, F)
mdfGibbs_chol = fdlmGibbs_chol(1000, 1000, 1, modeloNA, init.val, T, F)
mdfGibbsNA$gibbs$dur
mdfGibbs_chol$gibbs$dur
# funcao de verossimilhanca marginal
sourceCpp(paste("src/fdlmLogLik.cpp"))
sourceCpp(paste("src/fdlmLogLik_chol.cpp"))
sourceCpp(paste("src/fdlmLogLik_chol_mc.cpp"))
sourceCpp(paste("src/fdlmLogLik_mc.cpp"))
system.time(mdfLL <- fdlmLogLik(Y = mdfGibbsNA$mod$y, xFix = mdfGibbsNA$mod$xFixReg, 
           xDyn = mdfGibbsNA$mod$xDynReg, nFactors = mdfGibbsNA$mod$nFactors,
           hyperparms = list(m0 = mdfGibbsNA$mod$m0, C0 = mdfGibbsNA$mod$C0),
           Y_NA = mdfGibbsNA$values$Y_NA, whichNA = mdfGibbsNA$values$yNA,
           Beta = mdfGibbsNA$values$Beta, Lambda = mdfGibbsNA$values$Lambda, 
           psi = mdfGibbsNA$values$psi, discW = mdfGibbsNA$mod$discW))
system.time(mdfLL_chol <- 
              fdlmLogLik_chol(Y = mdfGibbs_chol$mod$y, xFix = mdfGibbs_chol$mod$xFixReg, 
                              xDyn = mdfGibbs_chol$mod$xDynReg, nFactors = mdfGibbs_chol$mod$nFactors,
                              hyperparms = list(m0 = mdfGibbs_chol$mod$m0, C0 = mdfGibbs_chol$mod$C0),
                              Y_NA = mdfGibbs_chol$values$Y_NA, whichNA = mdfGibbs_chol$values$yNA,
                              Beta = mdfGibbs_chol$values$Beta, Lambda = mdfGibbs_chol$values$Lambda, 
                              psi = mdfGibbs_chol$values$psi, discW = mdfGibbs_chol$mod$discW))


par(mfrow = c(1,2)); hist(mdfLL); hist(mdfLL_chol)
plot(density(mdfLL))
points(density(mdfLL_chol), type = "l", col = 2)
mean(mdfLL)
mean(mdfLL_chol)

# analise de convergencia das cadeias

# comunalidades
comunalMCMC = apply(mdfGibbsNA$values$Lambda, 3, function(x) diag(tcrossprod(x)))
comunalMCMC_chol = apply(mdfGibbs_chol$values$Lambda, 3, function(x) diag(tcrossprod(x)))
for(i in 1:q){
  plot(comunalMCMC[i,], type = "l", col = rgb(0,0,0,0.5), bty = "l",
       ylim = range(comunalMCMC[i,], comunalMCMC_chol[i,], comunal[i]))
  points(comunalMCMC_chol[i,], type = "l", col = rgb(0,0,0,0.5))
  abline(h = comunal[i], lwd = 2, lty = 2, col = "darkblue")
}

# variancias idiossincraticas
for(i in 1:q){
  plot(mdfGibbsNA$values$psi[,i], type = "l", col = rgb(0,0,0,0.5), bty = "l",
       ylim = range(mdfGibbsNA$values$psi[,i], mdfGibbs_chol$values$psi[,i], psi[i]))
  points(mdfGibbs_chol$values$psi[,i], type = "l", col = rgb(0,0,0,0.5))
  abline(h = psi[i], lwd = 2, lty = 1, col = "darkblue")
}

# parametros da regressao estatica
for(i in 1:q){
  plot(mdfGibbsNA$values$Beta[1,i,], type = "l", col = rgb(0,0,0,0.5), bty = "l",
       ylim = range(mdfGibbsNA$values$Beta[1,i,], mdfGibbs_chol$values$Beta[1,i,], 
                    parmsFix[i]))
  points(mdfGibbs_chol$values$Beta[1,i,], type = "l", col = rgb(0,0,0,0.5))
  abline(h = parmsFix[i], lwd = 2, lty = 1, col = "darkblue")
}

# relacao entre parametros de estado e os fixos
for(i in 1:q){
  plot(mdfGibbsNA$values$Beta[1,i,], mdfGibbsNA$values$theta[i,TT+1,], type = "p", 
       col = rgb(0,0,0,0.5), bty = "l", pch = 19,
       ylim = range(mdfGibbsNA$values$theta[i,TT+1,], mdfGibbs_chol$values$theta[i,TT+1,], 
                    mdfSim$th[i,TT+1]),
       xlim = range(mdfGibbsNA$values$Beta[1,i,], mdfGibbs_chol$values$Beta[1,i,], 
                    parmsFix[i]))
  points(parmsFix[i], mdfSim$th[i,TT+1], type = "p", col = "red", lwd = 2, pch = 10)
}

# media dinamica (mu_t)
par(mfrow = c(3,3))
for(i in 1:q){
  plot(mdfGibbsNA$values$Mu[TT,i,], type = "l", 
       col = rgb(0,0,0,0.5), bty = "l", 
       ylim = range(mdfGibbsNA$values$Mu[TT,i,], mdfGibbs_chol$values$Mu[TT,i,], 
                    mdfSim$y[TT,i]))
  points(mdfGibbs_chol$values$Mu[TT,i,], type = "l", 
       col = rgb(0,0,0,0.5))
  abline(h = mdfSim$y[TT,i], col = "red", lwd = 2)
}

# predicao dos valores faltantes
par(mfrow = c(1,1))
plot(modelo$y[is.na(modeloNA$y)], colMeans(mdfGibbsNA$values$Y_NA))
points(modelo$y[is.na(modeloNA$y)], colMeans(mdfGibbs_chol$values$Y_NA), col = 2, pch = 2)
abline(coef = c(0,1))

# matriz de cargas
par(mfrow = c(q,k))
for(i in 1:q){
  for(j in 1:k){
    plot(density(mdfGibbsNA$values$Lambda[i,j,]), col = rgb(0,0,0,0.5), bty = "l", 
         xlim = range(mdfGibbsNA$values$Lambda[i,j,], mdfGibbs_chol$values$Lambda[i,j,], 
                      Lambda[i,j]))
    points(density(mdfGibbs_chol$values$Lambda[i,j,]), type = "l", col = rgb(0,0,0,0.5))
    abline(v = Lambda[i,j], lwd = 2, col = 2)
  }
}

# aplicando WOP
LambdaWOP = wop.fdlm(mdfGibbsNA$values$Lambda)$Lambda
LambdaWOP_chol = wop.fdlm(mdfGibbs_chol$values$Lambda)$Lambda

par(mfrow = c(q,k))
for(i in 1:q){
  for(j in 1:k){
    plot(density(LambdaWOP[i,j,]), col = rgb(0,0,0,0.5), bty = "l", 
         xlim = range(LambdaWOP[i,j,], LambdaWOP_chol[i,j,], 
                      Lambda[i,j]))
    points(density(LambdaWOP_chol[i,j,]), type = "l", col = rgb(0,0,0,0.5))
    abline(v = Lambda[i,j], lwd = 2, col = 2)
  }
}

# aplicando PLT
LambdaPLT = plt.fdlm(LambdaWOP)$Lambda
LambdaPLT_chol = plt.fdlm(LambdaWOP_chol)$Lambda

par(mfrow = c(q,k))
for(i in 1:q){
  for(j in 1:k){
    plot(density(LambdaPLT[i,j,]), col = rgb(0,0,0,0.5), bty = "l", 
         xlim = range(LambdaPLT[i,j,], LambdaPLT_chol[i,j,], 
                      Lambda[i,j]))
    points(density(LambdaPLT_chol[i,j,]), type = "l", col = rgb(0,0,0,0.5))
    abline(v = Lambda[i,j], lwd = 2, col = 2)
  }
}
LambdaPLT[1,2,]
summary.fdlm <- function(object, sig.level = 0.95){
  Mu.mean = apply(object$values$Mu, c(1, 2), mean)
  Mu.lwr = apply(object$values$Mu, c(1, 2), quantile, probs = (1-sig.level)/2)
  Mu.upr = apply(object$values$Mu, c(1, 2), quantile, probs = (1+sig.level)/2)
  Mu = array(c(Mu.mean, Mu.lwr, Mu.upr), c(dim(Mu.mean), 3))
  dimnames(Mu)[[3]] = c("mean", "lwr", "upr")
  Mu
}
str(mdfGibbsNA)
Mu.summary = summary.fdlm(mdfGibbsNA)
Mu_chol.summary = summary.fdlm(mdfGibbs_chol)
par(mfrow = c(3,3), mar = c(2.1, 2.1, 0.1, 0.1))
for (i in 1:9){
  plot(Mu_chol.summary[,,"mean"][,i], type = "l", ylim = range(Mu.summary[,i,], Mu_chol.summary[,i,]))
  points(Mu.summary[,,"mean"][,i], type = "l", col = 2, lty = 2)
}

for (i in 1:9){
  plot(Mu_chol.summary[,,"lwr"][,i], type = "l", ylim = range(Mu.summary[,i,], Mu_chol.summary[,i,]))
  points(Mu.summary[,,"lwr"][,i], type = "l", col = 2, lty = 2)
}

for (i in 1:9){
  plot(Mu_chol.summary[,,"upr"][,i], type = "l", ylim = range(Mu.summary[,i,], Mu_chol.summary[,i,]))
  points(Mu.summary[,,"upr"][,i], type = "l", col = 2, lty = 2)
}

for (i in 1:9){
  plot(mdfSim$y[,i], type = "p", ylim = range(Mu.summary[,,"mean"][,i], mdfSim$y[,i]))
  points(Mu.summary[,,"mean"][,i], type = "l", col = 2)
}

res = mdfSim$y - Mu.mean
for (i in 1:9){
  plot(res[,i], type = "l")
}

for (i in 1:9){
  hist(res[,i])
}

th.mean = apply(mdfGibbs$gibbs$theta, c(1, 2), mean)
th.lwr = apply(mdfGibbs$gibbs$theta, c(1, 2), quantile, probs = 0.05)
th.upr = apply(mdfGibbs$gibbs$theta, c(1, 2), quantile, probs = 0.95)
for (i in 1:9){
  plot(mdfSim$th[i,], type = "p", ylim = range(th.upr[i,], th.lwr[i,],
                                               th.mean[i,], mdfSim$th[i,]))
  points(th.mean[i,], type = "l", col = 2)
  points(th.upr[i,], type = "l", col = 2)
  points(th.lwr[i,], type = "l", col = 2)
}

str(mdfGibbs$gibbs$Beta)
par(mfrow = c(3,3))
invisible(apply(mdfGibbs$gibbs$Beta, c(1,2), hist))
invisible(apply(mdfGibbs$gibbs$Beta, c(1,2), plot, type = "l"))
for(i in 1:q){
  plot(mdfGibbs$gibbs$Beta[1,i,], mdfGibbs$gibbs$theta[i, 50,])
  points(mdfSim$mod$parmsFixReg[,i], mdfSim$th[i, 50], pch = 19, col = 2)
}

apply(mdfGibbs$gibbs$psi, 2, hist)
psi


wop.fdlm <- function(Lambda, max.iter = 100){
  q = dim(Lambda)[1]
  k = dim(Lambda)[2]
  N = dim(Lambda)[3]
  Lambda.star = Lambda[,,N]
  Lambda.0 = array(0, c(q, k))
  it = 0
  D = array(NA, c(k, k, N))
  LambdaWOP = array(NA, c(q, k, N))
  while(sum((Lambda.star-Lambda.0)^2)>1e-9){
    Lambda.0 = Lambda.star
    for (r in 1:N){
      Sr = t(Lambda[,,r]) %*% Lambda.0
      Sr.svd = svd(Sr)
      D[,,r] = Sr.svd$u %*% t(Sr.svd$v)
      LambdaWOP[,,r] = Lambda[,,r] %*% D[,,r]
    }
    Lambda.star = apply(LambdaWOP, c(1,2), mean)
    it = it+1
    if(it>max.iter)
      break
  }
  return(list(LambdaWOP = LambdaWOP, D = D))  
}
res.wop = wop.fdlm(mdfGibbs)
par(mfcol = c(q, k))
apply(res.wop$LambdaWOP, c(1,2), plot, type = "l")
apply(mdfGibbs$gibbs$Lambda, c(1,2), plot, type = "l")

plt.fdlm <- function(Lambda){
  N = dim(Lambda)[3]
  LambdaPLT = Lambda
  k = dim(Lambda)[2]
  Q.Lambda = array(NA, c(k, k, N))
  for (r in 1:N){
    qr.Lambda = qr(t(Lambda[,,r]))
    LambdaPLT[,,r] = t(qr.R(qr.Lambda))
    reflexion = diag(sign(diag(LambdaPLT[,,r])))
    LambdaPLT[,,r] = LambdaPLT[,,r] %*% reflexion
    Q.Lambda[,,r] = reflexion %*% t(qr.Q(qr.Lambda))
  }
  return(list(LambdaPLT = LambdaPLT, Q = Q.Lambda))
}

res.plt = plt.fdlm(res.wop)

apply(res.plt$LambdaPLT, c(1,2), plot, type = "l")
apply(res.plt$LambdaPLT, c(1,2), hist, col = "darkblue", breaks = 50, xlim = range(res.plt$LambdaPLT))
Lambda
