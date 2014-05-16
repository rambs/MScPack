source("c:/users/u3rb/Desktop/Programações em R/Modelo dinamico fatorial/R/mdfve funcoes pacote.R")

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
tmp = diag(tcrossprod(Lambda) + diag(psi))
diag(tcrossprod(Lambda))/tmp
mdfSim = mdfDiscW.sim(xFix, array(parmsFix, c(1, q)), xDyn, m0, C0, 0.95, Lambda, psi)
par(mfrow = c(3, 3), mar = c(2.1, 2.1, 0.1, 0.1))
apply(mdfSim$y, 2, plot, type = "l")
apply(exp(mdfSim$y), 2, plot, type = "l")


library(MScPack)
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
modeloNA = modelo
modeloNA$y[c(50, 100, 150, 200, 250), 1:3] = NA
modeloNA$y[c(80, 180, 280, 380, 480), 6:9] = NA
mdfGibbs = fdlmGibbsNA(10, 10, 1, modelo, init.val, T, F)
str(mdfGibbs)
mdfGibbsNA = fdlmGibbsNA(1000, 1000, 1, modeloNA, init.val, T, F)
str(mdfGibbsNA)

plot(modelo$y[is.na(modeloNA$y)], colMeans(mdfGibbsNA$values$Y_NA))
abline(coef = c(0,1))
class(mdfGibbsNA$Y_NA)
summary.fdlm <- function(object, sig.level = 0.95){
  Mu.mean = apply(object$gibbs$Mu, c(1, 2), mean)
  Mu.lwr = apply(object$gibbs$Mu, c(1, 2), quantile, probs = (1-sig.level)/2)
  Mu.upr = apply(object$gibbs$Mu, c(1, 2), quantile, probs = (1+sig.level)/2)
  Mu = array(c(Mu.mean, Mu.lwr, Mu.upr), c(dim(Mu.mean), 3))
  dimnames(Mu)[[3]] = c("mean", "lwr", "upr")
  Mu
}
tmp = summary.fdlm(mdfGibbs)
str(tmp)
tmp[,,"mean"]

Mu.mean = apply(mdfGibbs$gibbs$Mu, c(1, 2), mean)
for (i in 1:9){
  plot(mdfSim$y[,i], type = "p", ylim = range(Mu.mean[,i], mdfSim$y[,i]))
  points(Mu.mean[,i], type = "l", col = 2)
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


wop.fdlm <- function(object, max.iter = 100){
  q = ncol(object$model$y)
  k = object$model$nFactors
  N = object$gibbs$N
  Lambda.star = object$gibbs$Lambda[,,N]
  Lambda.0 = array(0, c(q, k))
  it = 0
  D = array(NA, c(k, k, N))
  LambdaWOP = array(NA, c(q, k, N))
  while(sum((Lambda.star-Lambda.0)^2)>1e-9){
    Lambda.0 = Lambda.star
    for (r in 1:N){
      Sr = t(object$gibbs$Lambda[,,r]) %*% Lambda.0
      Sr.svd = svd(Sr)
      D[,,r] = Sr.svd$u %*% t(Sr.svd$v)
      LambdaWOP[,,r] = object$gibbs$Lambda[,,r] %*% D[,,r]
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

plt.fdlm <- function(object){
  N = dim(object$LambdaWOP)[3]
  LambdaPLT = object$LambdaWOP
  Q.Lambda = object$D
  for (r in 1:N){
    qr.Lambda = qr(t(object$LambdaWOP[,,r]))
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
