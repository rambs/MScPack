# MScPack
# Descricao: Funcao fdlm iteracao
# Autor: Rafael Barcellos
# Data: 09/05/2014
# R 3.0.2

#' fdlm - iteracao do Gibbs
#' 
#' Funcao a ser utilizada para se ter o resultado de uma iteracao do algoritmo de Gibbs.
#' Com seu uso e' possivel aplicar o modelo fdlm a muitas series a alocar as amostras da
#' posteriori a um objeto 'ff'.
#' @param model objeto tipo 'list' que contem os dados de entrada, as variaveis exogenas e os hiperparametros das prioris.
#' @param initVal valores iniciais para a cadeia.
#' 

library(MScPack)
tmpfdlmIter <- function(model, init.val){
  local{
    initVal = init.val
    function(iterVal = NULL){
      if(!is.null(iterVal))
        initVal = iterVal
      fdlmGibbs(1, 0, 1, model, initVal, FALSE, TRUE)  
    }  
  }
  
}

tmpfdlmIterMemo <- local({
  initVal = NULL
  function(model, init.val){
    initVal <<- init.val
    function(iterVal = NULL){
      if(!is.null(iterVal))
        initVal = iterVal
      fdlmGibbs(1, 0, 1, model, initVal, FALSE, TRUE)  
    }  
  }
})

source("c:/users/u3rb/desktop/Programações em R/Modelo dinamico fatorial/R/mdfve funcoes pacote.R")

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
fdlmInitQuant(modelo, F)
str(mdfSim)
str(modelo)
init.val = list(Lambda = mdfSim$mod$Lambda, psi = as.matrix(mdfSim$mod$psi), Beta = mdfSim$mod$parmsFixReg)
init.val

system.time(tmpb <- fdlmGibbs(1, 10000, 1, modelo, init.val, FALSE, FALSE))
tmp = tmpfdlmIter(modelo)

system.time(for(i in 1:1000){
  tmpa<-tmp(init.val)
})
str(tmpa)
