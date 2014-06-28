# MScPack
# Description: WOP and PLT functions
# Author: Rafael Barcellos
# Last updated 26th June, 2014
# R 3.1.0

#' Matriz de cargas via WOP
#' 
#' Aplica metodo WOP para garantir unicidade da matriz de cargas.
#' @param Lambda array oriundo do MCMC;
#' @param max.iter maximo de iteracoes ate convergencia.
#' @return Lista com matriz ortogonal \code{D} e cargas tendo aplicado o WOP.
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

#' Matriz de cargas via PLT
#' 
#' Aplica restricao PLT na matriz de cargas do MCMC.
#' @param Lambda array com as simulacoes da matriz de cargas.
#' @return Lista com array de cargas sujeitos a PLT e matriz ortogonal da restricao.
plt.fdlm <- function(Lambda){
  N = dim(Lambda)[3]
  LambdaPLT = Lambda
  k = dim(Lambda)[2]
  D.plt = array(NA, c(k, k, N))
  for (r in 1:N){
    qr.Lambda = qr(t(Lambda[,,r]))
    LambdaPLT[,,r] = t(qr.R(qr.Lambda))
    reflexion = diag(sign(diag(LambdaPLT[,,r])))
    LambdaPLT[,,r] = LambdaPLT[,,r] %*% reflexion
    D.plt[,,r] = qr.Q(qr.Lambda) %*% reflexion
  }
  # Lambda' = Q %*% R
  # Lambda = R' %*% Q'
  # Lambda %*% Q = R'
  # Lambda %*% Q %*% reflexion = R' %*% reflexion = Lambda_PLT
  # :. D.plt = Q %*% reflexion
  return(list(LambdaPLT = LambdaPLT, D.plt = D.plt))
}