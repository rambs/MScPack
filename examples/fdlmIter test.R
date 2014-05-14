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
tmp = fdlmInitQuant(modelo, FALSE)
modeloNA = modelo
modeloNA$y[c(50, 100, 150, 200, 250), 1:3] = NA
modeloNA$y[c(80, 180, 280, 380, 480), 6:9] = NA
tmpf = fdlmGibbsIter(modeloNA)
str(tmpf)
ls(environment(tmpf))
str(tmpf())
str(tmpf(F))
str(get("ival", envir = environment(tmpf)))
get("ival", envir = environment(tmpf))[["Y_NA"]]
tmp1 = tmpf()
tmp1_ = get("ival", envir = environment(tmpf))[["theta"]]
tmp2 = tmpf()
tmp2_ = get("ival", envir = environment(tmpf))[["theta"]]
str(tmp1)
all.equal(tmp1_, tmp2_)
system.time(for (i in 1:1000){
  tmpf()
})
str(x)
print(all.equal(tmp1_, tmp2_))
9.05*100/60
system.time(fdlmGibbs(1,1000,1,modelo, get("ival", envir = environment(tmpf)), T, T))

fdlmGibbsIter <- function(model){
    Y = model[["y"]]
    xFix = model[["xFixReg"]]
    xDyn = model[["xDynReg"]]
    k = model[["nFactors"]]
    
    LL0 = model[["L0"]] #hiperparametros de Lambda
    HH0 = model[["H0"]] 
    
    bb0 = model[["b0"]] #hiperparametros de Beta
    BB0 = model[["B0"]]
    
    mm0 = model[["m0"]] #hiperparametro de Theta
    
    nn0 = model[["n0"]] #hiperparametros da variancia idiossincratica
    ss0 = model[["s0sq"]]
    
    delta = modelo[["discW"]] #fator de desconto
    
    initQuant = fdlmInitQuant(model, FALSE)
    for (xn in names(initQuant)) assign(xn, initQuant[[xn]])
    
    q = ncol(Y) # matriz de dados deve entrar tempo x variavel
    TT = nrow(Y)
    pFix = ncol(xFix)
    pDyn = ncol(xDyn)
    rDyn = pDyn*q;
    
    #valores iniciais da cadeia
    ival = list(theta = array(mm0, c(rDyn, TT+1)),
                    MuDyn = array(NA, c(TT, q)), #media dinamica
                    Lambda = array(LL0, c(q, k)),
                    Factors = array(0, c(TT, k)),
                    psi = array(ss0, c(q, 1)),
                    Beta = array(model[["b0"]], c(pFix, q)),
                    MuFix = array(xFix %*% model[["b0"]], c(TT, q)), #media fixa
                    ZW = array(NA, c(rDyn, rDyn)))
    function(runOnly = TRUE){
      # simulacao dos estados
      E = Y - ival[["MuFix"]] - ival[["Factors"]] %*% t(ival[["Lambda"]])
      thList = MScPack:::.drmFFBSdiscW(E, xDyn, diag(ival[["psi"]][,1]), delta, mm0, ZC0)
      ival[["theta"]] <<- thList[["th"]]
      ival[["MuDyn"]] <<- thList[["Mu"]]
      ival[["ZW"]] <<- thList[["ZW"]]
      
      # simulacao da reg estatica
      E = Y - ival[["MuDyn"]] - ival[["Factors"]] %*% t(ival[["Lambda"]]);
      ival[["Beta"]] <<- MScPack:::.parmsFixRegSim(E, xFix, B1, ZB1, ival[["psi"]], B0Invb0);
      ival[["MuFix"]] <<- xFix %*% ival[["Beta"]];
      
      # simulacao dos fatores
      E = Y - ival[["MuFix"]] - ival[["MuDyn"]]
      ival[["Factors"]] <<- MScPack:::.FactorSim(E, ival[["Lambda"]], ival[["psi"]])
      
      # simulacao da matriz de cargas
      ival[["Lambda"]] <<- MScPack:::.LambdaSimMV(E, ival[["Factors"]], ival[["psi"]], 
                                                L0H0Inv, H0Inv)
      
      # simulacao das variancias idiossincraticas
      E = E - ival[["Factors"]] %*% t(ival[["Lambda"]])
      ival[["psi"]] <<- MScPack:::.psiSim(E, ival[["Beta"]], bb0, B0Inv, 
                                        ival[["Lambda"]], LL0, H0Inv, nn0, ss0)
      
      if(!runOnly)
        return(ival)
    }
}

    
    
 
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
