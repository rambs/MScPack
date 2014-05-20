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
#' @param model objeto tipo 'list' que contem os dados de entrada, as variaveis exogenas 
#' e os hiperparametros das prioris.
#' @details Nessa funcao e' necessario especificar o objeto \code{model} com os 
#' seguintes nomes: 
#' \itemize{
#' \item \code{y}: matriz de dados (tempo x nVars);
#' \item \code{xFixReg}: matriz de variaveis explicativas da regressao estatica;
#' \item \code{xDynReg}: matriz de variaveis explicativas da regressao dinamica;
#' \item \code{nFactors}: numero de fatores latentes;
#' \item \code{L0}: media \emph{a priori} da matriz (nVars x nFactors) de cargas;
#' \item \code{H0}: matriz (nFactors x nFactors) de covariancias das 
#' linhas da matriz de cargas;
#' \item \code{m0}: media \emph{a priori} do vetor (nVars*nxDyn x 1) dos 
#' parametros de estado;
#' \item \code{C0}: variancia inicial dos parametros de estado;
#' \item \code{b0}: matriz de media dos parametros da regressao estatica;
#' \item \code{B0}: matriz (nxFix x nxFix) de covariancias dos parametros da regressao estatica;
#' \item \code{n0}: graus de liberdade \emph{a priori} das variancias idiossincraticas;
#' \item \code{s0sq}: vetor (nVars x 1) de medias harmonicas das variancias idiossincratica;
#' \item \code{discW}: fator de desconto dos parametros de estado.
#' }
#' @return function cujo o parametro e' um valor logico indicando se e para retornar 
#' os valores da cadeia ou somente avancar nas iteracoes.
#' 
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
  
  # tratamento de missings
  Y_NA = NULL
  yNA = which(is.na(Y))
  yNAsize = length(yNA)
  Ys = Y
  if(yNAsize>0){
    Y_NA = rep(0, yNAsize)
  }
  
  #valores iniciais da cadeia
  ival = list(theta = array(mm0, c(rDyn, TT+1)),
              MuDyn = array(NA, c(TT, q)), #media dinamica
              ZW = array(NA, c(rDyn, rDyn)),
              Lambda = array(LL0, c(q, k)),
              Factors = array(0, c(TT, k)),
              psi = array(ss0, c(q, 1)),
              Beta = array(model[["b0"]], c(pFix, q)),
              MuFix = array(xFix %*% model[["b0"]], c(TT, q)), #media fixa
              Mu = array(NA, c(TT, q)),
              Y_NA = Y_NA)
  function(runOnly = TRUE){
    if(yNAsize>0){
      Ys[yNA] <<- ival[["Y_NA"]]
    }
    # simulacao dos estados
    E = Ys - ival[["MuFix"]] - ival[["Factors"]] %*% t(ival[["Lambda"]])
    thList = MScPack:::.drmFFBSdiscW(E, xDyn, diag(ival[["psi"]][,1]), delta, mm0, ZC0)
    ival[["theta"]] <<- thList[["th"]]
    ival[["MuDyn"]] <<- thList[["Mu"]]
    ival[["ZW"]] <<- thList[["ZW"]]
    
    # simulacao da reg estatica
    E = Ys - ival[["MuDyn"]] - ival[["Factors"]] %*% t(ival[["Lambda"]]);
    ival[["Beta"]] <<- MScPack:::.parmsFixRegSim(E, xFix, B1, ZB1, ival[["psi"]], B0Invb0);
    ival[["MuFix"]] <<- xFix %*% ival[["Beta"]];
    
    # simulacao dos fatores
    E = Ys - ival[["MuFix"]] - ival[["MuDyn"]]
    ival[["Factors"]] <<- MScPack:::.FactorSim(E, ival[["Lambda"]], ival[["psi"]])
    
    # simulacao da matriz de cargas
    ival[["Lambda"]] <<- MScPack:::.LambdaSimMV(E, ival[["Factors"]], ival[["psi"]], 
                                                L0H0Inv, H0Inv)
    
    # simulacao das variancias idiossincraticas
    E = E - ival[["Factors"]] %*% t(ival[["Lambda"]])
    ival[["psi"]] <<- MScPack:::.psiSim(E, ival[["Beta"]], bb0, B0Inv, 
                                        ival[["Lambda"]], LL0, H0Inv, nn0, ss0)
    
    ival[["Mu"]] <<- ival[["MuFix"]] + ival[["MuDyn"]] + 
      ival[["Factors"]] %*% t(ival[["Lambda"]])
    if(yNAsize>0){
      psiMat = t(array(ival[["psi"]], c(TT, q)))
      ival[["Y_NA"]] <<- rnorm(yNAsize, ival[["Mu"]][yNA], sqrt(psiMat[yNA]))
    }
    if(!runOnly)
      return(ival)
  }
}

#' Simula F-DLM
#' 
#' Funcao que simula dados de um modelo dinamico fatorial
#' @param xFixReg matriz de delinameanto da regressao estatica;
#' @param parmsFixReg matriz de parametros da regressao estatica;
#' @param xDynReg matriz de delineamento da regressao dinamica;
#' @param m0,C0 media e variancia \emph{a priori} dos estados;
#' @param discW fator de desconto;
#' @param Lambda matriz de cargas fatoriais;
#' @param psi vetor de variancias idiossincraticas;
#' @param seed semente do gerador aleatorio.
#' @return Lista com os parametros de entrada, os parametros simulados e os dados artificiais.
mdfDiscW.sim <- function(xFixReg, parmsFixReg, xDynReg, m0, C0,
                         discW = 0.95, Lambda, psi, seed = as.numeric(Sys.time()))
{  
  if(!is.null(xFixReg))
    xFixReg = as.matrix(xFixReg) #matriz de exógenas T x pFix
  
  if(is.null(xFixReg) | length(xFixReg)==0)
    stop("Modelo deve conter algum regressor estatico.")
  
  if(!is.null(xDynReg)){
    xDynReg = as.matrix(xDynReg) #matriz de exógenas T x pDyn
    if(nrow(xFixReg) != nrow(xDynReg))
      stop("'xDynReg' deve conter mesmo número de linhas que 'xFixReg'.")
    #removendo exogenas que equivalem a intercepto
    is.xreg.intercept = apply(xDynReg, 2, function(x) all(x==x[1]))
    xDynReg = xDynReg[, !is.xreg.intercept, drop = F]
    if(any(is.xreg.intercept))
      warning("Em 'xDynReg': removidas as exogenas que equivalem a intercepto.")
  } 
  
  if(is.null(xDynReg) | length(xDynReg)==0)
    stop("Modelo deve conter algum regressor dinamico")
  
  
  #condicoes necessarias
  if(ncol(parmsFixReg) != nrow(Lambda))
    stop("'ncol(parmsFixReg)' deve ser igual a 'nrow(Lambda)'.")
  
  if(nrow(parmsFixReg) != ncol(xFixReg))
    stop("'nrow(parmsFixReg)' deve ser igual a 'ncol(xFixReg)'.")
  
  if(length(psi) != nrow(Lambda))
    stop("'length(psi)' deve ser igual a 'nrow(Lambda)'.")
  
  TT = nrow(as.matrix(xFixReg))
  q = nrow(Lambda)
  pFix = ncol(xFixReg)
  pDyn = ncol(xDynReg)
  k = ncol(Lambda)
  
  # geracao dos valores aleatorios
  set.seed(seed)
  h1 = TT*k
  h2 = pDyn*q*(TT+1)
  h3 = q*TT
  normsim = rnorm(h1+h2+h3)
  interval1 = 1:h1
  interval2 = (h1+1):(h1+h2)
  interval3 = (h1+h2+1):(h1+h2+h3)
  
  Factors = array(normsim[interval1], c(TT, k))
  th = array(normsim[interval2], c(pDyn*q, TT+1))
  E = diag(sqrt(psi)) %*% array(normsim[interval3], c(q, TT))
  y = t(E)
  #return(list(Factors = Factors, th = th, E = E))
  
  ZCsvd = svd(C0)
  ZC = ZCsvd$u %*% diag(sqrt(ZCsvd$d))
  
  th[,1] = m0 + ZC %*% th[,1]
  
  mm = m0
  for(i in 1:TT){
    aa = mm
    ZR = sqrt(1/discW)*ZC
    RR = tcrossprod(ZR)
    ZW = sqrt((1-discW)/discW)*ZC
    th[,i+1] = th[,i] + ZW %*% th[,i+1]
    muFix = t(parmsFixReg) %*% xFixReg[i,]
    FF = kronecker(diag(1, q), xDynReg[i,])
    muDyn = t(FF) %*% th[,i+1]
    y[i,] =  muFix + muDyn + Lambda %*% Factors[i,] + E[,i]
    ff = t(FF) %*% aa
    QQ = t(FF) %*% RR %*% FF + diag(psi)
    ee = y[i,] - muFix - Lambda %*% Factors[i,] - ff
    AA = RR %*% FF %*% solve(QQ)
    mm = aa + AA %*% ee
    L = t(ZR) %*% FF
    L.eig = eigen(L %*% diag((1/psi)) %*% t(L), symmetric = TRUE);
    #cat(paste("\n", all.equal(diag(1, nrow(L)), L.eig$u %*% t(L.eig$v)), "\n"))
    ZC = ZR %*% L.eig$vec %*% diag(1.0/sqrt(L.eig$val + 1.0));
  }
  mod = mget(names(formals(mdfDiscW.sim)))
  return(list(mod = mod, y = y, th = th, Factors = Factors))
}


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