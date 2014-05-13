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
