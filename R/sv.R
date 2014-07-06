# MScPack
# Description: Get mixture of normals approximation parameters
# Author: Rafael Barcellos
# Last updated 3rd July, 2014
# R 3.1.0

#' @title Mixture of normals approximation
#' @description Get gaussian parameters which approximates the log-\eqn{\chi^2} 
#'   distribution. This approximation is made considering 10 gaussian kernels 
#'   such as those presented in Omori \emph{et al.} (2007).
#' @param y data matrix that follows a log-\eqn{\chi^2} distribution. The matrix
#'   must be entered with variables in columns.
#' @return A list containing mean and variance parameters of the normal
#'   distributions.
GetLogX2ApproxParms <- function(y){
  # Parms from Omori et al. (2007)
  p = c(0.00609, 0.04775, 0.13057, 0.20674, 0.22715,
        0.18842, 0.12047, 0.05591, 0.01575, 0.00115)
  
  m = c(1.92677, 1.34744, 0.73504, 0.02266, -0.85173,
        -1.97278, -3.46788, -5.55246, -8.68384, -14.65000)
  
  v2 = c(0.11265, 0.17788, 0.26768, 0.40611, 0.62699,
         0.98583, 1.57469, 2.54498, 4.16591, 7.33342)
  
  v = sqrt(v2)
  zSeq = 1:10
  
  z = sapply(1:ncol(y), function(k){
    sapply(y[,k], function(x){
      ppost = dnorm(x, m, v, TRUE) + log(p)
      ppost = exp(ppost-max(ppost)+700)
      sample(zSeq, 1, prob = ppost/sum(ppost))  
    })
  })
  M = array(m[z], dim(y))
  S = array(v2[z], dim(y))
  list(M = M, S = S)
}
