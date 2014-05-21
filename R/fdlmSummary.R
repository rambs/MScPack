# Pacote MScPack
# Descricao: Funcoes resumo para fdlm
# Autor: Rafael Barcellos
# Data: 20/05/2014
# R 3.0.2

#' Medidas resumo parametros de estado
#' 
#' Calcula a media e os quantis de interesse para os parametros de estado simulados.
#' @param theta array com as matrizes de estados simuladas via MCMC;
#' @param sig.level nivel de credibilidade do intervalo
#' @return \code{array} cujas laminas sao media, limite inferior e limite superior.
theta.summary <- function(theta, sig.level = 0.95){
  th.mean = apply(theta, c(1, 2), mean)
  th.lwr = apply(theta, c(1, 2), quantile, probs = (1-sig.level)/2)
  th.upr = apply(theta, c(1, 2), quantile, probs = (1+sig.level)/2)
  th = array(c(th.mean, th.lwr, th.upr), c(dim(th.mean), 3))
  dimnames(th)[[3]] = c("mean", "lwr", "upr")
  th
}
