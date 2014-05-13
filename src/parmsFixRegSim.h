#ifndef parmsFixRegSim_H
#define parmsFixRegSim_H

//[[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

//funcao para simulacao dos parametros de regressao estaticos
arma::mat parmsFixRegSim(arma::mat Y, arma::mat X, arma::mat B1, arma::mat ZB1, 
arma::vec psi, arma::mat B0Invb0);

#endif