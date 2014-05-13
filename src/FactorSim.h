
#ifndef FactorSim_H
#define FactorSim_H

//[[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

using namespace Rcpp;

//funcao para simular fatores
arma::mat FactorSim(arma::mat Y, arma::mat Lambda, arma::vec psi);

#endif