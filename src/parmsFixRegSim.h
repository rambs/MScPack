#ifndef parmsFixRegSim_H
#define parmsFixRegSim_H

//[[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

using namespace Rcpp;
using namespace arma;

//funcao para simulacao dos parametros de regressao estaticos
mat parmsFixRegSim(mat Y, mat X, mat B1, mat ZB1, vec psi, mat B0Invb0);

#endif