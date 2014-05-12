
#ifndef FactorSim_H
#define FactorSim_H

//[[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

using namespace Rcpp;
using namespace arma;

//funcao para simular fatores
mat FactorSim(mat Y, mat Lambda, vec psi);

#endif