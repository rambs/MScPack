
#ifndef LambdaSimMV_H
#define LambdaSimMV_H

//[[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

using namespace Rcpp;
using namespace arma;

//funcao para simular Lambda Matriz-Variada
mat LambdaSimMV(mat Y, mat Factors, vec psi, mat L0H0Inv, mat H0Inv);

#endif