
#ifndef drmFFBSdiscW_H
#define drmFFBSdiscW_H

//[[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

// funcao FFBS para regressao com V constante e fator de desconto
using namespace Rcpp;
using namespace arma;

// funcao FFBS para regressao com V constante e fator de desconto
Rcpp::List drmFFBSdiscW(mat Y, mat X, mat dV, double discW, vec m0, mat ZC0);

#endif