
#ifndef drmFFBSdiscW_chol_H
#define drmFFBSdiscW_chol_H

//[[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

// funcao FFBS para regressao com V constante e fator de desconto
using namespace arma;

// funcao FFBS para regressao com V constante e fator de desconto utilizando Choleski
Rcpp::List drmFFBSdiscW_chol(arma::mat Y, arma::mat X, arma::mat dV, double discW, 
arma::vec m0, arma::mat C0_chol, arma::mat C0Inv);

#endif