
#ifndef drmOnly1FFBSdiscW_H
#define drmOnly1FFBSdiscW_H

//[[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

// funcao FFBS para regressao com V constante e fator de desconto
using namespace arma;

// funcao FFBS para regressao com V constante e fator de desconto
Rcpp::List drmOnly1FFBSdiscW(arma::mat Y, arma::mat X, arma::mat dV, double discW, 
arma::vec m0, arma::mat ZC0);

#endif