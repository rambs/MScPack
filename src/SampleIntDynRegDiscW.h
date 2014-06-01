
#ifndef SampleIntDynRegDiscW_H
#define SampleIntDynRegDiscW_H

//[[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

// funcao FFBS para regressao com V constante e fator de desconto
using namespace arma;

// amosta intecepto e parametros de estado no FFBS
Rcpp::List SampleIntDynRegDiscW(arma::mat Y, arma::mat X, arma::mat dV, 
  double discW, arma::vec m0, arma::mat C0, arma::vec a0, arma::mat A0);

#endif