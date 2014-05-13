
#ifndef LambdaSimMV_H
#define LambdaSimMV_H

//[[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

//funcao para simular Lambda Matriz-Variada
arma::mat LambdaSimMV(arma::mat Y, arma::mat Factors, arma::vec psi, 
arma::mat L0H0Inv, arma::mat H0Inv);

#endif