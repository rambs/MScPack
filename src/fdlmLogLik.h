
#ifndef fdlmLogLik_H
#define fdlmLogLik_H

//[[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace arma;

// funcao que calcula a verossimilhanca marginal integrando analiticamente theta e Factors
arma::vec fdlmLogLik(arma::mat Y, arma::mat X, arma::cube Lambda, arma::mat psi, 
double discW);

#endif
