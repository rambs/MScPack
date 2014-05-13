#ifndef psiSim_H
#define psiSim_H

//[[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

// simulacao das variancias idiossincraticas
arma::vec psiSim(arma::mat Y, arma::mat Beta, arma::mat b0, arma::mat B0Inv, 
arma::mat Lambda, arma::mat L0, arma::mat H0Inv, double n0, arma::vec s0);

#endif