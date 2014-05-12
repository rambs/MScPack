#ifndef psiSim_H
#define psiSim_H

//[[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

using namespace Rcpp;
using namespace arma;

// simulacao das variancias idiossincraticas
vec psiSim(mat Y, mat Beta, mat b0, mat B0Inv, mat Lambda, mat L0, mat H0Inv, 
            double n0, vec s0);

#endif