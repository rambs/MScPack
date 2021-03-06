#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;
// simulacao das variancias idiossincraticas
//[[Rcpp::export(".psiSim")]]
arma::vec psiSim(arma::mat Y, arma::mat Beta, arma::mat b0, arma::mat B0Inv, 
arma::mat Lambda, arma::mat L0, arma::mat H0Inv, double n0, arma::vec s0){
  Environment stats("package:stats"); //chamando a função rgamma do "stats"
  Function rgammaR = stats["rgamma"];
  
  int T = Y.n_rows;
  int q = Y.n_cols;
  int k = Lambda.n_cols;
  int p = Beta.n_rows;
  vec Spsi = sum(Y%Y).t() + trace((Beta-b0).t()*B0Inv*(Beta-b0)) + 
                trace((Lambda-L0)*H0Inv*(Lambda-L0).t());
  double n1 = n0 + T + k + p;
  vec n1s1 = (Spsi + n0*s0);
  vec psi = 1/as<arma::vec>(rgammaR(q, n1/2, n1s1/2));
  return psi;
}