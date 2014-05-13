#include <RcppArmadillo.h>

using namespace arma;
//funcao para simulacao dos parametros de regressao estaticos
//[[Rcpp::export(".parmsFixRegSim")]]
arma::mat parmsFixRegSim(arma::mat Y, arma::mat X, arma::mat B1, arma::mat ZB1, 
arma::vec psi, arma::mat B0Invb0)
{ 
  mat b1 = B1*(X.t()*Y + B0Invb0);
  mat beta = b1 + ZB1*randn(X.n_cols, Y.n_cols)*diagmat(sqrt(psi));
  return beta;
}