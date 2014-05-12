#include <RcppArmadillo.h>

using namespace Rcpp;
using namespace arma;
//funcao para simulacao dos parametros de regressao estaticos
mat parmsFixRegSim(mat Y, mat X, mat B1, mat ZB1, vec psi, mat B0Invb0)
{ 
  mat b1 = B1*(X.t()*Y + B0Invb0);
  mat beta = b1 + ZB1*randn(X.n_cols, Y.n_cols)*diagmat(sqrt(psi));
  return beta;
}