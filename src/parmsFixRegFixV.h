
#include <Rcpp.h>
//[[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
//funcao para simulacao dos parametros de regressao estaticos
mat parmsFixRegFixdV(mat Y, mat X, mat B1, mat ZB1, vec psi, mat B0Invb0)
{ 
  mat b1 = B1*(X.t()*Y + B0Invb0);
  mat beta = b1 + ZB1*randn(X.n_cols, Y.n_cols)*diagmat(sqrt(psi));
  return beta;
}