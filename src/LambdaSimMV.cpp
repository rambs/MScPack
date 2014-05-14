#include <RcppArmadillo.h>
using namespace arma;
//funcao para simular Lambda Matriz-Variada
//[[Rcpp::export(".LambdaSimMV")]]
arma::mat LambdaSimMV(arma::mat Y, arma::mat Factors, arma::vec psi, 
arma::mat L0H0Inv, arma::mat H0Inv){
  int q = Y.n_cols;
  int k = Factors.n_cols;
  
  arma::mat H1Inv = (Factors.t()*Factors + H0Inv);
  arma::mat Eigvec;
  arma::vec eigval;
  eig_sym(eigval, Eigvec, H1Inv);
  arma::mat L1 = (Y.t()*Factors + L0H0Inv)*Eigvec*diagmat(1/eigval)*Eigvec.t();
  
  arma::mat ZHt = diagmat(1/sqrt(eigval))*Eigvec.t();
  arma::mat Lambda = L1 + diagmat(sqrt(psi))*randn(q, k)*ZHt;
  return Lambda;
}