#include <RcppArmadillo.h>
using namespace arma;
//funcao para simular Lambda Matriz-Variada
//[[Rcpp::export(".LambdaSimMV")]]
arma::mat LambdaSimMV(arma::mat Y, arma::mat Factors, arma::vec psi, 
arma::mat L0H0Inv, arma::mat H0Inv){
  int q = Y.n_cols;
  int k = Factors.n_cols;
  //arma::mat H1 = inv(Factors.t()*Factors + H0Inv);
  arma::mat H1Inv = Factors.t()*Factors + H0Inv;
  arma::mat rootH1 = arma::trans(arma::inv(arma::trimatu(arma::chol(H1Inv))));
  //arma::mat H1Inv = symmatu(Factors.t()*Factors + H0Inv);
  //arma::mat Eigvec;
  //arma::vec eigval;
  //eig_sym(eigval, Eigvec, H1Inv);
  
  arma::mat L1 = (Y.t()*Factors + L0H0Inv)*rootH1.t()*rootH1;
  /*
  arma::mat U, V;
  arma::vec s;
  svd(U, s, V, H1, "standard");
  arma::mat ZHt = diagmat(sqrt(s))*V.t();*/
  arma::mat Lambda = L1 + diagmat(sqrt(psi))*randn(q, k)*rootH1;//ZHt;
  return Lambda;
}