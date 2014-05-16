#include <RcppArmadillo.h>

using namespace arma;
//funcao para simular fatores
//[[Rcpp::export(".FactorSim")]]
arma::mat FactorSim(arma::mat Y, arma::mat Lambda, arma::vec psi){
  int k = Lambda.n_cols;
  int T = Y.n_rows;
  arma::mat G1 = inv(Lambda.t()*diagmat(1/psi)*Lambda + eye(k, k));
  //arma::mat Eigvec;
  //arma::vec eigval;
  //eig_sym(eigval, Eigvec, Lambda.t() * diagmat(1/psi) * Lambda);
  
  arma::mat F1 = (Y*diagmat(1/psi)*Lambda)*G1;
  
  arma::mat UG, VG, svdG;
  arma::vec sG;
  svd(UG, sG, VG, G1, "standard");
  svdG = diagmat(sqrt(sG))*VG;
  arma::mat Factors = F1 + randn(T, k)*svdG;
  return Factors;
}