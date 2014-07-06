/***
 * MScPack
 * Description: Sampling FSV loads in PLT form
 * Author: Rafael Barcellos
 * Last updated 6th July, 2014
 * R 3.1.0
 */

//[[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;

//' @title Sample FSV loads
//' @description Given the data, the factors and the idiosyncratic variances, 
//'   this function sample the loadings matrix with PLT prior constraints.
//' @param y data matrix which follow the FSV model.
//' @param factors matrix of factors \eqn(T \times k).
//' @param psi idiosyncratic variances.
//' @param m0 prior mean of the loads.
//' @param C0 prior variance of the loads;
//' @return A matrix with the loads.
// [[Rcpp::export]]
arma::mat SampleFsvLoads(arma::mat y, arma::mat factors, arma::vec psi, 
  double m0, double C0){
  int k = factors.n_cols;
  int q = psi.n_rows;
  arma::mat Fi, Ci, invCi, rootInvCi, rootCi;
  arma::vec mi;
  
  // simulacao da matriz com restricoes
  arma::mat Lambda_1(k, k);
  Lambda_1.eye();
  arma::mat FtF = factors.t() * factors;
  
  for (int ik = 1; ik < k; ik++){
    Fi = factors.cols(0, ik-1);
    invCi = (1/psi(ik, 0)) * FtF.submat(0, 0, ik-1, ik-1) + (1/C0);
    rootInvCi = arma::chol(arma::symmatu(invCi));
    rootCi = arma::trans(arma::inv(arma::trimatu(rootInvCi)));
    mi = rootCi.t() * rootCi * ((1/psi(ik, 0)) * Fi.t() * y.col(ik) + (m0/C0));
    Lambda_1.submat(ik, 0, ik, ik-1) = mi.t() + randn(1, ik)*rootCi;  
  }
    
  // simulacao da matriz sem restricoes
  arma::mat Lambda_2(q-k, k);
  Lambda_2.zeros();
  
  for (int ik = k; ik < q; ik++){
    invCi = (1/psi(ik, 0)) * FtF + (1/C0);
    rootInvCi = arma::chol(arma::symmatu(invCi));
    rootCi = arma::trans(arma::inv(arma::trimatu(rootInvCi)));
    mi = rootCi.t() * rootCi * 
      ((1/psi(ik, 0)) * factors.t() * y.col(ik) + (m0/C0));

    Lambda_2.row(ik-k) = mi.t() + randn(1, k)*rootCi;
  }
      
  arma::mat Lambda = join_cols(Lambda_1, Lambda_2); 
  return Lambda;
}
