/**
 * MScPack
 * Description: Sampling VAR(h) parms
 * Author: Rafael Barcellos
 * Last updated 4 June 2014
 * R 3.0.2
 */
 
//[[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;

// this function sample the VAR(h) parameters considering the covariance evolution matrix
//  being identity.
// [[Rcpp::export]]
arma::mat SampleVarParms(arma::mat Y, int h, double c0)
{
  /** The argument 'h' means the order of the VAR(h).
   *  The argument 'c0' is the prior variance, that's a constant times the 
   *   identity matrix.
   *  The prior mean is set to zero for all VAR parameters.
   *
   */
   int T = Y.n_rows;
   int k = Y.n_cols;
   arma::mat Yphi = Y.rows(h, T-1);
   arma::mat Xphi;
   Xphi = Y.rows(h-1, T-2);
   if (h > 1){
     for (int i = 2; i <= h; i++){
     Xphi = arma::join_rows(Xphi, Y.rows(h-i, T-i-1));
     }
   }
   arma::mat I_k = eye(k, k);
   arma::mat Eigvec;
   arma::vec eigval;
   
   eig_sym(eigval, Eigvec, Xphi.t() * Xphi);   
   arma::mat Z = diagmat(sqrt(c0/(c0 * eigval + 1.0))) * Eigvec.t();
   
   arma::vec yStar(Yphi.begin(), (T-h)*k, false);
   arma::vec g1 = (kron(I_k, Z.t() * Z * Xphi.t())*yStar);
   
   arma::mat gBar(g1.begin(), h*k, k, false);
   arma::mat PhiBar = gBar.t() + randn(k, h*k) * Z;
   return PhiBar;  
}
