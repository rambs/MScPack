/**
 * MScPack
 * Description: Sampling loadings in vectorized form
 * Author: Rafael Barcellos
 * Last updated 4 June 2014
 * R 3.0.2
 */
 
//[[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;

// this function sample the loading parameters considering constant 
//  observational variance.
// [[Rcpp::export]]
arma::mat SampleDfmLoads(arma::mat y, arma::mat factors, arma::vec psi, 
  int s, double c0)
{
  /**
   * The data matrix 'y' is (T \times q).
   * Note that the factors matrix is (T+h \times k). In the models considered
   *  in the dissertation, h > s. This means that we always have the 
   *  VAR order, 'h', greater than the number of dynamic factors, 's', in the 
   *  observational equation. 
   * The argument 's' is the order of the dynamic factors present in the
   *  observational equation.
   * The prior mean is set to zero.
   */
   
  int T = y.n_rows;
  int q = y.n_cols;;
  int k = factors.n_cols;
  int TpH = factors.n_rows;
  int h = TpH - T;
  if (h < s+1){
    throw std::range_error("Argument 's' greater than or equal to factors' VAR order in SampleDfmLoads.cpp");
  }
  arma::mat I_q = eye(q, q);
//  std::cout << "\nh = " << h << "\n";
  arma::mat Xlambda;
  Xlambda = factors.rows(h, TpH-1);
//  std::cout << "\nNumber of rows in Xlambda is " << Xlambda.n_rows << "\n";
  if (s > 0){
    for (int i = 1; i <= s; i++){
    Xlambda = arma::join_rows(Xlambda, factors.rows(h-i, TpH-i-1));
    }
  }
  
  // posterior variance
  arma::mat Eigvec;
  arma::vec eigval;
  
  arma::eig_sym(eigval, Eigvec, arma::symmatu(Xlambda.t() * Xlambda));
  
  arma::mat kronEig = kron(I_q, Eigvec);  
  arma::vec L1eigval = 1.0/(arma::kron((1.0/psi), eigval) + 1/c0);
  arma::mat L1 = kronEig * diagmat(L1eigval) * kronEig.t();
  
  // posterior mean
  arma::vec yStar(y.begin(),  T*q);
  arma::vec l1 = L1 * (arma::kron(arma::diagmat(1.0/psi), Xlambda.t()) * yStar);
  arma::mat l1Mat(l1.begin(), k*(s+1), q); // mean in matrix form
  
  // random generation
  arma::vec vecnorm = sqrt(L1eigval) % randn(q*k*(s+1), 1);  
  arma::mat LambdaBar(vecnorm.begin(), k*(s+1), q, false);
  LambdaBar = l1Mat + Eigvec * LambdaBar;
  
  return LambdaBar.t();
}
