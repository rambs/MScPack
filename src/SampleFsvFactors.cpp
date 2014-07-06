/***
 * MScPack
 * Description: Sampling factors from a FSV
 * Author: Rafael Barcellos
 * Last updated 6th July, 2014
 * R 3.1.0
 */


//[[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;

//' @title Sample factors of FSV model
//' @description Given the data and the log-volatilities, this function sample
//'   the factors from the posterior distribution.
//' @param y data matrix which follow the FSV model.
//' @param Lambda loadings matrix.
//' @param eta matrix of log-volatilities (\eqn{k \times T}).
//' @param psi vector of idiosyncratic variances.
//' @return A matrix of factors (\eqn{k \times T}).
// [[Rcpp::export]]
arma::mat SampleFsvFactors(arma::mat y, arma::mat Lambda, arma::mat eta, 
  arma::vec psi)
  {
    // This function sample the factors from a FSV model.    
    
    // y is the data matrix which follows the SV model.
    // y must enter T \times q (time by variables)
    int T = y.n_rows;
    int r = eta.n_rows;
    
    // eta is the matrix of log-variances entered as r \times T (variables by 
    //   time).
    if (eta.n_cols != y.n_rows){
      throw std::range_error("Arguments 'y' and 'eta' mismatched.");
    }
    
    arma::mat factors(r, T);
    factors.randn();
    arma::mat invPsi = arma::diagmat(1/psi);
    arma::mat LtInvPsiL = Lambda.t() * invPsi * Lambda;
    for (int t = 0; t < T; t++){
      arma::mat invH0 = arma::diagmat(1/exp(eta.col(t)));
      arma::mat invH = LtInvPsiL + invH0;
      arma::mat rootInvH = arma::chol(arma::symmatu(invH));
      arma::mat rootH = arma::trans(arma::inv(arma::trimatu(rootInvH)));
      arma::vec f1 = rootH.t() * rootH * 
        (Lambda.t()*invPsi*arma::trans(y.row(t)));
      factors.col(t) = f1 + rootH.t() * factors.col(t);
    }
    
    return factors;
  }