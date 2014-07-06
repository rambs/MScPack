/***
 * MScPack
 * Description: Sampling SV mean
 * Author: Rafael Barcellos
 * Last updated 5th July, 2014
 * R 3.1.0
 */

//[[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;

//' @title Sample SV mean
//' @description Given the stochastic volatilities and other parameters, this 
//'   function sample the mean parameter from the posterior distribution.
//' @param phi starting values for the Metropolis-Hastings algorithm.
//' @param eta matrix of log-volatilities in the form \eqn{r \times T} 
//'   (variable by time).
//' @param mu mean of the autoregressive process.
//' @param U covariance matrix of the autoregressive process.
//' @details The stochastic volatility model considers the evolution equation 
//'   such as \eqn{\eta[t] - \mu = \phi(\eta[t-1] - \mu) + \omega[t]}, where 
//'   \eqn{0 < \phi < 1} to guarantee stationarity. To sample \eqn{\phi} it is 
//'   needed to use a Metropolis-Hastings algorithm because the prior 
//'   distribution of \eqn{\eta[1]} is N(\mu, U/(1-\phi\phi')). 
//' @return A vector of values sampled.
// [[Rcpp::export]]
arma::vec SampleSvMean(arma::mat eta, arma::vec phi, arma::mat U, 
  double m0, double C0){
  
  int T = eta.n_cols;
  int r = eta.n_rows;
  arma::mat Phi = arma::diagmat(phi);
  arma::mat I_r = arma::eye(r, r);
  
  arma::mat Gamma = eta.cols(1, T-1) - Phi * eta.cols(0, T-2);  
  
  arma::mat W = U/(1 - phi*phi.t());
  arma::mat invW = arma::inv(arma::symmatu(W));  
  arma::mat invU = arma::inv(arma::symmatu(U));
  
  arma::mat invC1 = arma::eye(r, r)/C0 + invW + 
    (T-1) * (I_r - Phi) * invU * (I_r - Phi);
  arma::mat rootInvC1 = arma::chol(arma::symmatu(invC1));
  arma::mat rootC1 = arma::trans(arma::inv(arma::trimatu(rootInvC1)));
  arma::vec m1 = rootC1.t() * rootC1 * (m0/C0 + invW * eta.col(0) + 
    (I_r - Phi) * invU * sum(Gamma, 1));
  
  arma::vec mu = m1 + rootC1.t() * randn(r, 1);
  
  return mu;
}