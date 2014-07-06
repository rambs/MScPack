/***
 * MScPack
 * Description: Sampling SV variance
 * Author: Rafael Barcellos
 * Last updated 5th July, 2014
 * R 3.1.0
 */

//[[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;

// funcao necessaria para simular WI
  Environment MCMCpack("package:MCMCpack");
  Function riwishR = MCMCpack["riwish"];

// auxiliary function for Metropolis  
double logC(arma::vec phi, arma::vec e, arma::mat U){
  arma::mat WW = U/(1 - phi*phi.t());
  double logDetWW;
  double sign;
  log_det(logDetWW, sign, WW);
  arma::vec x = -0.5*logDetWW - 0.5 * e.t() * arma::inv(arma::symmatu(WW)) * e;
  return x(0);
}

//' @title Sample SV variance
//' @description Given the stochastic volatilities and other parameters, this 
//'   function sample the variance parameter from the posterior distribution.
//' @param phi autoregressive parameters.
//' @param eta matrix \eqn{r \times T} (variable by time) of log-volatilities.
//' @param mu mean of the autoregressive process.
//' @details The stochastic volatility model considers the evolution equation 
//'   such as \eqn{\eta[t] - \mu = \phi(\eta[t-1] - \mu) + \omega[t]}, where 
//'   \eqn{0 < \phi < 1} to guarantee stationarity.
//' @return A sampled covariance matrix.
// [[Rcpp::export]]
arma::mat SampleSvVariance(arma::mat U, arma::mat eta, arma::vec phi, 
  arma::vec mu, double v0, arma::mat D0){
  int T = eta.n_cols;
  arma::mat E = eta - repmat(mu, 1, T);

  arma::mat E2 = E.cols(1, T-1) - arma::diagmat(phi) * E.cols(0, T-2);
  arma::mat G = E2*E2.t();
  arma::mat D1 = D0 + G;
  double v1 = v0 + T-1;

  arma::mat UCond = as<arma::mat>(riwishR(v1, D1));
  
// argumentos para Metropolis
//mat UOut;
  Rcpp::NumericVector MHlog(2, 0.0);
  double log_a, log_u;

  MHlog(1) = logC(phi, E.col(0), UCond) - logC(phi, E.col(0), U);
  log_a = min(MHlog);
  log_u = log(randu(1))[0];

  if(log_u <= log_a){
    return UCond;
    //UOut = UCond;
  //n.aceita[i] = n.aceita[i]+1
  } else {
    return U;
    //UOut = U;
  }
}
