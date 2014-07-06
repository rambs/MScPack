/***
 * MScPack
 * Description: Sampling SV VAR parameters
 * Author: Rafael Barcellos
 * Last updated 5th July, 2014
 * R 3.1.0
 */

//[[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;

// calling 'rtruncnorm'
Environment truncnorm("package:truncnorm");
Function rtruncnormR = truncnorm["rtruncnorm"];

// auxiliary function for Metropolis  
double logC(arma::vec phi, arma::vec e, arma::mat U){
  arma::mat WW = U/(1 - phi*phi.t());
  double logDetWW;
  double sign;
  log_det(logDetWW, sign, WW);
  arma::vec x = -0.5*logDetWW - 0.5 * e.t() * arma::inv(arma::symmatu(WW)) * e;
  return x(0);
}

//' @title Sample SV VAR parameters
//' @description Given the stochastic volatilities and other parameters, this 
//'   function sample the autoregressive parameters that defines the evolution
//'   of the log-volatilities. For further information, see the details below.
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
arma::vec SampleSvPhi(arma::vec phi, arma::mat eta, arma::vec mu, arma::mat U){
  int r = eta.n_rows;
  int T = eta.n_cols;
  arma::mat E = eta - repmat(mu, 1, T);

  arma::mat Uinv = arma::inv(U);
  arma::mat Binv = Uinv % (E.cols(0, T-2) * arma::trans(E.cols(0, T-2)));
  arma::mat B = arma::inv(Binv);
  arma::mat Binv_b = sum(E.cols(0, T-2) % (Uinv * E.cols(1, T-1)), 1);
  arma::vec b = B * Binv_b;

  arma::vec phiOut = phi;
  arma::vec muCond, sigCond;
  arma::vec xCond = phi;

  // vector of indices
  Rcpp::IntegerVector ind = seq_len(r)-1;
  arma::uvec indices = as<arma::uvec>(ind);
  arma::uvec eqI;
  arma::uvec neqI;

  // objects for Metropolis
  double log_a, log_u;
  Rcpp::NumericVector MHlog(2, 0.0);
  
  for (int i = 0; i < r; i++){
    eqI = find(indices == i);
    neqI = find(indices != i);
    muCond = b(i) + B.submat(eqI, neqI) * arma::inv(B.submat(neqI, neqI)) *
      (xCond.elem(neqI) - b.elem(neqI));
    sigCond = B.submat(eqI, eqI) - B.submat(eqI,neqI) * 
      arma::inv(B.submat(neqI, neqI)) * B.submat(neqI, eqI);
    double x = as<double>(rtruncnormR(1, 0.0, 1.0, muCond, sqrt(sigCond))); 
    xCond(i) = x;

    MHlog(1) = logC(xCond, E.col(0), U) - logC(phiOut, E.col(0), U);
    log_a = min(MHlog);
    log_u = log(randu(1))[0];
    if(log_u <= log_a){
      phiOut(i) = x;
      //n.aceita[i] = n.aceita[i]+1
    } else {
      xCond(i) = phi(i);
    }
  }
  return phiOut;
}
