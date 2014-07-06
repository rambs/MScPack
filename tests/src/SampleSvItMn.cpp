/***
 * MScPack
 * Description: Sampling multivariate stochastic volatility
 * Author: Rafael Barcellos
 * Last updated 2nd July, 2014
 * R 3.1.0
 */


//[[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;

// calling R function
Environment MScPack("package:MScPack");
Function aproxMN = MScPack["GetLogX2ApproxParms"]; 

//' @title Sample stochastic volatilities
//' @description Given the parameters of the autoregressive process of the 
//'   log-volatilities, this function uses the FFBS algorithm to get a sample 
//'   from the posterior distribution.
//' @param y data matrix which follow the SV model.
// [[Rcpp::export]]
arma::mat SampleSvItMn(arma::mat y, arma::mat eta, arma::vec mu, arma::vec phi,
  arma::mat U)
  {
    // This function sample the stochastic volatilities according to the results
    //   found by Kim, Shephard and Chib (1998).
    // The difference is that it considers the mixture of normals approximation
    //   presented in Omori et al. (2007).
    
    
    // y is the data matrix which follows the SV model.
    int r = y.n_cols; // y must enter T \times r (time by variables)
    int T = y.n_rows;
    
    // eta is the matrix of log-variances entered as r \times T (variables by 
    //   time).
    if (eta.n_rows != y.n_cols || eta.n_cols != y.n_rows){
      throw std::range_error("Arguments 'y' and 'eta' mismatched.");
    }

    arma::mat ystar = log(y % y);

    Rcpp::List parmsMN = as<Rcpp::List>(aproxMN(arma::trans(ystar.row(0))-eta.col(0)));

    arma::vec M = parmsMN["M"];
    arma::vec V2 = parmsMN["S"];
    
    arma::vec logx2 = arma::trans(ystar.row(0)) - M;
    arma::vec Ew = (eye(r, r) - arma::diagmat(phi)) * mu;
    
    // Forward filter
    
    // F is an identity matrix
    arma::mat G = arma::diagmat(phi);
    arma::mat V = arma::diagmat(V2);
    arma::mat invV = arma::diagmat(1/V2);

    arma::mat m(r, T);
    arma::cube C(r, r, T);
    arma::mat a(r, T);
    arma::cube R(r, r, T);
    arma::mat rootR, rootInvR;
    arma::cube invR(r, r, T);
    invR.zeros();
    arma::vec f(r);
    arma::mat Q(r, r);
    arma::mat A;
    arma::vec e;

    // Note that the volatilities start at t = 1, with prior distribution for
    //  eta_1|D_0 following the VAR(1) stationary distribution. 
    // So, to start the forward filter we must calculate the posterior of 
    //  eta_1|D_1.
    
    arma::mat W = U/(1-phi*phi.t());
    arma::mat invW = arma::inv(arma::symmatu(W));
    C.slice(0) = arma::inv(arma::symmatu(invV + invW));
    m.col(0) = C.slice(0) * (invV * logx2 + invW * mu); // centered: mean = 0
    arma::mat rootC = arma::chol(arma::symmatu(C.slice(0)));
    
    // objects needed for the backward sampling
    arma::mat L, Eigvec, S, rootH;
    arma::vec eigval;
    arma::mat invU = arma::inv(arma::symmatu(U));
    
    for (int k = 1; k < T; k++){
      // quantities for backward sampling
      L = rootC * G;
      arma::eig_sym(eigval, Eigvec, L * invU * L.t());
      S = arma::diagmat(1.0/(eigval + 1.0));
      rootH = arma::diagmat(sqrt(S.diag())) * Eigvec.t() * rootC;
      eta.col(k-1) = rootH.t() * randn(r);
      
      // mixture of normals approximation
      Rcpp::List parmsMN = as<Rcpp::List>(aproxMN(arma::trans(ystar.row(k))-eta.col(k)));
      arma::vec M = parmsMN["M"];
      arma::vec V2 = parmsMN["S"];
      arma::vec logx2 = arma::trans(ystar.row(k)) - M;
    
      // evolution
      a.col(k) = G * m.col(k-1) + Ew; // see W&H (1997), p. 583
      R.slice(k) = symmatu(G * C.slice(k-1) * G.t() + U);
      rootR = arma::chol(R.slice(k));
      rootInvR = arma::trans(arma::inv(arma::trimatu(rootR)));
      invR.slice(k) = rootInvR.t() * rootInvR;
      
      // prediction
      f = a.col(k);
      e = logx2 - f; // algorithm starts at t = 2
      
      V = arma::diagmat(V2);      
      invV = arma::diagmat(1/V2);
      // L = rootR * F.t() = rootR 
      // ==> L * invV * L.t() = rootR * invV * rootR.t()
      arma::eig_sym(eigval, Eigvec, rootR * invV * rootR.t());
      arma::mat S = arma::diagmat(1.0/(eigval + 1.0));
      
      arma::mat rootC = arma::diagmat(sqrt(S.diag())) * Eigvec.t() * rootR;
      C.slice(k) = rootC.t() * rootC;
      
      A = C.slice(k) * invV;
      m.col(k) = a.col(k) + A*e;
      
    }

    // Backward sampling
    
    eta.col(T-1) = m.col(T-1) + rootC.t() * randn(r);
    arma::vec h;
  
    for (int k = 1; k < T; k++){
      h = m.col(T-1-k) + 
        C.slice(T-1-k) * G * invR.slice(T-k) * (eta.col(T-k) - a.col(T-k));
      eta.col(T-1-k) += h;
    }

    // adding the mean
    //eta += repmat(mu, 1, T);

    return eta;
  }