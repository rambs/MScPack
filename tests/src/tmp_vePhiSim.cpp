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

// auxiliary function for Metropolis  
double logC(arma::vec phi, arma::vec e, arma::mat U){
  arma::mat WW = U/(1 - phi*phi.t());
  double logDetWW;
  double sign;
  log_det(logDetWW, sign, WW);
  arma::vec x = -0.5*logDetWW - 0.5 * e.t() * arma::inv(arma::symmatu(WW)) * e;
  return x(0);
}

// [[Rcpp::export]]
vec vePhiSim(vec phi, mat eta, vec mu, mat U){
  int r = eta.n_rows;
  int T = eta.n_cols;
mat E = eta - repmat(mu, 1, T);

mat Uinv = inv(U);
mat Binv = Uinv%(E.cols(0, T-2)*E.cols(0, T-2).t());
mat B = inv(Binv);
mat Binv_b = sum(E.cols(0, T-2)%(Uinv*E.cols(1, T-1)), 1);
vec b = B*Binv_b;

vec phiOut = phi;
vec muCond, sigCond;
vec xCond = phi;

// vetor de indices
Rcpp::IntegerVector ind = seq_len(r)-1;
uvec indices = as<arma::uvec>(ind);
uvec eqI;
uvec neqI;

// argumentos para Metropolis
double log_a, log_u;

// funcao necessaria para simular normal truncada
Environment truncnorm("package:truncnorm");
Function rtruncnormR = truncnorm["rtruncnorm"];
Rcpp::NumericVector MHlog(2, 0.0);

for (int i=0; i<r; i++){
  eqI = find(indices==i);
  neqI = find(indices!=i);
  muCond = b(i) + B.submat(eqI, neqI)*inv(B.submat(neqI, neqI))*
           (xCond.elem(neqI) - b.elem(neqI));
  sigCond = B.submat(eqI, eqI) - B.submat(eqI,neqI) * inv(B.submat(neqI, neqI)) * 
            B.submat(neqI, eqI);
  double x = as<double>(rtruncnormR(1, 0.0, 1.0, muCond, sqrt(sigCond))); 
  xCond(i) = x;

  MHlog(1) = logC(xCond, E.col(0), U)-logC(phiOut, E.col(0), U);
  log_a = min(MHlog);
  log_u = log(randu(1))[0];
  if(log_u<=log_a){
  phiOut(i) = x;
  //n.aceita[i] = n.aceita[i]+1
  } else {
  xCond(i) = phi(i);
  }
}
return phiOut;
}
