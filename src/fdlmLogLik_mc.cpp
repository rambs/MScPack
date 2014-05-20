
#include <RcppArmadillo.h>
#include <omp.h>
using namespace Rcpp;

const double log2pi = std::log(2.0 * M_PI);

//[[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::plugins("openmp")]]
// [[Rcpp::export]]
arma::vec fdlmLogLik_mc(arma::mat Y, arma::mat xFix, arma::mat xDyn, double nFactors, 
Rcpp::List hyperparms, 
arma::mat Y_NA, IntegerVector whichNA, NumericVector Beta, NumericVector Lambda, 
arma::mat psi, double discW, int cores)
{
  omp_set_num_threads(cores);
  int T = Y.n_rows;
  int q = Y.n_cols;
  int pDyn = xDyn.n_cols;
  int pFix = xFix.n_cols;
  int r = pDyn*q;
  int N = psi.n_rows;
  
  arma::cube Beta_(Beta.begin(), pFix, q, N, false);
  arma::cube Lambda_(Lambda.begin(), q, nFactors, N, false);
  std::cout << "\nrodou ate Beta_ e Lambda_\n";
  
  double constants = -(q/2)*log2pi;
  arma::vec nloglik(N);
  
  arma::mat V(q, q);
  arma::mat VInv(q, q);
  std::cout << "\nrodou ate depois de V e VInv\n";
  
  arma::vec mm(r);
  arma::mat CC(r, r);
  std::cout << "\nrodou ate depois de mm e CC.\n";
  arma::vec aa(r);
  arma::mat RR(r, r);
  arma::vec ff(q);
  arma::mat QQ(q, q);
  std::cout << "\nrodou ate depois de ff e QQ.\n";
  arma::mat FF;
  arma::vec ee;
  arma::mat AA;
  arma::mat I_q = arma::eye(q, q);
  std::cout << "\nrodou ate depois de diagmat.\n";
  
  //arma::vec z;
  //double rootisum;
  arma::mat Eigvec;
  arma::vec eigval;
  mm = as<arma::vec>(hyperparms["m0"]);
  CC = as<arma::mat>(hyperparms["C0"]);
  
  
  arma::mat Yn(T, q);
  arma::mat En(T, q);
  Yn = Y;
  std::cout << "\nrodou ate antes do for.\n";
   
  for(int n = 0; n < N; n++){
    double loglik = 0.0;
    
    if(whichNA.size()>0){
      Yn.elem(as<arma::uvec>(whichNA)) = arma::trans(Y_NA.row(n));  
    }
    
    En = Yn - xFix * Beta_.slice(n);
    V = Lambda_.slice(n) * arma::trans(Lambda_.slice(n)) + arma::diagmat(psi.row(n));
    VInv = arma::inv(symmatu(V));
    
    #pragma omp parallel for schedule(static)    
    for(int k = 1; k < T+1; k++){
    //evolucao
    aa = mm;
    RR = CC/discW;
    
    //predicao
    FF = kron(I_q, xDyn.row(k-1));
    ff = arma::trans(FF) * aa;
    QQ = arma::symmatu(arma::trans(FF) * RR * FF + V);
    
    //verossimilhanca
    arma::mat rooti = arma::trans(arma::inv(trimatu(arma::chol(QQ))));
    double rootisum = arma::sum(log(rooti.diag()));
    arma::vec z = rooti * arma::trans( En.row(k-1) - ff.t()) ;    
    loglik += constants - 0.5 * arma::sum(z%z) + rootisum; 
    /*
    arma::eig_sym(eigval, Eigvec, QQ);
    rootisum = arma::sum(0.5*log(1.0/eigval));
    z = arma::trans( (En.row(k-1) - ff.t()) * arma::diagmat(1.0/eigval)*arma::trans(Eigvec) );
    loglik += constants - 0.5*arma::sum(z%z) + rootisum;
    */
    //atualizacao
    //CC = RR - RR * FF * arma::inv(QQ) * arma::trans(FF) * RR;
    CC = arma::inv(FF * VInv * arma::trans(FF) + arma::inv(RR));
    //arma::eig_sym(eigval, Eigvec, FF * VInv * arma::trans(FF));
    AA = CC * FF * VInv;
    mm = aa + AA*(arma::trans(En.row(k-1) - ff.t()));
    }
    nloglik(n) = loglik;
  }
  return nloglik;
}
