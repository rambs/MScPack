/**
 * MScPack
 * Description: Sampling dynamic factors with variance being identity
 * Author: Rafael Barcellos
 * Last updated 3 June 2014
 * R 3.0.2
 */
 
//[[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;

// It's important to note that the prior distribution of every factor before 
//  t = 1 was chosen so as to avoid the use of Metropolis step within the Gibbs 
//  sampler.
// A simple way to do that is to consider those factors following a normal
//  distribution with mean zero and variance equals to identity.

// [[Rcpp::export]]
arma::mat SampleDynFactors(arma::mat Y, arma::mat LambdaBar, arma::mat PhiBar,
  arma::vec psi)
  {
    int T = Y.n_rows; //Y is (T \times q)
    int q = Y.n_cols;
    int k = PhiBar.n_rows; // PhiBar is (k \times k*h)
    int h = PhiBar.n_cols/k; // PhiBar must enter correctly
    int s = LambdaBar.n_cols/k - 1; //LambdaBar is (q \times k*(s+1))
    int r = k*h; // vector length of state parms
    
    // definition of FLambda
    arma::mat FLambda;
    if (h == s+1){
      FLambda = LambdaBar;
    } else {
      FLambda = arma::join_rows(LambdaBar, zeros(q, k*(h-s-1)));
    }
    
    // definition of GPhi
    arma::mat GPhi;
    if (h == 1){
      GPhi = PhiBar;
    } else {
      GPhi = join_cols(PhiBar, eye(r-k, r));
    }
    
    // definition of V and InvV
    arma::mat V = arma::diagmat(psi);
    arma::mat InvV = arma::diagmat(1.0/psi);
    
    // definition of Wphi
    arma::mat Wphi(r, r);
    Wphi.zeros();
    Wphi.submat(0, 0, k-1, k-1) = eye(k, k);
    
    arma::vec a;
    arma::mat R;
    arma::mat rootR;
    
    arma::vec f;
    arma::vec e;
    arma::mat A;
    
    arma::mat m(r, T+1);
    arma::cube C(r, r, T+1);
    arma::mat rootC;
    
    // definition of prior hyperparms
    m.zeros();
    C.slice(0) = 1e2 * eye(r, r);
    
    // Forward filter
    // The filter is done considering the factorisation in square root form.
    // To do so we need to specify:
    arma::mat L;
    arma::vec eigval;
    arma::mat Eigvec;
    
    for (int t = 1; t < T+1; t++){
      
      a = GPhi * m.col(t-1);
      R = arma::symmatu(GPhi * C.slice(t-1) * GPhi.t() + Wphi);
      rootR = arma::chol(R);
      //std::cout << "\nCholeski factorisation of R at iteration " << t << " done.\n";
            
      f = FLambda * a;
      e = arma::trans(Y.row(t-1)) - f;
      
      L = rootR * FLambda.t();
      eig_sym(eigval, Eigvec, L * InvV * L.t());
      
      arma::mat S = arma::diagmat(1.0/(eigval + 1.0));
      
      rootC = arma::diagmat(sqrt(S.diag())) * Eigvec.t() * rootR;
      C.slice(t) = rootC.t() * rootC;
      
      A = C.slice(t) * FLambda.t() * InvV;
      
      m.col(t) = a + A * e;
      
    }
    
    // Backward sampling
    
    // generating normal random variables
    arma::vec phi = m.col(T) + rootC.t() * randn(r, 1);
    arma::mat phiMat(phi.begin(), k, h, false);
    
    // reordering the matrix
    IntegerVector colOrder = seq_len(h)-1;
    std::reverse(colOrder.begin(), colOrder.end());
    // it must be checked if this works fine
    // std::reverse checked in MScPackExamples/examples/rcppRev.R and worked ok.
    
    // since the factos before t=1 are N(0,I), they must be sampled backward
    arma::mat factors(k, T+h); 
    factors.cols(T, T+h-1) = phiMat.cols(as<arma::uvec>(colOrder));
    
    arma::mat C_t, C11, C22, C12, InvC11, C1, InvC1, rootInvHStar, rootHStar;
    arma::vec m1, m2, hStar;
    
    arma::mat PhiM, PhiH;
    if (h == 1){
      PhiH = PhiBar;
    } else {
      PhiM = arma::join_rows(eye(k, k), -PhiBar.cols(0, r-k-1));
      PhiH = PhiBar.cols(r-k, r-1);
    }
    
    if (r == k){
      for (int t = T-1; t >= 0; t--){
        C1 = C.slice(t);
        m1 = m.col(t);
        
        InvC1 = arma::inv(arma::symmatu(C1));
        rootInvHStar = arma::chol(arma::symmatu(InvC1 + PhiH.t() * PhiH));
        rootHStar = arma::trans(arma::inv(arma::trimatu(rootInvHStar)));
        m2 = factors.col(t+1);
        hStar = rootHStar.t() * rootHStar * (InvC1 * m1 + PhiH.t() * m2);
        factors.col(t) = hStar + rootHStar.t() * randn<vec>(k);
      }
    }
    
    if (r > k){
      for (int t = T-1; t >= 0; t--){
        C_t = C.slice(t);
        C11 = C_t.submat(0, 0, r-k-1, r-k-1);
        C22 = C_t.submat(r-k, r-k, r-1, r-1);
        C12 = C_t.submat(0, r-k, r-k-1, r-1);
        InvC11 = arma::inv(C11);
        
        C1 = C22 - C12.t() * InvC11 * C12;
        m1 = m.submat(r-k, t, r-1, t) + 
          C12.t() * InvC11 *(phi.subvec(k, r-1) - m.submat(0, t, r-k-1, t));
        
        InvC1 = arma::inv(arma::symmatu(C1));
        
        m2 = PhiM * phi;
        rootInvHStar = arma::chol(arma::symmatu(InvC1 + PhiH.t() * PhiH));
        rootHStar = arma::trans(arma::inv(arma::trimatu(rootInvHStar)));      
        hStar = rootHStar.t() * rootHStar * (InvC1 * m1 + PhiH.t() * m2);
        
        factors.col(t) = hStar + rootHStar.t() * randn<vec>(k);
        
        phi = arma::join_cols(phi.subvec(k, r-1), factors.col(t));
      }
    }
    return factors;
      // Is it possible to do this conditioning faster?
      // Maybe using the results presented in Bai and Wang (2012) p.47
      
  }
