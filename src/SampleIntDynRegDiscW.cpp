#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;
//amostra intercepto e parametros de estado

Rcpp::List SampleIntDynRegDiscW(arma::mat Y, arma::mat X, arma::mat dV, 
  double discW, arma::vec m0, arma::mat C0, arma::vec a0, arma::mat A0){

// quantidades para iteracoes  
int T = Y.n_rows; // ATENCAO: matriz de dados deve entrar T x q!!
int q = Y.n_cols; 
int r = q*X.n_cols; // ATENCAO: matriz de exogenas deve entrar sendo T x p

arma::mat I_q = eye(q, q);
arma::mat I_2 = eye(2, 2);

// FORWARD FILTERING
arma::mat mm(q+r, T+1); //vetor m contem intercepto e parms de estado
arma::cube CC(q+r, q+r, T+1); //matriz C e' cov. intercepto + parms de estado
arma::mat CCInv(q+r, q+r);
arma::mat rootCC(q+r, q+r);
arma::mat rootCCInv(q+r, q+r);

arma::mat aa(q + r, T + 1);
arma::mat RR(q + r, q + r);
arma::cube RRInv(q + r, q + r, T + 1);
arma::mat rootRRInv(q + r, q + r);

arma::vec ff(q); //nao e' preciso armazenar ff

arma::mat AA;
arma::vec ee;
arma::mat FF(q + r, q);
FF.submat(r, 0, r + q - 1, q - 1) = I_q;

//inversa da matriz de covariancias
arma::mat dVInv = arma::inv(arma::diagmat(dV));

//valores iniciais de m e C
mm.col(0) = arma::join_cols(m0, a0);
CC.slice(0) = arma::join_cols(arma::kron(I_2.row(0), C0), arma::kron(I_2.row(1), A0));

//std::cout << "\nRodou ate antes do forward filter\n";
for (int k = 1; k < (T+1); k++){
  //evolucao
  aa.col(k) = mm.col(k - 1);
  RR = arma::symmatu(CC.slice(k - 1));
  RR.submat(0, 0, r - 1, r - 1) += (1 - discW)/discW * RR.submat(0, 0, r - 1, r - 1);
  
  //predicao
  FF.submat(0, 0, r - 1, r - 1) = arma::kron(I_q, arma::trans(X.row(k - 1)));
  ff = FF.t()*aa.col(k);
  //nao e' preciso calcular QQ

  //atualizacao
  rootRRInv = arma::trans(arma::inv(arma::trimatu(arma::chol(RR))));
  RRInv.slice(k) = rootRRInv.t() * rootRRInv;
  CCInv = arma::symmatu(FF * dVInv * FF.t() + RRInv.slice(k)); 
  rootCCInv = arma::trimatu(arma::chol(CCInv));
  rootCC = arma::trans(arma::inv(rootCCInv));
  CC.slice(k) = rootCC.t() * rootCC;
  
  AA = CC.slice(k)*FF*dVInv; //p.105 West & Harrison (1997)
  ee = arma::trans(Y.row(k-1)) - ff;
  mm.col(k) = aa.col(k) + AA * ee;
  
}

//std::cout << "\n Rodou ate antes do backward sampler\n";
// BACKWARD SAMPLING
//valores aleatorios gerados
arma::vec alpha = arma::randn(q);
arma::mat th = arma::randn(r, T + 1);

alpha = mm.submat(r, T, q+r-1, T) + 
  arma::trans(rootCC.cols(r, q+r-1)) * arma::join_cols(th.col(T), alpha);
  
th.col(T) = mm.submat(0, T, r-1, T) + 
  arma::trans(rootCC.submat(0, 0, r-1, r-1)) * th.col(T);

arma::mat Mu(q, T); 
Mu = repmat(alpha, 1, T);

arma::vec hh;

for(int k = 0; k < T; k++){
  arma::mat HHAlpha = CC.slice(T-1-k) - CC.slice(T-1-k) * RRInv.slice(T-k) * CC.slice(T-1-k);
  arma::mat rootHH = arma::chol(symmatu(HHAlpha.submat(0, 0, r-1, r-1)));
  //matriz para amostrar somente state parms
  arma::mat P = CC.slice(T-1-k) * RRInv.slice(T-k);
  hh = mm.col(T-1-k) + 
    P * (arma::join_cols(th.col(T-k), alpha) - aa.col(T-k));

  th.col(T-1-k) = hh.subvec(0, r-1) + rootHH.t() * th.col(T-1-k);
  FF = arma::kron(I_q, arma::trans(X.row(T-1-k)));
  Mu.col(T-1-k) += FF.t() * th.col(T-k);
}

return Rcpp::List::create(Named("alpha") = alpha, Named("th") = th, Named("Mu") = Mu.t(), 
Named("rootW") = sqrt((1-discW)/discW)*rootCC.submat(0, 0, r - 1, r - 1));
}
