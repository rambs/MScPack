#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;
// funcao FFBS para regressao com V constante e fator de desconto
// CC e' calculada via fatoracao de Choleski
Rcpp::List drmFFBSdiscW_chol(arma::mat Y, arma::mat X, arma::mat dV, double discW, 
arma::vec m0, arma::mat C0_chol, arma::mat C0Inv){

// quantidades para iteracoes  
int T = Y.n_rows; // ATENCAO: matriz de dados deve entrar T x q!!
int q = Y.n_cols; 
int r = q*X.n_cols; // ATENCAO: matriz de exogenas deve entrar sendo T x p

arma::mat I_q = eye(q, q);

// FORWARD FILTERING
arma::mat mm(r, T+1);
arma::mat CC(r, r);
arma::mat CCInv(r, r);
arma::mat rootCC(r, r);
arma::mat rootCCInv(r, r);

arma::vec aa(r);//nao e' preciso armazenar aa_t
arma::mat RR(r, r);//nao e' preciso armazenar RR_t

arma::vec ff(q); //nao e' preciso armazenar ff
arma::mat QQ(q, q); //nao e' preciso armazenar seus valores

arma::mat AA;
arma::mat ee;
arma::mat FF;

//inversa da matriz de covariancias
arma::mat dVInv = inv(diagmat(dV));

// constante para a suavizacao
double sqrt1mDiscW = sqrt(1-discW);

//valores aleatorios gerados
arma::mat th = randn(r, T+1);

mm.col(0) = m0;
rootCC = C0_chol;
CC = rootCC.t() * rootCC; // faltava inicializacao de CC. 
CCInv = C0Inv;
th.col(0) = sqrt1mDiscW * rootCC.t() * th.col(0);

for (int k = 1; k < (T+1); k++){
  //evolucao
  aa = mm.col(k-1);
  RR = CC/discW;
  
  //predicao
  FF = arma::kron(I_q, X.row(k-1)).t();
  ff = FF.t()*aa;
  QQ = arma::symmatu(FF.t() * RR * FF + dV);

  //atualizacao
  CCInv = FF * dVInv * FF.t() + discW*CCInv; 
  rootCCInv = arma::trimatu(arma::chol(arma::symmatu(CCInv)));
  rootCC = arma::trans(arma::inv(rootCCInv));
  CC = rootCC.t()*rootCC;
  
  AA = CC*FF*dVInv; //p.105 West & Harrison (1997)
  ee = arma::trans(Y.row(k-1)) - ff;
  mm.col(k) = aa + AA*ee;
  
  //simulacao retrospectiva
  th.col(k) = sqrt1mDiscW*rootCC.t()*th.col(k);
  mm.col(k-1) -= discW*mm.col(k-1); //mm agora representa mm*
}

// BACKWARD SAMPLING
arma::mat Mu(q, T);
th.col(T) += mm.col(T);
for(int k = 0; k < T; k++){
  th.col(T-1-k) += mm.col(T-1-k) + discW*th.col(T-k);
  Mu.col(T-1-k) = arma::kron(I_q, X.row(T-1-k))*th.col(T-k);
}

return Rcpp::List::create(Named("th") = th, Named("Mu") = Mu.t(), 
Named("ZW") = sqrt((1-discW)/discW)*rootCC);
}
