#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;
// funcao FFBS para regressao com V constante e fator de desconto
//[[Rcpp::export(".drmFFBSdiscW")]]
Rcpp::List drmFFBSdiscW(arma::mat Y, arma::mat X, arma::mat dV, double discW, 
arma::vec m0, arma::mat ZC0){
int T = Y.n_rows; // ATENCAO: matriz de dados deve entrar T x q!!
int q = Y.n_cols; 
int r = q*X.n_cols; // ATENCAO: matriz de exogenas deve entrar sendo T x p

mat Yt = Y.t();
mat I_q = diagmat(ones(q));

// FORWARD FILTERING

mat mm(r, T+1);
vec aa(r);//nao e' preciso armazenar aa_t
mat RR(r, r);//nao e' preciso armazenar RR_t
vec ff(q); //nao e' preciso armazenar ff
mat QQ(q, q); //nao e' preciso armazenar seus valores
mat AA;
mat ee;
mat FF;

//para previsao, e' preciso armazenar W_{T+1}
mat dVInv = diagmat(1/dV.diag());

//parametros para forma quadrada
mat ZR;
cube ZC(r, r, T+1); //utilizado para fazer a amostragem retrospectiva

//objetos para decomposicao em valores singulares e espectral
mat U, V, L, Eigvec;
vec s, eigval; 

double sqrtDiscW = sqrt(discW);

mm.col(0) = m0;
ZC.slice(0) = ZC0;

for (int k = 1; k < (T+1); k++){
  //evolucao
  aa = mm.col(k-1);
  ZR = ZC.slice(k-1)/sqrtDiscW;
  RR = ZR * ZR.t();
  
  //predicao
  FF = kron(I_q, X.row(k-1)).t();
  ff = FF.t()*aa;
  QQ = symmatu(FF.t() * RR * FF + dV);

  //atualizacao
  L = ZR.t()*FF; //forma raiz quadrada
  eig_sym(eigval, Eigvec, L * dVInv * L.t());
  ZC.slice(k) = ZR * Eigvec * diagmat(1.0/sqrt(eigval + 1.0));
  
  //AA = RR * FF * inv(QQ); //p.104 West & Harrison (1997)
  AA = ZC.slice(k)*(ZC.slice(k)).t()*FF*dVInv; //p.105 West & Harrison (1997)
  ee = Yt.col(k-1) - ff;
  
  mm.col(k) = aa + AA*ee;
}

// BACKWARD SAMPLING

//valores aleatorios gerados
mat th = randn(r, T+1);

//simulacao via forma raiz quadrada
th.col(T) = mm.col(T) + ZC.slice(T) * th.col(T);

double sqrt1mDiscW = sqrt(1-discW);
mat hh;
mat ZH;
mat Mu(q, T);
for (int k = 0; k < T; k++){
  //inv(RR_t) = delta*inv(CC_{t-1}) 
  //==> CC.slice(T-1-k) * RRinv.slice(T-k) = delta * I_n
  //aa_t = mm_{t-1}
  hh = (1-discW)*mm.col(T-1-k) + discW * th.col(T-k);
  ZH = sqrt1mDiscW*ZC.slice(T-1-k); //forma raiz quadrada
  th.col(T-1-k) = hh + ZH*th.col(T-1-k);
  Mu.col(T-1-k) = kron(I_q, X.row(T-1-k))*th.col(T-k);
}

return Rcpp::List::create(Named("th") = th, Named("Mu") = Mu.t(), 
Named("ZW") = sqrt((1-discW)/discW)*ZC.slice(T));
}
