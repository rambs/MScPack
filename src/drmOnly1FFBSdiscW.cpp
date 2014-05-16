#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;
// funcao FFBS para regressao com V constante e fator de desconto
Rcpp::List drmOnly1FFBSdiscW(arma::mat Y, arma::mat X, arma::mat dV, double discW, 
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
mat dVInv = inv(dV);

//parametros para forma quadrada
mat ZR;
mat ZC(r, r); //utilizado para fazer a amostragem retrospectiva

//objetos para decomposicao em valores singulares e espectral
mat L, Eigvec;
vec eigval; 

double sqrtDiscW = sqrt(discW);
double sqrt1mDiscW = sqrt(1-discW);

//para evitar o loop da amostragem retrospectiva, recorri a uma propriedade
//do FFBS com o uso do fator de desconto.
//O que observei e' que 
//th_{t-k} = (1-delta) \sum delta^j m_{t-k+j} + delta^k m_T + \sum delta^j e_{T-k+j}
//IntegerVector seq0toT = seq_len(T+1)-1;
//vec powDiscW = exp(log(discW)*as<NumericVector>(seq0toT));
//mat Dmat(T+1, T+1); Dmat.zeros();

//valores aleatorios gerados
mat th = randn(r, T+1);

mm.col(0) = m0;
ZC = ZC0;
th.col(0) = sqrt1mDiscW*ZC*th.col(0);
//Dmat.submat(0, 0, T, 0) = powDiscW.subvec(0, T);

for (int k = 1; k < (T+1); k++){
  //evolucao
  aa = mm.col(k-1);
  ZR = ZC/sqrtDiscW;
  RR = ZR * ZR.t();
  
  //predicao
  FF = kron(I_q, X.row(k-1)).t();
  ff = FF.t()*aa;
  QQ = symmatu(FF.t() * RR * FF + dV);

  //atualizacao
  
  L = ZR.t()*FF; //forma raiz quadrada
  eig_sym(eigval, Eigvec, L * dVInv * L.t());
  ZC = ZR * Eigvec * diagmat(1.0/sqrt(eigval + 1.0));

  AA = ZC*ZC.t()*FF*dVInv; //p.105 West & Harrison (1997)
  ee = Yt.col(k-1) - ff;
  
  mm.col(k) = aa + AA*ee;
  
  th.col(k) = sqrt1mDiscW*ZC*th.col(k);
  mm.col(k-1) -= discW*mm.col(k-1); //mm agora representa mm*
  //Dmat.submat(k, k, T, k) = powDiscW.subvec(0, T-k);
}

// BACKWARD SAMPLING
//simulacao dos estados sem loop
/*
th = (mm+th)*Dmat;
rowvec oneP = ones(X.n_cols);
mat Mu = kron(I_q, oneP) * (th.cols(1, T) % repmat(X.t(), q, 1));
*/
mat Mu(q, T);
th.col(T) += mm.col(T);
for(int k = 0; k < T; k++){
  th.col(T-1-k) += mm.col(T-1-k) + discW*th.col(T-k);
  Mu.col(T-1-k) = kron(I_q, X.row(T-1-k))*th.col(T-k);
}
return Rcpp::List::create(Named("th") = th, Named("Mu") = Mu.t(), 
Named("ZW") = sqrt((1-discW)/discW)*ZC);
}
