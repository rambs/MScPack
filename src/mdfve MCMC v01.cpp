/***
 * Descricao: Gibbs do modelo dinamico fatorial WOP W via desconto 
 *            FFBS na forma raiz quadrada e outras otimizacoes do codigo
 * Autor: Rafael Barcellos
 * Data: 29/04/2014
 * R 3.0.1
 * 
 */

#include <RcppArmadillo.h>
using namespace arma;

// leitura dos argumentos da funcao
int N = as<int>(n);
int brn = as<int>(burn);
int thn = as<int>(thin);

mat Y = as<arma::mat>(y);
mat xFix = as<arma::mat>(xFixReg); //exogenas reg estatica
mat xDyn = as<arma::mat>(xDynReg); //exogenas reg dinamica

mat LL0 = as<arma::mat>(L0);//hiperparametros de Lambda
mat HH0 = as<arma::mat>(H0);

vec mm0 = as<arma::vec>(m0);//hiperparametros de Theta
mat CC0 = as<arma::mat>(C0);

mat bb0 = as<arma::mat>(b0);//hiperparametros da reg estatica
mat BB0 = as<arma::mat>(B0);

int k = as<int>(nFactors);

double nn0 = as<double>(n0);
vec ss0 = as<arma::vec>(s0sq);

double delta = as<double>(discW);

int q = Y.n_cols; // matriz de dados deve entrar tempo x variavel
int T = Y.n_rows;
int pFix = xFix.n_cols;
int pDyn = xDyn.n_cols;
int rFix = pFix*q;
int rDyn = pDyn*q;

// Algoritmo de Gibbs

mat ths(rDyn, T+1);
mat Ls = LL0;
mat Fs(T, k);
Fs.zeros();
vec ps = ss0;
mat Bs = bb0;
mat ZWs(rDyn, rDyn);

// Amostras finais
cube th(rDyn, T+1, N);
cube Lambda(q, k, N);
cube FF(T, k, N);
mat psi(N, q);
cube Beta(pFix, q, N);
cube ZW(rDyn, rDyn, N);
cube Mu(T, rDyn, N);

//parametros iniciais que aceleram MCMC
mat U, V, ZC0;
vec s;
svd(U, s, V, CC0, "standard");
ZC0 = U * diagmat(sqrt(s));

mat B0Inv = inv(BB0);
mat B0Invb0 = inv(BB0)*bb0;
mat B1 = inv(xFix.t()*xFix + B0Inv);
svd(U, s, V, B1, "standard");
mat ZB1 = U*diagmat(sqrt(s));

mat H0Inv = inv(HH0);
mat L0H0Inv = LL0*H0Inv;

// iterador da barra de progresso
int it = 0;
int itTot = brn + thn*N;

mat E; //matriz que recebe valores corrigidos a media condicional
Rcpp::List thList;
mat MuFix = xFix*Bs;
mat MuDyn(T, q);

int i, j;
//burn-in
if (brn>0){
  for (j = 0; j < brn; j++){
  // simulacao dos estados
  E = Y - MuFix - Fs*Ls.t();
  thList = drmFFBSfixVdiscW(E, xDyn, diagmat(ps), delta, mm0, ZC0);
  ths = as<arma::mat>(thList["th"]);
  MuDyn = as<arma::mat>(thList["Mu"]);
  
  // simulacao da reg estatica
  E = Y - MuDyn - Fs*Ls.t();
  Bs = parmsFixRegFixdV(E, xFix, B1, ZB1, ps, B0Invb0);
  MuFix = xFix*Bs;
  
  // simulacao dos fatores
  E = Y - MuFix - MuDyn;
  Fs = FactorSim(E, Ls, ps);
  
  // simulacao da matriz de cargas
  Ls = LambdaSimMV(E, Fs, ps, L0H0Inv, H0Inv);
  
  // simulacao das variancias idiossincraticas
  E = E - Fs*Ls.t();
  ps = psiSim(E, Bs, bb0, B0Inv, Ls, LL0, H0Inv, nn0, ss0);
  
  it++;
  progress_bar(it, itTot);
  }
} //fecha if (brn>0)


// Gibbs
if (thn<=0) thn = 1;
for (j=0; j<N; j++){
  for (i=0; i<thn; i++){
  // simulacao dos estados
  E = Y - MuFix - Fs*Ls.t();
  thList = drmFFBSfixVdiscW(E, xDyn, diagmat(ps), delta, mm0, ZC0);
  ths = as<arma::mat>(thList["th"]);
  MuDyn = as<arma::mat>(thList["Mu"]);
  
  // simulacao da reg estatica
  E = Y - MuDyn - Fs*Ls.t();
  Bs = parmsFixRegFixdV(E, xFix, B1, ZB1, ps, B0Invb0);
  MuFix = xFix*Bs;
  
  // simulacao dos fatores
  E = Y - MuFix - MuDyn;
  Fs = FactorSim(E, Ls, ps);
  
  // simulacao da matriz de cargas
  Ls = LambdaSimMV(E, Fs, ps, L0H0Inv, H0Inv);
  
  // simulacao das variancias idiossincraticas
  E = E - Fs*Ls.t();
  ps = psiSim(E, Bs, bb0, B0Inv, Ls, LL0, H0Inv, nn0, ss0);
  
  it++;
  progress_bar(it, itTot);
  }
  th.slice(j) = ths;
  Beta.slice(j) = Bs;
  Lambda.slice(j) = Ls;
  FF.slice(j) = Fs;
  psi.row(j) = ps.t();
  ZW.slice(j) = as<arma::mat>(thList["ZW"]);
  Mu.slice(j) = MuDyn;
}
printf("\n Algoritmo concluido. \n");

Rcpp::List mod = Rcpp::List::create(Named("n") = N, Named("burn") = brn, Named("thin") = thn,
Named("y") = Y, Named("xFix") = xFix, Named("xDyn") = xDyn, Named("nFactors") = k, 
Named("m0") = mm0, Named("C0") = CC0, Named("L0") = LL0, Named("G0") = HH0, 
Named("b0") = bb0, Named("B0") = BB0, Named("n0") = nn0, Named("s0sq") = ss0, 
Named("discW") = delta);
return Rcpp::List::create(Named("mod") = mod, Named("th") = th, Named("Beta") = Beta, 
Named("Lambda") = Lambda, 
Named("Factors") = FF, Named("psi") = psi, Named("ZW") = ZW, Named("MuDyn") = Mu);
//
