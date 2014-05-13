
//[[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <Rcpp/Benchmark/Timer.h>

using namespace Rcpp;
using namespace arma;

#include "drmFFBSdiscW.h"
#include "parmsFixRegSim.h"
#include "LambdaSimMV.h"
#include "FactorSim.h"
#include "psiSim.h"
#include "progress_bar.h"

//' Algoritmo de Gibbs para F-DLM
//' 
//' Estima as quantidades do modelo dinamico fatorial de acordo com os resultados presentes na dissertacao.
//' @param N tamanho final da amostra.
//' @param brn periodo de aquecimento do algoritmo.
//' @param model modelo F-DLM com valores da matriz de dados, de exogenas e hiperparametros da priori.
//' @param initVal valores iniciais para Beta, Lambda e psi.
//' @return Lista com modelo e valores simulados do Gibbs.
// [[Rcpp::export]]
Rcpp::List fdlmGibbs(int N, int brn, int thn, Rcpp::List model, Rcpp::List initVal, 
bool progressBar = true, bool onlyValues = false)
{
  Rcpp::Timer timer;
  // leitura dos argumentos da funcao
mat Y = as<arma::mat>(model["y"]);
mat xFix = as<arma::mat>(model["xFixReg"]); //exogenas reg estatica
mat xDyn = as<arma::mat>(model["xDynReg"]); //exogenas reg dinamica
int k = as<int>(model["nFactors"]);

mat LL0 = as<arma::mat>(model["L0"]);//hiperparametros de Lambda
mat HH0 = as<arma::mat>(model["H0"]);

vec mm0 = as<arma::vec>(model["m0"]);//hiperparametros de Theta
mat CC0 = as<arma::mat>(model["C0"]);

mat bb0 = as<arma::mat>(model["b0"]);//hiperparametros da reg estatica
mat BB0 = as<arma::mat>(model["B0"]);

double nn0 = as<double>(model["n0"]);//hiperparametros da variancia idiossincratica
vec ss0 = as<arma::vec>(model["s0sq"]);

double delta = as<double>(model["discW"]); //fator de desconto

int q = Y.n_cols; // matriz de dados deve entrar tempo x variavel
int T = Y.n_rows;
int pFix = xFix.n_cols;
int pDyn = xDyn.n_cols;
int rDyn = pDyn*q;

// Algoritmo de Gibbs

mat ths(rDyn, T+1);
mat MDs(T, q); //media dinamica
mat Ls(q, k);
mat Fs(T, k); Fs.zeros();
vec ps(q);
mat Bs(pFix, q);
mat MFs(T, q); //media fixa
mat ZWs(rDyn, rDyn);

// Pontos iniciais da cadeia
Ls = as<arma::mat>(initVal["Lambda"]);
ps = as<arma::vec>(initVal["psi"]);
Bs = as<arma::mat>(initVal["Beta"]);
MFs = xFix*Bs;
MDs.zeros();

// Amostras finais
cube theta(rDyn, T+1, N);
cube Lambda(q, k, N);
cube Factors(T, k, N);
mat psi(N, q);
cube Beta(pFix, q, N);
cube ZW(rDyn, rDyn, N);
cube MuDyn(T, q, N);
cube MuFix(T, q, N);
cube Mu(T, q, N);

//parametros iniciais que aceleram MCMC
mat U, V, ZC0;
vec s;
svd(U, s, V, CC0, "standard");
ZC0 = U * diagmat(sqrt(s));

mat B0Inv = inv(BB0);
mat B0Invb0 = inv(BB0)*bb0;
mat B1 = symmatu(inv(xFix.t()*xFix + B0Inv));
mat Eigvec;
vec eigval;
eig_sym(eigval, Eigvec, B1);
mat ZB1 = Eigvec*diagmat(sqrt(eigval));
//svd(U, s, V, B1, "standard");
//mat ZB1 = U*diagmat(sqrt(s));

mat H0Inv = inv(HH0);
mat L0H0Inv = LL0*H0Inv;

// iterador da barra de progresso
int it = 0;
int itTot = brn + thn*N;

mat E; //matriz que recebe valores corrigidos a media condicional
Rcpp::List thList;


int i, j;
//burn-in
if (brn>0){
  for (j = 0; j < brn; j++){
  // simulacao dos estados
  E = Y - MFs - Fs*Ls.t();
  thList = drmFFBSdiscW(E, xDyn, diagmat(ps), delta, mm0, ZC0);
  ths = as<arma::mat>(thList["th"]);
  MDs = as<arma::mat>(thList["Mu"]);
  
  // simulacao da reg estatica
  E = Y - MDs - Fs*Ls.t();
  Bs = parmsFixRegSim(E, xFix, B1, ZB1, ps, B0Invb0);
  MFs = xFix*Bs;
  
  // simulacao dos fatores
  E = Y - MFs - MDs;
  Fs = FactorSim(E, Ls, ps);
  
  // simulacao da matriz de cargas
  Ls = LambdaSimMV(E, Fs, ps, L0H0Inv, H0Inv);
  
  // simulacao das variancias idiossincraticas
  E = E - Fs*Ls.t();
  ps = psiSim(E, Bs, bb0, B0Inv, Ls, LL0, H0Inv, nn0, ss0);
  
  if(progressBar){
    it++;
    progress_bar(it, itTot);  
  }
  
  }
} //fecha if (brn>0)


// Gibbs
if (thn<=0) thn = 1;
for (j=0; j<N; j++){
  for (i=0; i<thn; i++){
  // simulacao dos estados
  E = Y - MFs - Fs*Ls.t();
  thList = drmFFBSdiscW(E, xDyn, diagmat(ps), delta, mm0, ZC0);
  ths = as<arma::mat>(thList["th"]);
  MDs = as<arma::mat>(thList["Mu"]);
  
  // simulacao da reg estatica
  E = Y - MDs - Fs*Ls.t();
  Bs = parmsFixRegSim(E, xFix, B1, ZB1, ps, B0Invb0);
  MFs = xFix*Bs;
  
  // simulacao dos fatores
  E = Y - MFs - MDs;
  Fs = FactorSim(E, Ls, ps);
  
  // simulacao da matriz de cargas
  Ls = LambdaSimMV(E, Fs, ps, L0H0Inv, H0Inv);
  
  // simulacao das variancias idiossincraticas
  E = E - Fs*Ls.t();
  ps = psiSim(E, Bs, bb0, B0Inv, Ls, LL0, H0Inv, nn0, ss0);
  
  if(progressBar){
    it++;
    progress_bar(it, itTot);
  }
  
  }
  theta.slice(j) = ths;
  Beta.slice(j) = Bs;
  Lambda.slice(j) = Ls;
  Factors.slice(j) = Fs;
  psi.row(j) = ps.t();
  ZW.slice(j) = as<arma::mat>(thList["ZW"]);
  MuDyn.slice(j) = MDs;
  MuFix.slice(j) = MFs;
  Mu.slice(j) = MFs + MDs;
}

if(progressBar){
  printf("\n Algoritmo concluido. \n");
}

timer.step("Time elapsed (in sec)");
Rcpp::NumericVector duration(timer);
duration[0] = duration[0]/1000000000;

Rcpp::List gibbs = Rcpp::List::create(Named("duration") = duration, Named("N") = N, Named("burn") = brn, 
Named("thin") = thn);

Rcpp::List values = Rcpp::List::create(Named("theta") = theta, Named("Beta") = Beta,
Named("Lambda") = Lambda, Named("Factors") = Factors, Named("psi") = psi, Named("ZW") = ZW,
Named("MuFix") = MuFix, Named("MuDyn") = MuDyn, Named("Mu") = Mu);

if(onlyValues){
  return values;
} else{
  return Rcpp::List::create(Named("model") = model, Named("gibbs") = gibbs, 
  Named("values") = values);  
}

}
