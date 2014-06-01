
//[[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <Rcpp/Benchmark/Timer.h>

using namespace Rcpp;
using namespace arma;

#include "SampleIntDynRegDiscW.h"
#include "parmsFixRegSim.h"
#include "LambdaSimMV.h"
#include "FactorSim.h"
//#include "psiSim.h"
#include "progress_bar.h"

//' Amostrador de Gibbs para F-DRM com intercepto
//' 
//' Roda o amostrador de Gibbs para F-DRM onde o intercepto
//'  e os parametros de estados sao amostrados conjuntamente.
//' @param N tamanho final da amostra.
//' @param brn periodo de aquecimento do algoritmo.
//' @param thn espacamento entre as extracoes.
//' @param model modelo F-DLM com valores da matriz de dados,
//'  de exogenas e hiperparametros da priori.
//' @param initVal valores iniciais para Beta, Lambda e psi.
//' @return Lista com modelo e valores simulados do Gibbs.
// [[Rcpp::export]]
Rcpp::List RunFdlmIntGibbs(int N, int brn, int thn, Rcpp::List model, Rcpp::List initVal, 
bool progressBar = true, bool onlyValues = false)
{
  //contador do tempo
  Rcpp::Timer timer;
  
  //chamando a função rgamma do "stats"
  Rcpp::Environment stats("package:stats"); 
  Rcpp::Function rgammaR2 = stats["rgamma"];

  // leitura dos argumentos da funcao
  mat Y = as<arma::mat>(model["y"]);
  
  mat xDyn = as<arma::mat>(model["xDynReg"]); //exogenas reg dinamica
  
  int k = as<int>(model["nFactors"]);
  
  mat LL0 = as<arma::mat>(model["L0"]);//hiperparametros de Lambda
  mat HH0 = as<arma::mat>(model["H0"]);
  
  vec mm0 = as<arma::vec>(model["m0"]);//hiperparametros de Theta
  mat CC0 = as<arma::mat>(model["C0"]);
  
  vec aa0 = as<arma::vec>(model["a0"]);//hiperparametros do intercepto
  mat AA0 = as<arma::mat>(model["A0"]);
  
  double nn0 = as<double>(model["n0"]);//hiperparametros da variancia idiossincratica
  vec ss0 = as<arma::vec>(model["s0sq"]);
  
  double delta = as<double>(model["discW"]); //fator de desconto
  
  int q = Y.n_cols; // matriz de dados deve entrar tempo x variavel
  int T = Y.n_rows;
  int pDyn = xDyn.n_cols;
  int rDyn = pDyn*q;
  
  std::cout << "\n Comecou\n";
  //tratamento de missings
  NumericVector y(Y.begin(), Y.end());//vetorizacao da matriz de dados
  int ysize = y.size();
  IntegerVector yseq = seq_len(ysize)-1;
  IntegerVector yNA = yseq[is_na(y)];
  uvec uyNA = as<arma::uvec>(yNA);
    
  // Algoritmo de Gibbs
  mat Ys(T, q);//objeto para lidar com missings
  vec ys(yNA.size());
  
  vec ints(q);
  mat ths(rDyn, T+1);
  mat MDs(T, q); //media dinamica
  mat Ls(q, k);
  mat Fs(T, k); 
  vec ps(q);
  mat psMat(T, q); //matriz de var. idio. para simulacao dos dados faltantes
  mat ZWs(rDyn, rDyn);
  mat Ms;

  // Pontos iniciais da cadeia
  Ys = Y;
  if(yNA.size()>0){
    NumericVector ysInit(yNA.size(), 0.0);
    Ys.elem(uyNA) = as<arma::vec>(ysInit);
    }
    
  Ls = as<arma::mat>(initVal["Lambda"]);
  Fs.zeros();
  ps = as<arma::vec>(initVal["psi"]);
  MDs.zeros();

  // Amostras finais
  mat alpha(N, q); //intercepto
  cube theta(rDyn, T+1, N);
  cube Lambda(q, k, N);
  cube Factors(T, k, N);
  mat psi(N, q);
  cube rootW(rDyn, rDyn, N);
  cube MuDyn(T, q, N);
  mat Y_NA;
  if(yNA.size()>0){
    Y_NA.resize(N, yNA.size());
    } 

  //parametros iniciais que aceleram MCMC  
  mat H0Inv = inv(HH0);
  mat L0H0Inv = LL0*H0Inv;
  
  // iterador da barra de progresso
  int it = 0;
  int itTot = brn + thn*N;
  
  mat E; //matriz que recebe valores corrigidos a media condicional completa
  Rcpp::List thList;
  
  std::cout << "Rodou ate antes do burn-in";
  int i, j;
  //burn-in
  if (brn > 0){
    for (j = 0; j < brn; j++){
    // simulacao dos estados
    E = Ys - Fs*Ls.t();
    thList = SampleIntDynRegDiscW(E, xDyn, diagmat(ps), delta, mm0, CC0, aa0, AA0);
    ints = as<arma::vec>(thList["alpha"]);
    ths = as<arma::mat>(thList["th"]);
    MDs = as<arma::mat>(thList["Mu"]);
  
    // simulacao dos fatores
    E = Ys - MDs;
    Fs = FactorSim(E, Ls, ps);
  
    // simulacao da matriz de cargas
    Ls = LambdaSimMV(E, Fs, ps, L0H0Inv, H0Inv);
  
    // simulacao das variancias idiossincraticas
    E = E - Fs*Ls.t();
    vec Spsi = sum(E%E).t() + trace((Ls-LL0)*H0Inv*(Ls-LL0).t());
    double n1 = nn0 + T + k;
    vec n1s1 = (Spsi + nn0*ss0);
    ps = 1/as<arma::vec>(rgammaR2(q, n1/2, n1s1/2));
  
    //simulacao dos valores faltantes
    if(yNA.size() > 0){
      Ms = MDs + Fs*Ls.t();
      psMat = repmat(ps.t(), T, 1);//variancia idiossincratica constante
      ys = Ms.elem(uyNA) + diagmat(sqrt(psMat.elem(uyNA)))*randn(yNA.size(), 1);
      Ys.elem(uyNA) = ys;    
    }
  
    if(progressBar){
      it++;
      progress_bar(it, itTot);  
    }
  
  }
} //fecha if (brn > 0)


// Gibbs
if (thn <= 0) thn = 1;
for (j = 0; j < N; j++){
  for (i = 0; i < thn; i++){
    // simulacao dos estados
    E = Ys - Fs*Ls.t();
    thList = SampleIntDynRegDiscW(E, xDyn, diagmat(ps), delta, mm0, CC0, aa0, AA0);
    ints = as<arma::vec>(thList["alpha"]);
    ths = as<arma::mat>(thList["th"]);
    MDs = as<arma::mat>(thList["Mu"]);
  
    // simulacao dos fatores
    E = Ys - MDs;
    Fs = FactorSim(E, Ls, ps);
  
    // simulacao da matriz de cargas
    Ls = LambdaSimMV(E, Fs, ps, L0H0Inv, H0Inv);
  
    // simulacao das variancias idiossincraticas
    E = E - Fs*Ls.t();
    vec Spsi = sum(E%E).t() + trace((Ls-LL0)*H0Inv*(Ls-LL0).t());
    double n1 = nn0 + T + k;
    vec n1s1 = (Spsi + nn0*ss0);
    ps = 1/as<arma::vec>(rgammaR2(q, n1/2, n1s1/2));
  
    //simulacao dos valores faltantes
    if(yNA.size() > 0){
      Ms = MDs + Fs*Ls.t();
      psMat = repmat(ps.t(), T, 1);//variancia idiossincratica constante
      ys = Ms.elem(uyNA) + diagmat(sqrt(psMat.elem(uyNA)))*randn(yNA.size(), 1);
      Ys.elem(uyNA) = ys;    
    }
  
    if(progressBar){
      it++;
      progress_bar(it, itTot);  
    }
  }
  alpha.row(j) = ints.t();
  theta.slice(j) = ths;
  Lambda.slice(j) = Ls;
  Factors.slice(j) = Fs;
  psi.row(j) = ps.t();
  rootW.slice(j) = as<arma::mat>(thList["rootW"]);
  MuDyn.slice(j) = MDs;
  if(yNA.size()>0){
    Y_NA.row(j) = ys.t();  
  }
}

timer.step("Time elapsed (in sec)");
Rcpp::NumericVector duration(timer);
duration[0] = duration[0]/1000000000;

if(progressBar){
  printf("\n Algoritmo concluido. \n");
}

Rcpp::List gibbs = Rcpp::List::create(Named("duration") = duration, Named("N") = N, Named("burn") = brn, 
Named("thin") = thn);

Rcpp::List values = Rcpp::List::create(Named("Y_NA") = Y_NA, Named("yNA") = yNA, 
  Named("alpha") = alpha, Named("theta") = theta, 
  Named("Lambda") = Lambda, Named("Factors") = Factors, 
  Named("psi") = psi, Named("rootW") = rootW, 
  Named("MuDyn") = MuDyn);

if(onlyValues){
  return values;
} else{
  return Rcpp::List::create(Named("model") = model, Named("gibbs") = gibbs, 
  Named("values") = values);  
}

}
