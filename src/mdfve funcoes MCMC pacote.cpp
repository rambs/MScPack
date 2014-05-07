/**
 * Pacote mdfve
 * Descricao: Funcoes para o algoritmo MCMC
 * Autor: Rafael Barcellos
 * Data: 27/04/2014
 * R 3.0.2
 */

// funcao FFBS para regressao com V constante e fator de desconto
using namespace arma;

Rcpp::List drmFFBSfixVdiscW(mat Y, mat X, mat dV, double discW, vec m0, mat ZC0){
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
  AA = RR * FF * inv(QQ);
  ee = Yt.col(k-1) - ff;
  
  mm.col(k) = aa + AA*ee;
  L = ZR.t()*FF; //forma raiz quadrada
  eig_sym(eigval, Eigvec, L * dVInv * L.t());
  ZC.slice(k) = ZR * Eigvec * diagmat(1.0/sqrt(eigval + 1.0));
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
}//

//funcao para simulacao dos parametros de regressao estaticos
mat parmsFixRegFixdV(mat Y, mat X, mat B1, mat ZB1, vec psi, mat B0Invb0)
{ 
  //B0Invb0 e' parcela da priori para media
  //B1 e' variancia a posteriori que nao altera de uma iteracao para outra do MCMC
  //ZB1 forma quadrada de B1
      
  //mat B1 = inv(XtX + B0Inv);
  //  mat U, V;
  //vec s;
  //svd(U, s, V, B1, "standard")
  //mat ZB1 = U*diagmat(sqrt(s));
  mat b1 = B1*(X.t()*Y + B0Invb0);
  mat beta = b1 + ZB1*randn(X.n_cols, Y.n_cols)*diagmat(sqrt(psi));
  return beta;
  /*
  int r = Y.n_cols*X.n_cols;
  mat Yt = Y.t();//matriz T x q
  mat Xt = X.t();//matriz T x p
  mat dVInv = inv(dV);
  mat B0Inv = inv(B0);
  //sendo dV constante
  mat XtX = X.t()*X;
  mat XtXInv = inv(XtX);
  mat B1Inv = B0Inv + kron(dVInv, XtXInv);
  
  vec b1 = kron(dVInv, Xt.col(0))*Yt.col(0) + B0Inv*b0;
  mat B1Inv = B0Inv + kron(dVInv, Xt.col(0)*X.row(0));
  for(int k = 1; k < T; k++){
    b1 = b1 + kron(dVInv, Xt.col(k))*Yt.col(k);
    B1Inv = B1Inv + kron(dVInv, Xt.col(k)*X.row(k));
  }
  
  mat B1 = inv(B1Inv);
  vec b1 = B1*b1;
  mat U, V;
  vec s;
  svd(U, s, V, B1, "standard");
  vec beta = b1 + U*sqrt(s)*randn(r, 1);
  return beta;
  */
}

//funcao para simular Lambda Matriz-Variada
mat LambdaSimMV(mat Y, mat Factors, vec psi, mat L0H0Inv, mat H0Inv){
  int q = Y.n_cols;
  int k = Factors.n_cols;
  
  mat H1 = inv(Factors.t()*Factors + H0Inv);
  mat L1 = (Y.t()*Factors + L0H0Inv)*H1;
  mat U, V;
  vec s;
  svd(U, s, V, H1, "standard");
  mat Lambda = L1 + diagmat(sqrt(psi))*randn(q, k)*diagmat(sqrt(s))*V;
  return Lambda;
}


//funcao para simular fatores
mat FactorSim(mat Y, mat Lambda, vec psi){
  int k = Lambda.n_cols;
  int T = Y.n_rows;
  mat G1 = inv(Lambda.t()*diagmat(1/psi)*Lambda + eye(k, k));
  mat F1 = (Y*diagmat(1/psi)*Lambda)*G1;
  
  mat UG, VG, svdG;
  vec sG;
  svd(UG, sG, VG, G1, "standard");
  svdG = diagmat(sqrt(sG))*VG;
  mat Factors = F1 + randn(T, k)*svdG;
  return Factors;
}

// simulacao da variancias idiossincraticas
vec psiSim(mat Y, mat Beta, mat b0, mat B0Inv, mat Lambda, mat L0, mat H0Inv, 
            double n0, vec s0){
  Environment stats("package:stats"); //chamando a função rgamma do "stats"
  Function rgammaR = stats["rgamma"];
  
  int T = Y.n_rows;
  int q = Y.n_cols;
  int k = Lambda.n_cols;
  int p = Beta.n_rows;
  vec Spsi = sum(pow(Y, 2)).t() + trace((Beta-b0).t()*B0Inv*(Beta-b0)) + 
                trace((Lambda-L0)*H0Inv*(Lambda-L0).t());
  double n1 = n0 + T + k + p;
  vec n1s1 = (Spsi + n0*s0);
  vec psi = 1/as<arma::vec>(rgammaR(q, n1/2, n1s1/2));
  return psi;
}

// funcao que insere barra de progresso
void progress_bar(double x, double N)
{
    // how wide you want the progress meter to be
    int totaldotz=40;
    double fraction = x / N;
    // part of the progressmeter that's already "full"
    int dotz = round(fraction * totaldotz);

    // create the "meter"
    int ii=0;
    
    printf("\r%3.0f%% [",fraction*100);
    // part  that's full already
    for ( ; ii < dotz;ii++) {
        printf("=");
    }
    // remaining part (spaces)
    for ( ; ii < totaldotz;ii++) {
        printf(" ");
    }
    // and back to line begin - do not forget the fflush to avoid output buffering problems!
    printf("]\r");
    fflush(stdout);
}
