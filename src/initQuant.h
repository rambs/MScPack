//' Quantidades iniciais para Gibbs F-DLM
//' 
//' Calculo de quantidades inicias para evitar repeticao de operacoes identicas a cada iteracao do algoritmo de Gibbs..
//' @param model modelo F-DLM com valores da matriz de dados, de exogenas e hiperparametros da priori.
//' @return Lista com modelo e valores iniciais para a cadeia do Gibbs.
// [[Rcpp::export]]
Rcpp::List fdlmInitQuant(Rcpp::List model, bool modelOut = true)
{
  // leitura dos argumentos da funcao
mat Y = as<arma::mat>(model["y"]);
mat xFix = as<arma::mat>(model["xFixReg"]); //exogenas reg estatica
mat xDyn = as<arma::mat>(model["xDynReg"]); //exogenas reg dinamica

mat LL0 = as<arma::mat>(model["L0"]);//hiperparametros de Lambda
mat HH0 = as<arma::mat>(model["H0"]);

vec mm0 = as<arma::vec>(model["m0"]);//hiperparametros de Theta
mat CC0 = as<arma::mat>(model["C0"]);

mat bb0 = as<arma::mat>(model["b0"]);//hiperparametros da reg estatica
mat BB0 = as<arma::mat>(model["B0"]);

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

mat H0Inv = inv(HH0);
mat L0H0Inv = LL0*H0Inv;

//Lista com quantidade iniciais
Rcpp::List initQuant = Rcpp::List::create(Named("ZC0") = ZC0, Named("B1") = B1,
Named("ZB1") = ZB1, Named("H0Inv") = H0Inv, Named("L0H0Inv") = L0H0Inv);
if(modelOut){
  return Rcpp::List::create(Named("model") = model, Named("initQuant") = initQuant);
} else{
  return initQuant;  
}

}