
//[[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;


arma::vec mat2vec(arma::mat y)
{
arma::vec b(y.begin(), y.size(), /*copy_aux_mem*/false, /*strict*/true);
return b;
}

//[[Rcpp::export]]
LogicalVector isNA(NumericVector x) {
  int n = x.size();
  LogicalVector out(n);
  
  for (int i = 0; i < n; ++i) {
    out[i] = NumericVector::is_na(x[i]);
  }
  return out;
}

//[[Rcpp::export]]
IntegerVector cleanNA(NumericVector x) {
  int n = x.size();
  IntegerVector y = seq_len(n);//(x.begin(), x.end());
  IntegerVector z = y[is_na(x)];
  return z;
}

//[[Rcpp::export]]
IntegerVector whichNA(arma::mat X) {
  NumericVector x(X.begin(), X.end());
  int n = x.size();
  IntegerVector y = seq_len(n);//(x.begin(), x.end());
  IntegerVector z = y[is_na(x)];
  if(z.size()==0)
    std::cout << "\nNenhum valor faltante.\n";
  return z;
}

//[[Rcpp::export]]
Rcpp::List changeNA(arma::mat X) {
  NumericVector x(X.begin(), X.end());
  int n = x.size();
  IntegerVector y = seq_len(n)-1;//(x.begin(), x.end());
  IntegerVector z = y[is_na(x)];
  if(z.size()==0){
    std::cout << "Nenhum valor faltante.\n";
    SEXP xNULL = R_NilValue;
    return Rcpp::List::create(xNULL);  
  } else{
    arma::uvec Xna = as<arma::uvec>(z);
    NumericVector bx(Xna.n_rows, 100.0);
    arma::vec b = as<arma::vec>(bx);
    X.elem(Xna) = b;
    return Rcpp::List::create(X);
  }  
}