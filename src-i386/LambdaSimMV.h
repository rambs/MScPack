//funcao para simular Lambda Matriz-Variada
mat LambdaSimMV(mat Y, mat Factors, vec psi, mat L0H0Inv, mat H0Inv){
  int q = Y.n_cols;
  int k = Factors.n_cols;
  
  mat H1 = inv(Factors.t()*Factors + H0Inv);
  mat L1 = (Y.t()*Factors + L0H0Inv)*H1;
  //mat U, V;
  //vec s;
  //svd(U, s, V, H1, "standard");
  mat Eigvec;
  vec eigval;
  eig_sym(eigval, Eigvec, H1);
  mat ZHt = diagmat(sqrt(eigval))*Eigvec.t();
  mat Lambda = L1 + diagmat(sqrt(psi))*randn(q, k)*ZHt;//diagmat(sqrt(s))*V;
  return Lambda;
}