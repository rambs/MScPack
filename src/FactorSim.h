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