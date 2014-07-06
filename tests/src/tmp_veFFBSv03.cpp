//[[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
mat veFFBSv03(mat Y, mat eta, vec mu, vec phi, mat U){
//reparar que Y deve ser a matriz de dados T x p que seguem MVE
// eta deve ser p x T

int r = Y.n_cols;  
int T = Y.n_rows;

Environment MScPack("package:MScPack");
Function aproxMN = MScPack["GetLogX2ApproxParms"]; 

mat Ystar = log(pow(Y, 2));

Rcpp::List parmsMN = as<Rcpp::List>(aproxMN(Ystar-eta.t()));

mat M = parmsMN["M"];
mat V2 = parmsMN["S"];
mat E = Ystar - M;
E = E.t();

mat WW = U/(1-phi*phi.t());

mat CC0 = inv(diagmat(1/V2.row(0)) + inv(WW));
vec mm0 = CC0*(diagmat(1/V2.row(0))*E.col(0) + inv(WW)*mu);

// FORWARD FILTERING
// FF é matriz identidade
mat GG = diagmat(phi);
mat dV = diagmat(V2.row(0));
// incorporando a media de w_t ao FFBS de acordo com W&H 1997 pp. 583.
mat Ew = (eye(r, r) - diagmat(phi))*mu;

mat mm(r, T);
cube CC(r, r, T);
mat aa(r, T);
cube RR(r, r, T);
vec ff(r);
mat QQ(r, r);
mat AA;
vec ee;

mm.col(0) = mm0;
CC.slice(0) = CC0;

int k;
for (k=1; k<T; k++){
aa.col(k) = GG*mm.col(k-1) + Ew; //W&H pp.583
RR.slice(k) = symmatu(GG*CC.slice(k-1)*GG.t() + U);

ff = aa.col(k);
dV = diagmat(V2.row(k));
QQ = symmatu(RR.slice(k) + dV);

AA = RR.slice(k)*inv(QQ);
ee = E.col(k) - ff; // filtro inicia em y_2, pois y_1 ja condicionado em m0, C0

mm.col(k) = aa.col(k) + AA*ee;
CC.slice(k) = symmatu(inv(inv(RR.slice(k)) + inv(dV)));
}

// BACKWARD SAMPLING

mat th = randn(r, T);
mat U1;
vec s1; 
mat V1;
svd(U1, s1, V1, CC.slice(T-1), "standard");
mat L = U1*diagmat(sqrt(s1)); //chol(CC.slice(T));

th.col(T-1) = mm.col(T-1) + L*th.col(T-1);

mat invRR;
vec hh;
mat HH;

for (int k=1; k<T; k++){
invRR = inv(RR.slice(T-k));
hh = mm.col(T-1-k) + CC.slice(T-1-k)* GG.t() * invRR * (th.col(T-k) - aa.col(T-k));
HH = symmatu(CC.slice(T-1-k) - CC.slice(T-1-k) * GG.t() * invRR * GG * CC.slice(T-1-k));
svd(U1, s1, V1, HH, "standard");
L = U1*diagmat(sqrt(s1));
th.col(T-1-k) = hh + L*th.col(T-1-k);
}

return th;

}