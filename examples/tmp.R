x = rnorm(5)
tmp = chol(tcrossprod(x), pivot = TRUE)
pivot = order(attr(tmp, "pivot"))
crossprod(tmp)
tcrossprod(x)[pivot, pivot]
tmp = bayesm::rwishart(1,diag(8))
crossprod(tmp$C)
tmp$W

library(MScPack)
cores <- parallel::detectCores()
set.seed(123)
sigma <- bayesm::rwishart(10000,diag(9))$IW
set.seed(1928)
Lambda = array(runif(q*k, 0, Lambda.lim), c(q, k))
Lambda[upper.tri(Lambda)] = 0
diag(Lambda) = c(0.99, 0.95, 0.9)*Lambda.lim
psi = c(0.02, 0.19, 0.36, 0.02, 0.02, 0.19, 0.19, 0.36, 0.36)*Lambda.lim/10

sigma <- tcrossprod(Lambda) + diag(psi)
means <- rnorm(9)
X     <- mvtnorm::rmvnorm(500, means, sigma)
system.time(LL <- sum(dmvnrm_arma_mc(X,means,sigma,T,cores)))
LL
