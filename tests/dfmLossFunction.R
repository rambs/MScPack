# MScPack
# Description: function to optimize DFM WOP quadratic loss
# Author: Rafael Barcellos
# Last updated 7 June 2014
# R 3.0.2

library(MScPack)
set.seed(8739)
parms <- runif(3, -pi, pi)
BuildOrthMat(parms)

# example optimization ----------------------------------------------------

q <- 700
k <- 10
h <- 4
LambdaStar <- array(rnorm(q*k), c(q, k))
Lambda <- array(rnorm(q*k), c(q, k))
Phi <- array(rnorm(k*k*h, 0, 1), c(k, k*h))
PhiStar <- array(rnorm(k*k*h, 0, 0), c(k, k*h))

# finding best orthogonal transformation
out.opt <- ExPostLossOptim(Lambda, LambdaStar, Phi, PhiStar)
par(mfrow = c(2, 2), mar = c(2.1, 2.1, 2.1, 0.5))
plot(LambdaStar[, 1], Lambda[, 1], pch = as.character(1:q))
plot(LambdaStar[, 2], Lambda[, 2], pch = as.character(1:q))
plot(LambdaStar[, 1], out.opt$Lambda.opt[, 1], pch = as.character(1:q))
plot(LambdaStar[, 2], out.opt$Lambda.opt[, 2], pch = as.character(1:q))

WopLoss <-  function (Lambda, parms, LambdaStar, PhiBar, PhiStar) {
  
  k <- ncol(Lambda)
  D <- list()
  length(D) <- 2
  D[[1]] <- BuildOrthMat(parms)
  if (nrow(D[[1]]) != k){
    stop("Number of elements in 'parms' misspecified.")
  }
  # reflexion <- diag(c(rep(1, k-1), -1), k)
  D.minus <- D[[1]]
  D.minus[k, ] <- -D.minus[k, ]
  D[[2]] <- D.minus
  
  h <- ncol(PhiBar)/k
  Phi <- PhiBar
  dim(Phi) <- c(k, k, h)
  
  Loss <- function (D) {
    Lambda.D <- Lambda %*% D
    L1 <- sum(colSums((Lambda.D - LambdaStar)^2))

    PhiR <- apply(Phi, 3, function(x) t(D) %*% x %*% D)
    if (nrow(PhiBar) != k) {
      stop("Number of rows in 'PhiBar' must be equal to the number of cols in
         'Lambda' and 'LambdaStar'.")
    }
    dim(PhiR) <- c(k, k*h)
    L2 <- sum(colSums((PhiR - PhiStar)^2))
    L <- L1 + L2
    return (L)
  }
  losses <- sapply(D, Loss)
  sign.min <- which.min(losses)
  out <- losses[sign.min]
  attr(out, "sign") <- c(1L, -1L)[sign.min]
  return (out)
}



parms <- runif(k*(k-1)/2, -pi, pi)



WopDynFactors <- function(Lambda, LambdaStar, PhiBar, PhiStar) {
  
  svd.S <- svd(t(Lambda) %*% LambdaStar)
  D.wop <- svd.S$u %*% t(svd.S$v)
  Lambda.wop <- Lambda %*% D.wop
  PhiBar.wop <- RotatePhi(PhiBar, D.wop)
  
  function(parms){
    ExPostLossFunction(Lambda = Lambda.wop, parms = parms, 
                       LambdaStar = LambdaStar,
                       PhiBar = PhiBar.wop, PhiStar = PhiStar)
  }
  
#   parms.opt <- optim(par = parms.init, fn = ExPostLossFunction, 
#                      Lambda = Lambda.wop, LambdaStar = LambdaStar,
#                      PhiBar = PhiBar.wop, PhiStar = PhiStar, 
#                      lower = -pi+1e-12, upper = pi-1e-12, 
#                      method = "L-BFGS-B")
#   return(list(D.wop = D.wop, optim = parms.opt))
}
WopDfm <- WopDynFactors(Lambda, LambdaStar, PhiBar, PhiStar)
WopDfmEnv <- environment(WopDfm)

wop.opt <- optim(rep(0, k*(k-1)/2), WopDfm, 
                 lower = -pi+1e-12, upper = pi-1e-12, 
                 method = "L-BFGS-B")
D.wop <- get("D.wop", environment(WopDfm))

D.sign <- attr(ExPostLossFunction(Lambda %*% D.wop, wop.opt$par, LambdaStar,
                                  RotatePhi(PhiBar, D.wop), PhiStar),
               "sign")

reflexion <- diag(c(rep(1, nrow(D.wop)-1), D.sign))

D <- D.wop %*% reflexion %*% BuildOrthMat(wop.opt$par)

parms.opt <- apply(parms.init, 1, optim, fn = WopDfm, lower = -pi+1e-12, upper = pi-1e-12, 
                   method = "L-BFGS-B")

N <- 100
set.seed(4812)
parms.init <- matrix(runif(k*(k-1)*N/2, -pi, pi), 
                    N, k*(k-1)/2, byrow = TRUE)
parms.init <- rbind(rep(0, k*(k-1)/2), parms.init)
system.time(
parms.opt <- apply(parms.init, 1, ExPostWopOptim, Lambda = Lambda, 
                   LambdaStar = LambdaStar, PhiBar = PhiBar, 
                   PhiStar = PhiStar)
)
plot(sapply(parms.opt, function(x) x$value))
it.min <- which.min(sapply(parms.opt, function(x) x$value))
parms.min <- parms.opt[[it.min]]$par
loss <- ExPostLossFunction(Lambda.wop, parms.min, LambdaStar, PhiBar.wop, PhiStar)
D <- diag(c(rep(1, k-1), attr(loss, "sign"))) %*% BuildOrthMat(parms.min)

all.equal(D.wop %*% D, D.wop)

sum((Lambda.wop %*% D - LambdaStar)^2)
sum((Lambda.wop - LambdaStar)^2)

Phi <- PhiBar.wop
dim(Phi) <- c(k, k, h)
PhiBar.num <- apply(Phi, 3, function(x) t(D) %*% x %*% D)
dim(PhiBar.num) <- c(k, k*h)

L.num <- sum((Lambda.wop %*% D - LambdaStar)^2) + sum(colSums((PhiBar.num - PhiStar)^2))
L.wop <- sum((Lambda.wop - LambdaStar)^2) + sum((PhiBar.wop - PhiStar)^2)
L.num; L.wop
(L.wop - L.num)/L.wop


install.packages("genalg")
install.packages("rgenoud")
library(genalg)
# optimize two values to match pi and sqrt(50)
evaluate <- function(string=c()) {
  returnVal = NA;
  if (length(string) == 2) {
    returnVal = abs(string[1]-pi) + abs(string[2]-sqrt(50));
  } else {
    stop("Expecting a chromosome of length 2!");
  }
  returnVal
}

monitor <- function(obj) {
  # plot the population
  xlim = c(obj$stringMin[1], obj$stringMax[1]);
  ylim = c(obj$stringMin[2], obj$stringMax[2]);
  plot(obj$population, xlim=xlim, ylim=ylim, xlab="pi", ylab="sqrt(50)");
}

rbga.results = rbga(c(1, 1), c(5, 10), monitorFunc=NULL, 
                    evalFunc=evaluate, verbose=TRUE, mutationChance=0.01)

plot(rbga.results)
plot(rbga.results, type="hist")
plot(rbga.results, type="vars")

WopEvalFunc <- function(parms){
  WopLossFunction(Lambda.wop, parms, LambdaStar, PhiBar.wop, PhiStar)
}

parms.rbga <- rbga(stringMin = rep(-pi+1e-12, k*(k-1)/2), 
                   stringMax = rep(pi-1e-12, k*(k-1)/2),
                   popSize = 20,
                   #suggestions = list(rep(0, k*(k-1)/2)),
                   evalFunc = WopEvalFunc)
str(parms.rbga)
parms.rbga$best
summary(parms.rbga, T)
pdf("R/rgba.pdf")
plot(parms.rbga)
dev.off()
par(mar = c(2.1, 2.1, 2.1, 0.5))
plot(parms.rbga, type="hist")
plot(parms.rbga, type="vars")


library(rgenoud)
nparms <- k*(k-1)/2
parms.genoud <- genoud(WopEvalFunc, nvars = k*(k-1)/2, max = FALSE,
                       starting.values = parms.opt[[it.min]]$par,
                       default.domains = pi-1e-15, boundary.enforcement = 2)
str(parms.genoud)
str(parms.opt[[it.min]])
plot(parms.genoud$par, parms.opt[[it.min]]$par)


