# MScPack
# Description: testing the mixture of normals approximation
# Author: Rafael Barcellos
# Last updated 3rd July, 2014
# R 3.1.0

TT <-  2.5e3
k <- 6
y <- array(log(rnorm(TT*k)^2), c(TT, k))

parms.mn <- GetLogX2ApproxParms(y)

y.approx <- array(rnorm(TT*k, parms.mn$M, sqrt(parms.mn$S)), c(TT, k))

par(mfrow = c(2, 3))
for (i in 1:k){
  hist(y[, i], main = "", breaks = 50, col = rgb(0, 0, 0.5, 0.3),
       xlim = range(y, y.approx), freq = FALSE, ylim = c(0, 0.3))
  par(new = TRUE)
  hist(y.approx[, i], main = "", breaks = 50, col = rgb(0.5, 0, 0, 0.3),
       xlim = range(y, y.approx), freq = FALSE, ylim = c(0, 0.3), axes = F,
       xlab = "", ylab = "")
}
