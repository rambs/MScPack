# MScPack
# Description: functions to optimize DFM ex-post quadratic loss
# Author: Rafael Barcellos
# Last updated 21st June, 2014
# R 3.1.0

#' @title Build orthogonal matrix
#' @description Build orthogonal matrix from the formula 24 presented in Assmann
#'   et al. (2014).
#' @param parms numeric vector of rotation angles.
#' @return Orthogonal matrix.
#' @references Assmann C, Boysen-Hogrefe J and Pape M (2014). "Bayesian 
#'   Analysis of Dynamic Factor Models: An Ex-Post Approach towards the 
#'   Rotation Problem." \emph{Technical Report 1902}, Kiel Institution 
#'   for the World Economy.

BuildOrthMat <- function(parms){
  
  #if (any(abs(parms) >= pi)){
  #  stop("Argument 'parms' must have elements in -pi < x < pi.")
  #}
  
  n <- length(parms)
  k <- (1 + sqrt(8*n + 1))/2
  
  if (k != floor(k)){
    k <- floor(k)
    n.used <- k*(k-1)/2
    warning(paste("Only the first", n.used, "elements in 'parms' were used."))
  }
  
  it <- 0
  list.x <- list()
  for (i in 1:(k-1)){
    for (j in (i+1):k){
      it <- it + 1
      th <- parms[it]
      x <- diag(1, k)    
      x[i, i] <- cos(th)
      x[i, j] <- sin(th)
      x[j, i] <- -sin(th)
      x[j, j] <- cos(th)      
      list.x[it] <- list(x)
    }
  }
  D.plus <- Reduce("%*%", list.x)
  return(D.plus)
}

# DFM Ex-post loss --------------------------------------------------------

#' @title Ex-post quadratic loss
#' @description Function that calculates the quadratic loss function of the
#'   ex-post approach to solve the rotation problem.
#' @param Lambda DFM loadings matrices piled in a single matrix.
#' @param LambdaStar objective piled loadings matrix.
#' @param PhiBar matrix of VAR parameters joined side by side.
#' @param PhiStar objetive joined matrices of VAR parameters.
#' @return The quadratic loss is returned.

ExPostLoss <- function(Lambda, LambdaStar, PhiBar, PhiStar){
  L1 <- sum((Lambda - LambdaStar)^2)
  L2 <- sum((PhiBar - PhiStar)^2)
  L1 + L2
}


#' @title Ex-post quadratic loss function.
#' @description This function returns the quadratic loss function 
#'   of an orthogonal transformation of the loadings and the VAR parameters
#'   in a DFM (\emph{dynamic factor model}).
#' @param Lambda matrix of loadings of contemporaneous and lagged factors.
#'   See details.
#' @param parms vector of rotation angles to build the orthogonal matrix.
#' @param LambdaStar matrix of loadings to be achieved by a loss minimization
#'   procedure.
#' @param PhiBar matrix concatenation of VAR parameter matrices. See details.
#' @param PhiStar VAR parameter to be achieved by a procedure of loss 
#'   minimization.
#' @details The arguments \code{Lambda} and \code{LambdaStar} must be 
#'   a row binding of all loadings matrices. For example, suppose we have
#'   \eqn{s = 2}. In this case, we have three matrices of loadings which are 
#'   \eqn{\Lambda_0}{\Lambda[0]}, \eqn{\Lambda_1}{\Lambda[1]} and 
#'   \eqn{\Lambda_2}{\Lambda[2]}. So, \code{Lambda} is equal to 
#'   \eqn{(\Lambda_0, \Lambda_1', \Lambda_2')'}{
#'   (\Lambda[0]^T,\Lambda[1]^T,\Lambda[2]^T)^T}.
#' @return Quadratic loss and, as an attribute, the best sign of the orthogonal
#'   matrix.
ExPostLossFunction <- function (Lambda, parms, LambdaStar, PhiBar, PhiStar) {
  
  k <- ncol(Lambda)
  D <- list()
  length(D) <- 2
  D[[1]] <- BuildOrthMat(parms)
  if (nrow(D[[1]]) != k){
    stop("Number of elements in 'parms' misspecified.")
  }

  D.minus <- D[[1]]
  D.minus[k, ] <- -D.minus[k, ]
  D[[2]] <- D.minus
  
  if (nrow(PhiBar) != k) {
    stop("The number of cols in 'Lambda' must be equal 
         to the number of rows in 'PhiBar'.")
  }
  h <- ncol(PhiBar)/k
  Phi <- PhiBar
  dim(Phi) <- c(k, k, h)
  
  Loss <- function (D) {
    Lambda.D <- Lambda %*% D
    L1 <- sum((Lambda.D - LambdaStar)^2)
    
    PhiR <- apply(Phi, 3, function(x) t(D) %*% x %*% D)
    if (nrow(PhiBar) != k) {
      stop("Number of rows in 'PhiBar' must be equal to the number of cols in
         'Lambda' and 'LambdaStar'.")
    }
    dim(PhiR) <- c(k, k*h)
    L2 <- sum((PhiR - PhiStar)^2)
    L <- L1 + L2
    return (L)
  }
  losses <- sapply(D, Loss)
  sign.min <- which.min(losses)
  out <- losses[sign.min]
  attr(out, "sign") <- c(1L, -1L)[sign.min]
  return (out)
}

#' @title Ex-post quadratic loss optimization
#' @description This function uses a heuristic optimization to the
#'   quadratic loss function of the ex-post approach towards the rotation
#'   problem.
#' @param parms vector of rotation angles that build an orthogonal matrix.
#' @return Output of the \code{optim()} function.
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
}

#' @title Rotate matrix of VAR parameters
#' @description Given an orthogonal matrix D, do right and left rotations to 
#'   the matrix of VAR parameters.
#' @param PhiBar matrix of binded matrices of VAR parameters.
#' @param D rotation matrix.
#' @details The rotation is of the form \eqn{D' \Phi[j] D}, where 
#' \eqn{j = 1, \ldots , h}.
#' @return Rotated \code{PhiBar}.
RotatePhi <- function (PhiBar, D) {
  if (nrow(PhiBar) != nrow(D)) {
    stop("'PhiBar' and 'D' must have the same number of rows.")
  }
  k <- nrow(PhiBar)
  h <- ncol(PhiBar)/k
  Phi <- PhiBar
  dim(Phi) <- c(k, k, h)
  Phi <- apply(Phi, 3, function(x) t(D) %*% x %*% D)
  dim(Phi) <- c(k, k*h)
  return(Phi)
}

# Rotated loss ------------------------------------------------------------

#' @title Build loss function of rotation matrix
#' @description This function builds the quadratic loss function of a
#'   rotation matrix using the \code{BuildOrthMat} function parameterisation.
#' @param Lambda loadings matrix.
#' @param LambdaStar objective loadings matrix.
#' @param PhiBar matrix of VAR parameters.
#' @param PhiStar objetive matrix of VAR parameters.
#' @return A function with two arguments: 
#' \itemize{
#'   \item \code{parms} (which contains the parameters that build the orthogonal 
#'     matrix); and
#'   \item \code{reflexion} (element which defines if the orthogonal matrix has
#'     positive or negative determinant).
#' }
BuildExPostRtdLoss <- function(Lambda, LambdaStar, PhiBar, PhiStar){
  function(parms, reflexion){
    if(reflexion != -1 & reflexion != 1){
      stop("'reflexion' must be equal to -1 or 1.")
    }
    if(any(abs(parms) >= pi)){
      return(Inf)
    }
    D <- BuildOrthMat(parms)
    D[nrow(D), ] <- reflexion * D[nrow(D), ]
    Lambda.rtd <- Lambda %*% D
    Phi.rtd <- RotatePhi(PhiBar, D)
    ExPostLoss(Lambda.rtd, LambdaStar, Phi.rtd, PhiStar)
  }  
}

# DFM quadratic loss optimisation -----------------------------------------

#' @title DFM quadratic loss optimisation
#' @description This function enables one to optimise the quadratic loss
#'   of a DFM.
#' @param Lambda loadings matrix. See details.
#' @param LambdaStar loadings matrix objective.
#' @param PhiBar VAR parameters defining the factors evolution.
#' @param PhiStar objective value of PhiBar.
#' @return A list which contains rotation matrix, rotated loadings, 
#'   and rotated PhiBar.
ExPostLossOptim <- function(Lambda, LambdaStar, PhiBar, PhiStar,
                            n.values = 0, seed = as.numeric(Sys.time())){
  # applying wop
  svd.S <- svd(t(Lambda) %*% LambdaStar)
  D.wop <- svd.S$u %*% t(svd.S$v)
  
  Lambda.wop <- Lambda %*% D.wop
  Phi.wop <- RotatePhi(PhiBar, D.wop)
  # building rotation loss
  RtdLoss <- BuildExPostRtdLoss(Lambda.wop, LambdaStar, Phi.wop, PhiStar)
  # optimizing loss for each reflexion
  reflexions <- c(-1, 1)
  opt.rtd <- lapply(reflexions, function(x){
    k <- nrow(D.wop)
    set.seed(seed)
    if(n.values > 0){
      init.values <- array(runif(n.values*k*(k-1)/2, -pi, pi), 
                           c(n.values, k*(k-1)/2))  
    } else {
      init.values <- NULL
    }
    init.values <- rbind(rep(0, k*(k-1)/2), init.values)
    
    opt <- apply(init.values, 1, function(init){
      optim(init, RtdLoss, reflexion = x, method = "L-BFGS-B", 
          lower = -pi+1e-14, upper = pi-1e-14)
    })
    posit.opt <- which.min(sapply(opt, function(z) z$value))
    return(opt[[posit.opt]])
  })
  which.reflexion <- which.min(sapply(opt.rtd, function(x) x$value))
  D <- BuildOrthMat(opt.rtd[[which.reflexion]]$par)
  D[nrow(D), ] <- reflexions[which.reflexion] * D[nrow(D), ]
  # after wop and numerical optimization
  D.opt <- D.wop %*% D
  Lambda.opt <- Lambda %*% D.opt  
  Phi.opt <- RotatePhi(PhiBar, D.opt)
  return(list(Lambda.opt = Lambda.opt, Phi.opt = Phi.opt, D.opt = D.opt))
}

# optmisation algorithm ---------------------------------------------------

#' @title Ex-post DFM loss optimisation
#' @description Apply the optimisation algorithm that solves the rotation
#'   problem for the dynamic factor model.
#' @param LambdaBar array of simulated loadings matrix.
#' @param PhiBar array of simulated VAR parameters.
#' @param tol tolerance level of the optimisation.
#' @param max.iter maximum number of iterations for the algorithm.
#' @return A list containing the rotated parameters as well as the objective 
#'   (best) values for them.
RunDfmExPostOptAlg <- function(LambdaBar, PhiBar, tol = 1e-7, max.iter = 10){
  start.time <- Sys.time()
  
  N <- dim(LambdaBar)[3]
  k <- nrow(PhiBar)
  q <- nrow(LambdaBar)
  
  # initial values
  Lstar <- LambdaBar[,, N]
  dim(Lstar) <- c(q, k, s+1)
  Lstar <- do.call("rbind", lapply(1:(s+1), function(i) Lstar[,,i]))
  Pstar <- PhiBar[,, N]
  
  LambdaD <- array(NA, c(dim(Lstar), N))
  Phi.rtd <- array(NA, dim(PhiBar))
  D <- array(NA, c(ncol(Lstar), ncol(Lstar), N))
  
  if (max.iter <= 0){
    stop("Argument 'max.iter' must be greater than zero.")
  }
  it <-0
  if (tol <= 0){
    stop("Argument 'tol' must be greater than zero.")
  }
  eps <- 2*tol
  
  pb <- txtProgressBar(style = 3)
  
  while (eps > tol){
    if (it >= max.iter){
      stop("Iteration limit reached without convergence.")
    }
    Lstar0 <- Lstar
    Pstar0 <- Pstar
    for (j in 1:N){
      Ls <- LambdaBar[,, j]
      dim(Ls) <- c(q, k, s+1)
      Ls <- do.call("rbind", lapply(1:(s+1), function(i) Ls[,, i]))
      Ps <- PhiBar[,, j]
      out.opt <- ExPostLossOptim(Ls, Lstar, Ps, Pstar)#, n.values = 0)
      LambdaD[,, j]  <- out.opt$Lambda.opt
      D[,, j] <- out.opt$D.opt
      Phi.rtd[,, j] <- out.opt$Phi.opt
      setTxtProgressBar(pb, j/N)
    }
    Lstar <- apply(LambdaD, c(1, 2), mean)
    Pstar <- apply(Phi.rtd, c(1, 2), mean)
    eps <- sum((Lstar-Lstar0)^2) + sum((Pstar-Pstar0)^2)
    it <- it + 1
    message(paste("Iteration", it, "is done."))
  }
  end.time <- Sys.time()
  duration <- end.time - start.time 
  message(paste("Time elapsed:", round(duration, 2), attr(duration, "units")))
  return(list(LambdaD = LambdaD, Phi.rtd = Phi.rtd, D = D, 
              Lstar = Lstar, Pstar = Pstar))
}

