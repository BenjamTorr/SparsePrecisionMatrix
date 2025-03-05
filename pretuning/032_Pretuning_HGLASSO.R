
##############################################
##############################################
###################################### HGLASSO
##############################################
##############################################


##############################################
## 
#################################################
#################################################
## hgl.select:
##  Finds a set of tuning parameters (rho1, rho2, rho3) 
##  that are optimal for a simulation setting. This
##  is used later as pretraining parameters to reduce
##  the computational time of the HGL in the systematic
##  simulations.
##
##  INPUTS
##    pm    : true underlying precision matrix. Used to
##              generate sample for pretraining.
##    p     : number of variables.
##    n     : sample size.
##    cov   : if FALSE, pretraining on correlation matrix.
##    
##  OUTPUT:
##    .rho  : tuning parameter value for which HWGL method
##              is diagonal.
## 
hgl.select <- function(
  pm, p, n, cval = 0.5, length.out = NULL, cov = TRUE) {
  
  if (is.null(length.out)) {
    .length.out <- 10
  } else {
    .length.out <- length.out
  }
  ######################
  ######################
  ## Step 1: Generate data:
  .sigma <- solve(pm)
  .X <- rmvnorm(n = n, sigma = .sigma, method = "svd")
  
  ######################
  ######################
  ## Step 2: find a local 
  ## optimal lambda1,lambda2
  .cov <- NULL
  if (cov) .cov <- cov(.X)
  if (!cov) .cov <- cor(.X)
  
  .HGL <- BIChgl.reduced2(
    mat = .cov, rho1length = 10, rho2length = 10, rho3length = 10,
    cval = cval, p = p, n = n)
  .optimal.rho1 <- .HGL$optimal.rho[1]
  .optimal.rho2 <- .HGL$optimal.rho[2]
  .optimal.rho3 <- .HGL$optimal.rho[3]
  

  .output <- list(
    rho1 = .optimal.rho1, rho2 = .optimal.rho2,  rho3 = .optimal.rho3)
  return(.output)
  
}


