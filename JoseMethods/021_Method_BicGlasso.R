
##############################################
## BIC function for HWGLASSO and GLASSO:
original <- TRUE


#################################################
#################################################
## .BICemp:
##   Bayesian Information criteria for GGM estimation
##    as in Gao et al. (2012)
##
##  INPUTS:
##      mat       : sample covariance/correlation matrix.
##      inv_est   : estimated precision matrix/inverse correlation.
##      p         : total dimension.
##      n         : sample size.
##      threshold : degree count in BIC only considers 
##                    absolute entries over threshold.
##
##  OUTPUTS: 
##      .bic : value of the BIC.
##
.BICemp <- function(mat, inv_est, p, n, threshold = 0){
  ## After symmetrization, there is a legitimate minimizer. Consider this...
  inv_est_symm  <- (inv_est + t(inv_est)) / 2 
  
  .eig          <- eigen(inv_est_symm)
  .eig_val      <- Re(.eig$values)
  .eig_val[.eig_val <= 0] <- 1e-5
  
  .bic <- 0
  .bic <- .bic - n * sum( log(.eig_val) )
  .bic <- .bic + n * Trace(inv_est %*% mat)
  
  .edgenum <- sum(abs(inv_est[upper.tri(inv_est, FALSE)]) > threshold)
  
  .bic <- .bic + log(n) * .edgenum
  return(.bic)
}


.Likemp <- function(mat, inv_est){
  ## After symmetrization, there is a legitimate minimizer. Consider this...
  inv_est_symm  <- (inv_est + t(inv_est)) / 2 
  
  .eig          <- eigen(inv_est_symm)
  .eig_val      <- Re(.eig$values)
  .eig_val[.eig_val <= 0] <- 1e-5
  
  .lik <- 0
  .lik <- .lik + sum( log(.eig_val) )
  .lik <- .lik - Trace(inv_est %*% mat)
  
  return(.lik)
}
