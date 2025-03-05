#################################################
#################################################
#################################################
#################################################
## Auxiliary functions useful for jumping from
## Cov, corr,PM and IC.


#################################################
#################################################
## .PMtoIC: 
##    Auxiliary function. Given a precision
##    matrix Theta, calculates its corresponding
##    inverse correlation matrix.
##
##  INPUTS
##    Theta     : pxp precision matrix.
##
##  OUTPUT
##    .inv_cor  : pxp inverse correlation matrix 
##                  associated with Theta.
##
.PMtoIC <- function(Theta) {
  .Sigma <- solve(Theta)
  .diag <- diag(.Sigma)
  .diaginvsq <- 1 / sqrt(.diag)

  .rho <- diag(.diaginvsq) %*% .Sigma %*% diag(.diaginvsq)
  .inv_cor <- solve(.rho)

  return(.inv_cor)
}

#################################################
#################################################
## .COVtoIC:
##    Auxiliary function. Given a covariance
##    matrix Sigma, calculates its corresponding
##    inverse correlation matrix.
##
##  INPUTS
##    Sigma     : pxp covariance matrix.
##
##  OUTPUT
##    .inv_cor  : pxp inverse correlation matrix 
##                  associated with Sigma.
##
.COVtoIC <- function(Sigma) {
  .diag <- diag(Sigma)
  .diaginvsq <- 1 / sqrt(.diag)

  .rho <- diag(.diaginvsq) %*% Sigma %*% diag(.diaginvsq)
  .inv_cor <- solve(.rho)

  return(.inv_cor)
}

#################################################
#################################################
## .COVtoCOR:
##    Auxiliary function. Given a covariance
##    matrix Sigma, calculates its corresponding
##    inverse correlation matrix.
##
##  INPUTS
##    Sigma : pxp covariance matrix.
##
##  OUTPUT
##    .rho  : pxp inverse correlation matrix 
##                  associated with Sigma.
##
.COVtoCOR <- function(Sigma) {
  .diag <- diag(Sigma)
  .diaginvsq <- 1 / sqrt(.diag)

  .rho <- diag(.diaginvsq) %*% Sigma %*% diag(.diaginvsq)
  return(.rho)
}


#################################################
#################################################
## Trace: 
##    Trace of a square matrix.
##
##  INPUTS
##    M : square matrix.
##
##  OUTPUT
##    x  : trace of M.
##
Trace = function(M){
  return(sum(diag(M)))
}