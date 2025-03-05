

##############################################
##############################################
#################################### HW.GLASSO
##############################################
##############################################
##
## The tuning parameter range can be estimated
## by using the fact that GLASSO = Diagmat
## if the penalty > max(|Sigma|).
##
##



#################################################
#################################################
## tp.hwgl:
##  Finds a tuning parameter for which the output of the
##  HWGL method is diagonal. This bounds the range
##  of values for training the HWGL method.  
##
##  INPUTS
##    X     : nxp data matrix.
##    pm    : true underlying precision matrix. Used to
##              generate sample for pretraining. Used 
##              only if X = NULL.
##    p     : number of variables.
##    n     : sample size.
##    cov   : if FALSE, pretraining on correlation matrix.
##    
##  OUTPUT:
##    .rho  : tuning parameter value for which HWGL method
##              is diagonal.
## 
tp.weightedGL_pre <- function(X = NULL, pm = NULL, .W ,p, n, cov = TRUE) {
  ##############################################
  ##############################################
  ## Step 1: Generate data/ use the given data:
  if (is.null(X)) {
    .sigma <- solve(pm)
    .X <- mvrnorm(n = n, mu = rep(0, p), Sigma = solve(pm))
  } else {
    .X <- X
  }
  
  ## Find covariance.
  .cov <- NULL
  .rho <- NULL
  if (cov) {
    .cov <- cov(.X)
  } else {
    .cov <- cor(.X)
  } 
  .rho <- 10^30

  ##############################################
  ##############################################
  ## Estimate the value for which the model is 
  ## first diagonal.
  .output.hwgl <- glasso(
    s = .cov, rho = .rho*.W, nobs = n,
    maxit = 150, penalize.diagonal = FALSE)
  .is.diagonal <- is.diagonal.matrix(.output.hwgl$wi)
  
  .count <- 1
  while (!.is.diagonal & .count < 5) {
    .count <- .count + 1
    .rho <- 10^5 * .rho
    if(!is.finite(.rho) | is.na(.rho)){
      return(3)
      break
    }else{
    .output.hwgl <- glasso(
      s = .cov, rho = .rho*.W, nobs = n,
      maxit = 150, penalize.diagonal = FALSE)
    .is.diagonal <- is.diagonal.matrix(.output.hwgl$wi)
    }
  }
  if (!.is.diagonal) {
    print("Warning: No diagonal/non-diagonal threshold was found: lower bound only.")
    return(3)
  }
  
  .count <- 1
  while (.is.diagonal & .count < 201) {
    ## Reduce the value of the tuning parameter
    .count <- .count + 1
    .rho <- .rho / 2
    if(!is.finite(.rho) | is.na(.rho)){
      return(3)
     break 
    }else{
      ## find the new model:
      .output.hwgl <- glasso(
        s = .cov, rho =  .rho * .W, nobs = n,
        maxit = 150, penalize.diagonal = FALSE)
      .is.diagonal <- is.diagonal.matrix(.output.hwgl$wi)
    }
  }   
  if (.is.diagonal) {
    print("Warning: No diagonal/non-diagonal threshold was found: upper bound only.")
    return(.rho)
  }  
  
  #print(paste("The diagonal/non-diagonal treshold is:", .rho))
  return(.rho)
}

tp.weightedGL <- function(X = NULL, pm = NULL, .W ,p, n, cov = TRUE){
  .rho =   tryCatch(
    tp.weightedGL_pre(X, pm, .W, p, n, cov),
    error = function(e) {
      return(3)
    }
  )
}



