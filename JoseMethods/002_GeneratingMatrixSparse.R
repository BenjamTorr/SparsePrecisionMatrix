#################################################
#################################################
#################################################
##
## In the following document, we will
## define a method for generating positive
## definite matrices with hubs, which is 
## based on first generating the eigen-
## structure.
##
#################################################
#################################################


#################################################
#################################################
## .adjmat: 
##    Auxiliary function. Generates a
##    random adjacency matrix of a network with 
##    hubs. Uses non-uniform Erdos Renyi model,
##    where hubs have higher connection
##    probability.
##
##  INPUTS
##    p     : Total dimension of pxp matrix.
##    T0    : Size of T0xT0 highly connected submatrix.
##    r     : number of hubs.
##    ph    : Probability of hub connectivity 
##    pnh   : Probability of connection for variables 
##              in the high-connection T0xT0 submatrix.
##    pneff : Probability of connection for all
##              variables outside of the T0xT0 
##              sumatrix.
##
##  OUTPUT
##    .A    : pxp matrix with {0,1} entries. 
##
.adjmat <- function(p, T0, r, ph, pnh, pneff) {
  .A <- matrix(rep(0, p * p), ncol = p)

  .A[] <- rbinom(p * p,  1, pneff)
  .A[1:T0, (r+1):T0] <- rbinom(T0 * (T0 - r),  1, pnh)
  .A[1:T0, 1:r] <- rbinom(T0 * r,  1, ph)
  .A[upper.tri(.A, TRUE)] <- 0
  .A <- .A + t(.A)
  return(.A)
}

.rsign <- function(n = 1) {
  return(2 * rbinom(n, 1, 0.5) - 1)
}

#################################################
#################################################
## .rsign_unif:
##    Auxiliary function. Generates a 
##    random variable with distribution:
##    Unif( [-max,-min]U[min,max] )
##
##  INPUTS
##    n       : size of output vector.
##    min     : minimum absolute value of variables. 
##    max     : maximum absolute value of variables.
##
##  OUTPUT
##    .output : length n vector with entries of 
##                distribution Unif([-max,-min]U[min,max])
## 
.rsign_unif <- function(n = 1 , min = 4, max = 5) {
  .cont <- runif(n, min, max)
  .sign <- .rsign(n)
  return(.cont * .sign)
}

#################################################
#################################################
## .rsymmmatrix: 
##    Auxiliary function. Generates a
##    random matrix with signed-uniform or normal
##    entries. The distribution of the entries 
##    depends on whether variables are hubs, 
##    highly connected or of low connection.
##
##  INPUTS
##    p       : Total dimension of pxp matrix.
##    T0      : Size of T0xT0 highly connected submatrix.
##    r       : number of hubs.
##    type    : whether the distribution signed-uniform
##                or zero-mean normal.
##    hmin    : Used when type = "unif". Minimum absolute 
##                value of hub entries.
##    hmax    : Used when type = "unif". Maximum absolute 
##                value of hub entries. 
##    nhmin   : Used when type = "unif". Minimum absolute 
##                value of entries in T0xT0 submatrix 
##                not related to hub variables.
##    nhmax   : Used when type = "unif". Maximum absolute 
##                value of entries in T0xT0 submatrix 
##                not related to hub variables.
##    neffmin : Used when type = "unif". Minimum absolute 
##                value of entries outside of the 
##                T0xT0 submatrix 
##    neffmin : Used when type = "unif". Minimum absolute 
##                value of entries outside of the 
##                T0xT0 submatrix 
##    hsd     : Used when type = "gaussian". Standard
##                deviation of normal hub entries.
##    nhsd    : Used when type = "gaussian". Standard
##                deviation of normal entries in T0xT0 
##                submatrix not related to hub variables.
##    nhsd    : Used when type = "gaussian". Standard
##                deviation of normal entries outside of 
##                the T0xT0 submatrix 
##
##  OUTPUT
##    .theta  : random pxp matrix with entrywise normal or
##                signed uniform distribution
##                depending on given parameters.     
##
.rsymmmatrix <- function(p, T0, r, 
                         type = c("unif", "gaussian"),
                         hmin = 0.5, hmax = 0.8,
                         nhmin = 0.5, nhmax = 0.8,
                         neffmin = 0.5, neffmax = 0.8,
                         hsd = 1, nhsd = 1, neffsd = 1){
  
  .theta <- matrix(rep(0, p * p), ncol = p)
  if (type == "unif") {
    .theta[] <- .rsign_unif(p * p, min = neffmin, max = neffmax)
    .theta[1:T0, (r+1):T0] <- .rsign_unif(T0 * (T0 - r), min = nhmin, max = nhmax)
    .theta[1:T0, 1:r] <- .rsign_unif(T0 * r,  min = hmin, max = hmax)
    .theta[upper.tri(.theta, TRUE)] <- 0 
  } 
  if(type == "gaussian"){
    .theta[] <- rnorm(p * p, sd = neffsd)
    .theta[1:T0, (r+1):T0] <- rnorm(T0 * (T0 - r), sd = nhsd)
    .theta[1:T0, 1:r] <- rnorm(T0 * r,  sd = hsd)
    .theta[upper.tri(.theta, TRUE)] <- 0
  }

  .theta <- .theta + t(.theta)
  return(.theta)
}

#################################################
#################################################
## .shufflemat: 
##    Auxiliary function. Shuffles the 
##    position of the hubs from the
##    first couple entries to any 
##    random entry.
##
##  INPUTS
##    p   : Total dimension of pxp matrix.
##    A   : pxp matrix.
##
##  OUTPUT
##    .M  : matrix with randomly shuffled columns/
##            rows of A.
##
.shufflemat <- function(A, p) {
  .neworder <- sample(1:p, p, replace = FALSE)
  .M <- A[.neworder, .neworder]
  return(.M)
}

#################################################
#################################################
## .rhubmat: 
##    Auxiliary function. Generates a 
##    random sparse symmetric matrix with hubs, where 
##    the non-zero entries are either signed-uniform 
##    or Gaussian. 
##    The matrix which contains a  T0xT0 matrix 
##    containing hub variables and high connectivity. 
##    The entries outside the T0xT0 submatrix are 
##    very sparse.
##
##  INPUTS
##    p       : Total dimension of pxp matrix.
##    T0      : Size of T0xT0 highly connected submatrix.
##    r       : number of hubs.
##    ph      : Probability of hub connectivity 
##    pnh     : Probability of connection for variables 
##                in the high-connection T0xT0 submatrix.
##    pneff   : Probability of connection for all
##                variables outside of the T0xT0 
##                sumatrix.
##    shuffle : If true, shuffles rows/columns to a random
##                position.
##    type    : whether the distribution signed-uniform
##                or zero-mean normal.
##    hmin    : Used when type = "unif". Minimum absolute 
##                value of hub entries.
##    hmax    : Used when type = "unif". Maximum absolute 
##                value of hub entries. 
##    nhmin   : Used when type = "unif". Minimum absolute 
##                value of entries in T0xT0 submatrix 
##                not related to hub variables.
##    nhmax   : Used when type = "unif". Maximum absolute 
##                value of entries in T0xT0 submatrix 
##                not related to hub variables.
##    neffmin : Used when type = "unif". Minimum absolute 
##                value of entries outside of the 
##                T0xT0 submatrix 
##    neffmin : Used when type = "unif". Minimum absolute 
##                value of entries outside of the 
##                T0xT0 submatrix 
##    hsd     : Used when type = "gaussian". Standard
##                deviation of normal hub entries.
##    nhsd    : Used when type = "gaussian". Standard
##                deviation of normal entries in T0xT0 
##                submatrix not related to hub variables.
##    nhsd    : Used when type = "gaussian". Standard
##                deviation of normal entries outside of 
##                the T0xT0 submatrix 
##
##  OUTPUT
##    .theta  : random sparse pxp matrix with hubs, and a T0xT0
##                highly connected matrix. Non-zero entries 
##                are either normally distributed or
##                signed uniform distribution depending on 
##                given parameters.     
##
.rhubmat <- function(p, T0, r, ph, pnh, pneff, shuffle = FALSE,
                     type = c("unif", "gaussian"),
                     hmin = 0.5, hmax = 0.8,
                     nhmin = 0.5, nhmax = 0.8,
                     neffmin = 0.5, neffmax = 0.8,
                     hsd = 1, nhsd = 1, neffsd = 1){
  .A <- .adjmat(p, T0, r, ph, pnh, pneff)
  .theta <- .A * .rsymmmatrix(p, T0, r,
                             type = type,
                             hmin = hmin, hmax = hmax,
                             nhmin = nhmin, nhmax = nhmax,
                             neffmin = neffmin, neffmax = neffmax,
                             hsd = hsd, nhsd = nhsd, neffsd = neffsd)
  if (shuffle) {
    .theta <- .shufflemat(.theta, p)
  }
  return(.theta)
}

#################################################
#################################################
## r.sparse.pdhubmat: 
##    Generates a random sparse positive definite
##    matrix with hubs, where the non-zero entries
##    are either signed-uniform or Gaussian. 
##    The matrix which contains a  T0xT0 matrix 
##    containing hub variables and high connectivity. 
##    The entries outside the T0xT0 submatrix are 
##    very sparse.
##
##  INPUTS
##    p       : Total dimension of pxp matrix.
##    T0      : Size of T0xT0 highly connected submatrix.
##    r       : number of hubs.
##    ph      : Probability of hub connectivity 
##    pnh     : Probability of connection for variables 
##                in the high-connection T0xT0 submatrix.
##    pneff   : Probability of connection for all
##                variables outside of the T0xT0 
##                sumatrix.
##    diagonal_shift : 
##                Size of minimum eigenvalue of the 
##                output matrix.
##    shuffle : If true, shuffles rows/columns to a random
##                position.
##    type    : whether the distribution signed-uniform
##                or zero-mean normal.
##    hmin    : Used when type = "unif". Minimum absolute 
##                value of hub entries.
##    hmax    : Used when type = "unif". Maximum absolute 
##                value of hub entries. 
##    nhmin   : Used when type = "unif". Minimum absolute 
##                value of entries in T0xT0 submatrix 
##                not related to hub variables.
##    nhmax   : Used when type = "unif". Maximum absolute 
##                value of entries in T0xT0 submatrix 
##                not related to hub variables.
##    neffmin : Used when type = "unif". Minimum absolute 
##                value of entries outside of the 
##                T0xT0 submatrix 
##    neffmin : Used when type = "unif". Minimum absolute 
##                value of entries outside of the 
##                T0xT0 submatrix 
##    hsd     : Used when type = "gaussian". Standard
##                deviation of normal hub entries.
##    nhsd    : Used when type = "gaussian". Standard
##                deviation of normal entries in T0xT0 
##                submatrix not related to hub variables.
##    nhsd    : Used when type = "gaussian". Standard
##                deviation of normal entries outside of 
##                the T0xT0 submatrix 
##
##  OUTPUT
##    .theta  : random sparse and positive definite pxp matrix 
##                with hubs, and a T0xT0
##                highly connected matrix. Non-zero entries 
##                are either normally distributed or
##                signed uniform distribution depending on 
##                given parameters.     
##
r.sparse.pdhubmat <- function(p, T0, r, ph, pnh, pneff,
                              diagonal_shift = 1, 
                              shuffle = FALSE,
                              type = c("unif", "gaussian"),
                              hmin = 0.5, hmax = 0.8,
                              nhmin = 0.5, nhmax = 0.8,
                              neffmin = 0.5, neffmax = 0.8,
                              hsd = 1, nhsd = 1, neffsd = 1,
                              verbose = FALSE, ...){
  .pm <- .rhubmat(p, T0, r, ph, pnh, pneff, shuffle = shuffle,
                 type = type,
                 hmin = hmin, hmax = hmax,
                 nhmin = nhmin, nhmax = nhmax,
                 neffmin = neffmin, neffmax = neffmax,
                 hsd = hsd, nhsd = nhsd, neffsd = neffsd)

  .lambda <- eigen(.pm)$value[p]
  .I <- (-.lambda + diagonal_shift) * diag(p)

  if (verbose) {
    print(paste("The dimension of PM is:", dim(.pm)[1], "x", dim(.pm)[2]))
    print(paste("The dimension of .I is:", dim(.I)[1], "x", dim(.I)[2]))
  }
  
  .pm <- .pm + .I
  
  return(.pm)
}

#################################################
#################################################
#################################################
#################################################

example <- FALSE

if (example) {
  A <- .rhubmat(
    10, 5, 2, 0, 0, 1, FALSE, "unif",
    1, 1.1,
    1, 1.1,
    1, 1.1)

  round(A, 0)

  A <- r.sparse.pdhubmat(
    10, 5, 2, 0, 0, 1, 20,   FALSE, "unif",
    1, 1.1,
    1, 1.1,
    1, 1.1)
  eigen(A)$values[10]
  round(A, 0)
  
}
rm(example)

