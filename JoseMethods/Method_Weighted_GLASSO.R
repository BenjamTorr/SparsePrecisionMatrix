
#################################################
#################################################
## BICweighted_glasso:
##  This function computes the HW-GLASSO estimator
##  for the given range of tuning parameters, 
##  and then compute the BIC of each estimator.
##  Outputs the fit with optimal BIC Tun. parameter.
##
##  INPUTS:
##      mat       : sample covariance/correlation matrix.
##      rho       : tuning parameter range vector.
##      .w        : Weights to be used in the fitting
##      p         : total dimension.
##      n         : sample size.
##      threshold : degree count in BIC only considers 
##                    absolute entries over threshold.
##      maxit     : Maximum number of iterations of GLASSO.
##      penalize.diagonal :
##                    if FALSE, off-diagonal L1 penalty.
##
##  OUTPUTS: 
##    OUTPUT        : list of GLASSO outputs for tun. pars. in rho.
##    BIC           : vector of BIC values for tun. par. in rho.
##    optimal.index : index corresponding to optimal rho.
##    optimal.rho   : optimal tuning parameter value. 
##    optimal.model : object glasso trained with optimal.rho
##    total.time    : total time for training.
##
BICweighted_glasso <- function(mat, rho, .W, p, n,
                        threshold = 0, maxit = 200,
                        penalize.diagonal = TRUE){
  pb <- txtProgressBar(min = 0, max = length(rho), style = 3)
  
  ## Fit weigthed GLASSO models:
  .rholength = length(rho)
  .BIC.HWGL = rep(NA, .rholength)
  .OUTPUT.HWGL = list()
  .start_time = Sys.time()
  for(.i in 1:.rholength){
    setTxtProgressBar(pb, .i)
    if(.i == 1){
      .OUTPUT.HWGL[[.i]] <- glasso(
        s = mat, rho = rho[.i]*.W,
        nobs = n, zero = NULL,
        thr = 1.0e-4, maxit = maxit,  approx = FALSE,
        penalize.diagonal = penalize.diagonal,
        start = "cold", w.init = NULL, wi.init = NULL,
        trace = FALSE)
      
    } else {
      .OUTPUT.HWGL[[.i]] <- glasso(
        s = mat, rho = rho[.i]*.W,
        nobs = n, zero = NULL, 
        thr = 1.0e-4, maxit = maxit,  approx = FALSE,
        penalize.diagonal = penalize.diagonal, 
        start = "warm",
        w.init  = .OUTPUT.HWGL[[.i-1]]$w, 
        wi.init = .OUTPUT.HWGL[[.i-1]]$wi,
        trace = FALSE)
      
    }
    .BIC.HWGL[.i] <- .BICemp(
      mat = mat, 
      inv_est = .OUTPUT.HWGL[[.i]]$wi, 
      n = n, p = p, threshold = threshold)
    
  }
  ## Select optimal model:
  .end_time = Sys.time()
  .total_time_HWGL <- .end_time - .start_time
  
  ## Return outputs:
  .OUTPUT = list(OUTPUT = .OUTPUT.HWGL, 
                 BIC = .BIC.HWGL, 
                 optimal.index = which.min(.BIC.HWGL),
                 optimal.rho = rho[which.min(.BIC.HWGL)],
                 optimal.model = .OUTPUT.HWGL[[which.min(.BIC.HWGL)]],
                 total.time = .total_time_HWGL)
  close(pb)
  return(.OUTPUT)
}




Lik_weighted_glasso <- function(mat_train, mat_val, rho, .W, p, n,
                               threshold = 0, maxit = 200,
                               penalize.diagonal = TRUE){
  
  pb <- txtProgressBar(min = 0, max = length(rho), style = 3)
  ## Fit weigthed GLASSO models:
  .rholength = length(rho)
  .BIC.HWGL = rep(NA, .rholength)
  .OUTPUT.HWGL = list()
  
  .start_time = Sys.time()
  for(.i in 1:.rholength){
    setTxtProgressBar(pb, .i)
    if(.i == 1){
      .OUTPUT.HWGL[[.i]] <- glasso(
        s = mat_train, rho = rho[.i]*.W,
        nobs = n, zero = NULL,
        thr = 1.0e-4, maxit = maxit,  approx = FALSE,
        penalize.diagonal = penalize.diagonal,
        start = "cold", w.init = NULL, wi.init = NULL,
        trace = FALSE)
      
    } else {
      .OUTPUT.HWGL[[.i]] <- glasso(
        s = mat_train, rho = rho[.i]*.W,
        nobs = n, zero = NULL, 
        thr = 1.0e-4, maxit = maxit,  approx = FALSE,
        penalize.diagonal = penalize.diagonal, 
        start = "warm",
        w.init  = .OUTPUT.HWGL[[.i-1]]$w, 
        wi.init = .OUTPUT.HWGL[[.i-1]]$wi,
        trace = FALSE)
      
    }
    .BIC.HWGL[.i] <- -.Likemp(
      mat = mat_val,
      inv_est = .OUTPUT.HWGL[[.i]]$wi)
    
  }
  ## Select optimal model:
  .end_time = Sys.time()
  .total_time_HWGL <- .end_time - .start_time
  
  ## Return outputs:
  .OUTPUT = list(OUTPUT = .OUTPUT.HWGL, 
                 BIC = .BIC.HWGL, 
                 optimal.index = which.min(.BIC.HWGL),
                 optimal.rho = rho[which.min(.BIC.HWGL)],
                 optimal.model = .OUTPUT.HWGL[[which.min(.BIC.HWGL)]],
                 total.time = .total_time_HWGL)
  close(pb)
  return(.OUTPUT)
}


BIC_2step_glasso <- function(mat, rho, p, n,
                               threshold = 0, maxit = 200,
                               penalize.diagonal = TRUE){
  pb <- txtProgressBar(min = 0, max = length(rho), style = 3)
  
  ## Fit weigthed GLASSO models:
  .rholength = length(rho)
  .BIC.HWGL = matrix(NA, nrow = .rholength, ncol = .rholength)
  .OUTPUT.HWGL = list()
  .start_time = Sys.time()
  for(.j in 1:.rholength){  
    setTxtProgressBar(pb, .j)
    .OUTPUT.HWGL[[.j]] = list()
    for(.i in 1:.rholength){
      W = matrix(rho[.j], args$p, args$p)
      W[1:args$r, ] = rho[.i]
      W[, 1:args$r] = rho[.i]
      diag(W) = 0
      if(.i == 1){
        .OUTPUT.HWGL[[.j]][[.i]] <- glasso(
          s = mat, rho = W,
          nobs = n, zero = NULL,
          thr = 1.0e-4, maxit = maxit,  approx = FALSE,
          penalize.diagonal = penalize.diagonal,
          start = "cold", w.init = NULL, wi.init = NULL,
          trace = FALSE)
        
      } else {
        .OUTPUT.HWGL[[.j]][[.i]] <- glasso(
          s = mat, rho = W,
          nobs = n, zero = NULL, 
          thr = 1.0e-4, maxit = maxit,  approx = FALSE,
          penalize.diagonal = penalize.diagonal, 
          start = "warm",
          w.init  = .OUTPUT.HWGL[[.j]][[.i-1]]$w, 
          wi.init = .OUTPUT.HWGL[[.j]][[.i-1]]$wi,
          trace = FALSE)
        
      }
      .BIC.HWGL[.j, .i] <- .BICemp(
        mat = mat, 
        inv_est = .OUTPUT.HWGL[[.j]][[.i]]$wi, 
        n = n, p = p, threshold = threshold)
      
    }
  }
  ## Select optimal model:
  .end_time = Sys.time()
  .total_time_HWGL <- .end_time - .start_time
  pos = which(.BIC.HWGL == min(.BIC.HWGL), arr.ind = TRUE)
  pos_matrix = matrix(pos, nrow = length(pos)/2, ncol= 2) #create matrix of rows and columns
  
  ## Return outputs:
  .OUTPUT = list(OUTPUT = .OUTPUT.HWGL, 
                 BIC = .BIC.HWGL, 
                 optimal.index = c(pos_matrix[1,1], pos_matrix[1,2]),
                 optimal.rho = c(rho[pos_matrix[1,1]], rho[pos_matrix[1,2]]),
                 optimal.model = .OUTPUT.HWGL[[pos_matrix[1,1]]][[pos_matrix[1,2]]],
                 total.time = .total_time_HWGL)
  close(pb)
  return(.OUTPUT)
}


Lik_2step_glasso <- function(mat, mat_val, rho, p, n,
                             threshold = 0, maxit = 200,
                             penalize.diagonal = TRUE){
  pb <- txtProgressBar(min = 0, max = length(rho), style = 3)
  
  ## Fit weigthed GLASSO models:
  .rholength = length(rho)
  .BIC.HWGL = matrix(NA, nrow = .rholength, ncol = .rholength)
  .OUTPUT.HWGL = list()
  .start_time = Sys.time()
  for(.j in 1:.rholength){
    setTxtProgressBar(pb, .j)
    .OUTPUT.HWGL[[.j]] = list()
    for(.i in 1:.rholength){
      W = matrix(rho[.j], args$p, args$p)
      W[1:args$r, ] = rho[.i]
      W[, 1:args$r] = rho[.i]
      diag(W) = 0
      if(.i == 1){
        .OUTPUT.HWGL[[.j]][[.i]] <- glasso(
          s = mat, rho = W,
          nobs = n, zero = NULL,
          thr = 1.0e-4, maxit = maxit,  approx = FALSE,
          penalize.diagonal = penalize.diagonal,
          start = "cold", w.init = NULL, wi.init = NULL,
          trace = FALSE)
        
      } else {
        .OUTPUT.HWGL[[.j]][[.i]] <- glasso(
          s = mat, rho = W,
          nobs = n, zero = NULL, 
          thr = 1.0e-4, maxit = maxit,  approx = FALSE,
          penalize.diagonal = penalize.diagonal, 
          start = "warm",
          w.init  = .OUTPUT.HWGL[[.j]][[.i-1]]$w, 
          wi.init = .OUTPUT.HWGL[[.j]][[.i-1]]$wi,
          trace = FALSE)
        
      }
      .BIC.HWGL[.j, .i] <- -.Likemp(
        mat = mat_val, 
        inv_est = .OUTPUT.HWGL[[.j]][[.i]]$wi)
      
    }
  }
  ## Select optimal model:
  .end_time = Sys.time()
  .total_time_HWGL <- .end_time - .start_time
  pos = which(.BIC.HWGL == min(.BIC.HWGL), arr.ind = TRUE)
  pos_matrix = matrix(pos, nrow = length(pos)/2, ncol= 2) #create matrix of rows and columns
  
  ## Return outputs:
  .OUTPUT = list(OUTPUT = .OUTPUT.HWGL, 
                 BIC = .BIC.HWGL, 
                 optimal.index = c(pos_matrix[1,1], pos_matrix[1,2]),
                 optimal.rho = c(rho[pos_matrix[1,1]], rho[pos_matrix[1,2]]),
                 optimal.model = .OUTPUT.HWGL[[pos_matrix[1,1]]][[pos_matrix[1,2]]],
                 total.time = .total_time_HWGL)
  close(pb)
  return(.OUTPUT)
}

# BICweighted_glasso_penHubs <- function(mat, rho, .W, p, n,
#                                threshold = 0, maxit = 200,
#                                penalize.diagonal = TRUE, penalize.threshold = 1e-5){
#   pb <- txtProgressBar(min = 0, max = length(rho), style = 3)
#   
#   ## Fit weigthed GLASSO models:
#   .rholength = length(rho)
#   .BIC.HWGL = rep(NA, .rholength)
#   .OUTPUT.HWGL = list()
#   .start_time = Sys.time()
#   
#   
#   for(.i in 1:.rholength){
#     setTxtProgressBar(pb, .i)
#     .W[.W <= penalize.threshold] = rho[.i]
#     if(.i == 1){
#       .OUTPUT.HWGL[[.i]] <- glasso(
#         s = mat, rho =.W,
#         nobs = n, zero = NULL,
#         thr = 1.0e-4, maxit = maxit,  approx = FALSE,
#         penalize.diagonal = penalize.diagonal,
#         start = "cold", w.init = NULL, wi.init = NULL,
#         trace = FALSE)
#       
#     } else {
#       .OUTPUT.HWGL[[.i]] <- glasso(
#         s = mat, rho = .W,
#         nobs = n, zero = NULL, 
#         thr = 1.0e-4, maxit = maxit,  approx = FALSE,
#         penalize.diagonal = penalize.diagonal, 
#         start = "warm",
#         w.init  = .OUTPUT.HWGL[[.i-1]]$w, 
#         wi.init = .OUTPUT.HWGL[[.i-1]]$wi,
#         trace = FALSE)
#       
#     }
#     .BIC.HWGL[.i] <- .BICemp(
#       mat = mat, 
#       inv_est = .OUTPUT.HWGL[[.i]]$wi, 
#       n = n, p = p, threshold = threshold)
#     
#   }
#   ## Select optimal model:
#   .end_time = Sys.time()
#   .total_time_HWGL <- .end_time - .start_time
#   
#   ## Return outputs:
#   .OUTPUT = list(OUTPUT = .OUTPUT.HWGL, 
#                  BIC = .BIC.HWGL, 
#                  optimal.index = which.min(.BIC.HWGL),
#                  optimal.rho = rho[which.min(.BIC.HWGL)],
#                  optimal.model = .OUTPUT.HWGL[[which.min(.BIC.HWGL)]],
#                  total.time = .total_time_HWGL)
#   close(pb)
#   return(.OUTPUT)
# }
# 
# Lik_weighted_glasso_penHubs <- function(mat_train, mat_val, rho, .W, p, n,
#                                 threshold = 0, maxit = 200,
#                                 penalize.diagonal = TRUE, penalize.threshold=1e-5){
#   
#   pb <- txtProgressBar(min = 0, max = length(rho), style = 3)
#   ## Fit weigthed GLASSO models:
#   .rholength = length(rho)
#   .BIC.HWGL = rep(NA, .rholength)
#   .OUTPUT.HWGL = list()
#   
#   .start_time = Sys.time()
#   for(.i in 1:.rholength){
#     setTxtProgressBar(pb, .i)
#     .W[.W <= penalize.threshold] = rho[.i]
#     if(.i == 1){
#       .OUTPUT.HWGL[[.i]] <- glasso(
#         s = mat_train, rho = .W,
#         nobs = n, zero = NULL,
#         thr = 1.0e-4, maxit = maxit,  approx = FALSE,
#         penalize.diagonal = penalize.diagonal,
#         start = "cold", w.init = NULL, wi.init = NULL,
#         trace = FALSE)
#       
#     } else {
#       .OUTPUT.HWGL[[.i]] <- glasso(
#         s = mat_train, rho = .W,
#         nobs = n, zero = NULL, 
#         thr = 1.0e-4, maxit = maxit,  approx = FALSE,
#         penalize.diagonal = penalize.diagonal, 
#         start = "warm",
#         w.init  = .OUTPUT.HWGL[[.i-1]]$w, 
#         wi.init = .OUTPUT.HWGL[[.i-1]]$wi,
#         trace = FALSE)
#       
#     }
#     .BIC.HWGL[.i] <- -.Likemp(
#       mat = mat_val,
#       inv_est = .OUTPUT.HWGL[[.i]]$wi)
#     
#   }
#   ## Select optimal model:
#   .end_time = Sys.time()
#   .total_time_HWGL <- .end_time - .start_time
#   
#   ## Return outputs:
#   .OUTPUT = list(OUTPUT = .OUTPUT.HWGL, 
#                  BIC = .BIC.HWGL, 
#                  optimal.index = which.min(.BIC.HWGL),
#                  optimal.rho = rho[which.min(.BIC.HWGL)],
#                  optimal.model = .OUTPUT.HWGL[[which.min(.BIC.HWGL)]],
#                  total.time = .total_time_HWGL)
#   close(pb)
#   return(.OUTPUT)
# }

BICweighted_glasso_penHubs <- function(mat, rho, p, n, r,
                                       threshold = 0, maxit = 200,
                                       penalize.diagonal = TRUE, penalize.threshold = 1e-5, const = 0.1){
  pb <- txtProgressBar(min = 0, max = length(rho), style = 3)
  
  ## Fit weigthed GLASSO models:
  .rholength = length(rho)
  .BIC.HWGL = rep(NA, .rholength)
  .OUTPUT.HWGL = list()
  .start_time = Sys.time()
  
  
  for(.i in 1:.rholength){
    setTxtProgressBar(pb, .i)
    .W = matrix(rho[.i], args$p, args$p)
    .W[1:r, ] = const * rho[.i]
    .W[, 1:r] = const * rho[.i]
    if(.i == 1){
      .OUTPUT.HWGL[[.i]] <- glasso(
        s = mat, rho =.W,
        nobs = n, zero = NULL,
        thr = 1.0e-4, maxit = maxit,  approx = FALSE,
        penalize.diagonal = penalize.diagonal,
        start = "cold", w.init = NULL, wi.init = NULL,
        trace = FALSE)
      
    } else {
      .OUTPUT.HWGL[[.i]] <- glasso(
        s = mat, rho = .W,
        nobs = n, zero = NULL, 
        thr = 1.0e-4, maxit = maxit,  approx = FALSE,
        penalize.diagonal = penalize.diagonal, 
        start = "warm",
        w.init  = .OUTPUT.HWGL[[.i-1]]$w, 
        wi.init = .OUTPUT.HWGL[[.i-1]]$wi,
        trace = FALSE)
      
    }
    .BIC.HWGL[.i] <- .BICemp(
      mat = mat, 
      inv_est = .OUTPUT.HWGL[[.i]]$wi, 
      n = n, p = p, threshold = threshold)
    
  }
  ## Select optimal model:
  .end_time = Sys.time()
  .total_time_HWGL <- .end_time - .start_time
  
  ## Return outputs:
  .OUTPUT = list(OUTPUT = .OUTPUT.HWGL, 
                 BIC = .BIC.HWGL, 
                 optimal.index = which.min(.BIC.HWGL),
                 optimal.rho = rho[which.min(.BIC.HWGL)],
                 optimal.model = .OUTPUT.HWGL[[which.min(.BIC.HWGL)]],
                 total.time = .total_time_HWGL)
  close(pb)
  return(.OUTPUT)
}

Lik_weighted_glasso_penHubs <- function(mat_train, mat_val, rho, p, n, r,
                                threshold = 0, maxit = 200,
                                penalize.diagonal = TRUE, const = 0.1){
  
  pb <- txtProgressBar(min = 0, max = length(rho), style = 3)
  ## Fit weigthed GLASSO models:
  .rholength = length(rho)
  .BIC.HWGL = rep(NA, .rholength)
  .OUTPUT.HWGL = list()
  
  .start_time = Sys.time()
  for(.i in 1:.rholength){
    setTxtProgressBar(pb, .i)
    .W = matrix(rho[.i], args$p, args$p)
    .W[1:r, ] = const * rho[.i]
    .W[, 1:r] = const * rho[.i]
    if(.i == 1){
      .OUTPUT.HWGL[[.i]] <- glasso(
        s = mat_train, rho = rho[.i]*.W,
        nobs = n, zero = NULL,
        thr = 1.0e-4, maxit = maxit,  approx = FALSE,
        penalize.diagonal = penalize.diagonal,
        start = "cold", w.init = NULL, wi.init = NULL,
        trace = FALSE)
      
    } else {
      .OUTPUT.HWGL[[.i]] <- glasso(
        s = mat_train, rho = rho[.i]*.W,
        nobs = n, zero = NULL, 
        thr = 1.0e-4, maxit = maxit,  approx = FALSE,
        penalize.diagonal = penalize.diagonal, 
        start = "warm",
        w.init  = .OUTPUT.HWGL[[.i-1]]$w, 
        wi.init = .OUTPUT.HWGL[[.i-1]]$wi,
        trace = FALSE)
      
    }
    .BIC.HWGL[.i] <- -.Likemp(
      mat = mat_val,
      inv_est = .OUTPUT.HWGL[[.i]]$wi)
    
  }
  ## Select optimal model:
  .end_time = Sys.time()
  .total_time_HWGL <- .end_time - .start_time
  
  ## Return outputs:
  .OUTPUT = list(OUTPUT = .OUTPUT.HWGL, 
                 BIC = .BIC.HWGL, 
                 optimal.index = which.min(.BIC.HWGL),
                 optimal.rho = rho[which.min(.BIC.HWGL)],
                 optimal.model = .OUTPUT.HWGL[[which.min(.BIC.HWGL)]],
                 total.time = .total_time_HWGL)
  close(pb)
  return(.OUTPUT)
}

# 
# BICweighted_glasso_frob <- function(Theta, mat, rho, .W, p, n,
#                                threshold = 0, maxit = 200,
#                                penalize.diagonal = TRUE){
#   pb <- txtProgressBar(min = 0, max = length(rho), style = 3)
#   
#   ## Fit weigthed GLASSO models:
#   .rholength = length(rho)
#   .BIC.HWGL = rep(NA, .rholength)
#   .OUTPUT.HWGL = list()
#   .start_time = Sys.time()
#   for(.i in 1:.rholength){
#     setTxtProgressBar(pb, .i)
#     if(.i == 1){
#       .OUTPUT.HWGL[[.i]] <- glasso(
#         s = mat, rho = rho[.i]*.W,
#         nobs = n, zero = NULL,
#         thr = 1.0e-4, maxit = maxit,  approx = FALSE,
#         penalize.diagonal = penalize.diagonal,
#         start = "cold", w.init = NULL, wi.init = NULL,
#         trace = FALSE)
#       
#     } else {
#       .OUTPUT.HWGL[[.i]] <- glasso(
#         s = mat, rho = rho[.i]*.W,
#         nobs = n, zero = NULL, 
#         thr = 1.0e-4, maxit = maxit,  approx = FALSE,
#         penalize.diagonal = penalize.diagonal, 
#         start = "warm",
#         w.init  = .OUTPUT.HWGL[[.i-1]]$w, 
#         wi.init = .OUTPUT.HWGL[[.i-1]]$wi,
#         trace = FALSE)
#       
#     }
#     .BIC.HWGL[.i] <- Frob_norm(Theta, .OUTPUT.HWGL[[.i]]$wi)
#     
#   }
#   ## Select optimal model:
#   .end_time = Sys.time()
#   .total_time_HWGL <- .end_time - .start_time
#   
#   ## Return outputs:
#   .OUTPUT = list(OUTPUT = .OUTPUT.HWGL, 
#                  BIC = .BIC.HWGL, 
#                  optimal.index = which.min(.BIC.HWGL),
#                  optimal.rho = rho[which.min(.BIC.HWGL)],
#                  optimal.model = .OUTPUT.HWGL[[which.min(.BIC.HWGL)]],
#                  total.time = .total_time_HWGL)
#   close(pb)
#   return(.OUTPUT)
# }
# 
