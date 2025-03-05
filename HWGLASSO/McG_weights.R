library(glasso)

hw_weights = function(Sigma){
  ## Calculate weighting matrix:
  p = nrow(Sigma)
  .inv_mat <- solve(Sigma + 0.1 * diag(p))
  .W  <- matrix(rep(0, p * p), ncol = p)
  influence = rep(0, p)
  influence[1] = sum(abs(.inv_mat[1, -1]))  
  
  for (.i in (2:p)) {
    for(.j in (1:(.i-1))) {
      .a  <- abs(.inv_mat[.i, .j])
      .ai <- sum(abs(.inv_mat[.i, -.i])) 
      .aj <- sum(abs(.inv_mat[.j, -.j])) 
      
      influence[.i] = .ai
      .W[.i, .j] <- (.i != .j) / (.a * .ai * .aj)
    }
  }
  .W <- .W + t(.W)
  
  .W = .W / quantile(.W, 0.95)
  return(list(w= .W,influence=influence))
}


adaGL_weights = function(Sigma){
  ## Calculate weighting matrix:
  p = nrow(Sigma)
  .inv_mat <- solve(Sigma + 0.1 * diag(p))
  .W = matrix(0, p, p)
  for (.i in (2:p)) {
    for(.j in (1:(.i-1))) {
      .a  <- abs(.inv_mat[.i, .j])
      .W[.i, .j] <- (.i != .j) / (.a)
    }
  }
  .W <- .W + t(.W)
  
  .W = .W / quantile(.W, 0.95)
  return(.W)
}
