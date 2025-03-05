library(pracma)
library(MASS)

sample_from_precision = function (n, p, Theta, train_split = 0.8) 
{
  Sigma = solve(Theta)
  X = mvrnorm(n = n, mu = rep(0, p), Sigma = Sigma)
  train_size = floor(n * train_split)
  X_train = X[1:train_size, ]
  X_val = X[(train_size + 1):n, ]
  return(list(X = X, X_train = X_train, X_val = X_val))
}


get_cov_cor = function (X) 
{
  Sigma_hat = cov(X)
  Corr_hat = .COVtoCOR(Sigma_hat)
  #Sigma_hat = sta_thresholding_perc(X, mat_type = "cov", mat = Sigma_hat, var_inds=1:nrow(X), true_mat = NULL)$mat
  
  return(list(Cov = Sigma_hat, Corr = Corr_hat))
}