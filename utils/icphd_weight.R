sig_prop = function(p, s){
  w = matrix(0, nrow = p, ncol = p)
  
  for(i in 2:p){
    for(j in 1:(i-1)){
      w[i,j] = s[i] * s[j]
    }
  }
  w = w + t(w)
  return(w)
}

rev_prop = function(p, w_inf){
  w = matrix(0, nrow = p, ncol = p)
  s = 1 - w_inf
  for(i in 2:p){
    for(j in 1:(i-1)){
      w[i,j] = s[i] * s[j]
    }
  }
  w = w + t(w)
  return(w)
}


inverse_prop = function(p, w_inf){
  w = matrix(0, nrow = p, ncol = p)
  
  for(i in 2:p){
    for(j in 1:(i-1)){
      w[i,j] = 1 / (w_inf[i] * w_inf[j])
    }
  }
  w = w + t(w)
  return(w)
}

inv_harmonic_weigh = function(p, w_inf){
  w = matrix(0, nrow = p, ncol = p)
  
  for(i in 2:p){
    for(j in 1:(i-1)){
      w[i,j] = (w_inf[i] + w_inf[j]) / ( 2 * w_inf[i] * w_inf[j])
    }
  }
  w = w + t(w)
  return(w)
}


inv_max = function(p, w_inf){
  w = matrix(0, nrow = p, ncol = p)
  
  for(i in 2:p){
    for(j in 1:(i-1)){
      w[i,j] = 1 / max(w_inf[i], w_inf[j])
    }
  }
  w = w + t(w)
  return(w)
}

rel_prod = function(p, w_inf){
  w = matrix(0, nrow = p, ncol = p)
  
  for(i in 2:p){
    for(j in 1:(i-1)){
      w[i,j] = 1 / (w_inf[i] * w_inf[j] * (1 - abs(w_inf[i] - w_inf[j]) / (w_inf[i] + w_inf[j])))
    }
  }
  w = w + t(w)
  w = w / quantile(w, 0.95)
  return(w)
}


rel_prod_sig = function(p, w_inf){
  w = rel_prod(p, w_inf)
  w = sigmoid(w, a = pi / sqrt(3 * sd(w) ** 2), b = mean(w))
  return(w)
}