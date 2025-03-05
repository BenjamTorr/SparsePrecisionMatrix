################################################################################

rm(list=ls())
setwd("D:/UNC/PhD UNC/Research/Sparse Precision Matrix/code")

invisible(lapply(list.files("JoseMethods", pattern = "\\.R$", full.names = TRUE), source))
invisible(lapply(list.files("utils", pattern = "\\.R$", full.names = TRUE), source))
invisible(lapply(list.files("HWGLASSO", pattern = "\\.R$", full.names = TRUE), source))

source('SimulationParameters.R')
################################################################################
library(ggplot2)
library(reshape2)
library(pracma)



test_func = function(p, w_inf){
  w = matrix(0, nrow = p, ncol = p)
  
  for(i in 2:p){
    for(j in 1:(i-1)){
      a  = min(w_inf[i], w_inf[j])
      b = max(w_inf[i], w_inf[j])
      w[i,j] = 1 / (a * b)
    }
  }
  w = w + t(w)
  
  return(w / max(w))
}


Theta = r.sparse.pdhubmat(
  args$p, T0 = args$T0, r = args$r, ph = args$ph, pnh = args$pnh, pneff = args$pneff, diagonal_shift = args$diagonal_shift, 
  shuffle = args$shuffle, type = "unif", hmin = args$hmin, hmax = args$hmax, nhmin = args$nhmin, nhmax = args$nhmax, 
  neffmin = args$neffmin, neffmax = args$neffmax)


sim_data = sample_from_precision(args$n, args$p, Theta)
X = sim_data$X

CoVR_X = get_cov_cor(X)

icphd_inf = sta_ipchd(X, mat_type = 'cor', mat = CoVR_X$Corr, var_inds = 1:args$p, overest_type = 'frac')
w_icphd = test_func(args$p, icphd_inf)

plot(icphd_inf, main = 'Influence measure of correlation', ylab = 'ICPHD Measure', xlab = 'Variable')

w_icphd = w_icphd / quantile(w_icphd[50, 1:args$p], 0.95)

plot(w_icphd[1, 1:args$p], main = 'Weights for IPCHD', xlab = 'Variable', ylab = 'Weight of of entry 1,j', 
     ylim = c(min(w_icphd[50, 1:args$p]), quantile(w_icphd[50, 1:args$p], 0.95)))
points(w_icphd[50, 1:args$p], col = 'blue')
legend("topright", legend = c("Hub node - 1", "Non-hub node 50"), col = c("black", "blue"), pch = 19)


#w = sigmoid(w_icphd, a = 0.0001*pi / sqrt(3 * sd(w_icphd) ** 2), b = 3*mean(w_icphd))
w= w_icphd

plot(w[1, 1:args$p], main = 'Weights for IPCHD sigmoid', xlab = 'Variable', ylab = 'Weight of of entry 1,j', 
     ylim = c(min(w[50, 1:args$p]), quantile(w[50, 1:args$p], 0.95)))
points(w[50, 1:args$p], col = 'blue')
legend("topright", legend = c("Hub node - 1", "Non-hub node 50"), col = c("black", "blue"), pch = 19)



find_param_sig = function(measure, n_hub, n_nh){
  a = pi / sqrt(3 * sd(measure) ** 2)
  pos = n_hub + 1 
  b = (s[order(s)[pos]] + s[order(s)[n_hub]]) /2
  factor = 1
  n_nh_hat = 0
  n_hub_hat = 0
  while(factor < 10 && n_nh >= n_nh_hat && n_hub >= n_hub_hat){
    transf = sigmoid(s, a * factor, b)
    n_nh_hat = sum(transf >= 1 - 0.05)
    n_hub_hat = sum(transf <=0.05)
    factor = factor + 0.01
  }
  return(list(a = a, b = b, s = transf))
}


linear_trans = function(measure, n_hub, n_nh){
  measure_order = order(measure)
  pos_hub = measure_order[n_hub]
  pos_nh = tail(measure_order, n_nh)[1]
  
  new_measure = (measure - measure[pos_hub]) / measure[pos_nh]
  new_measure[measure_order[1:n_hub]] = 0
  p = length(measure)
  new_measure[measure_order[(p - n_nh + 1):p]] = 1
  return(new_measure)
}

get_linear_tr = function(measure){
  GMM = Mclust(measure)
  optimal.size = min(GMM$G,3)
  if(optimal.size == 1){
    optimal.size = 2
  }
  fit = Mclust(measure, G=optimal.size)
  labels = fit$classification
  first_group = measure[labels == 1]
  second_group = measure[labels == 2]
  groups = list("1" = first_group, "2" = second_group)
  if(optimal.size == 2){
    means_order = order(c(mean(first_group), mean(second_group)))
    n_hub = length(groups[[means_order[1]]])
    n_nh = length(groups[[means_order[2]]])
  }else{
    third_group = measure[labels == 3]
    groups[[3]] = third_group
    means_order = order(c(mean(first_group), mean(second_group), mean(third_group)))
    n_hub = length(groups[[means_order[1]]])
    n_nh = length(groups[[means_order[3]]])
  }
  return(linear_trans(measure, n_hub, n_nh)$s)
}

get_sig_tr = function(measure){
  GMM = Mclust(measure)
  optimal.size = min(GMM$G,3)
  if(optimal.size == 1){
    optimal.size = 2
  }
  fit = Mclust(measure, G=optimal.size)
  labels = fit$classification
  first_group = measure[labels == 1]
  second_group = measure[labels == 2]
  groups = list("1" = first_group, "2" = second_group)
  if(optimal.size == 2){
    means_order = order(c(mean(first_group), mean(second_group)))
    n_hub = length(groups[[means_order[1]]])
    n_nh = length(groups[[means_order[2]]])
  }else{
    third_group = measure[labels == 3]
    groups[[3]] = third_group
    means_order = order(c(mean(first_group), mean(second_group), mean(third_group)))
    n_hub = length(groups[[means_order[1]]])
    n_nh = length(groups[[means_order[3]]])
  }
  return(find_param_sig(measure, n_hub, n_nh)$s)
}
  





s = 1 - icphd_inf
a = 4*pi / sqrt(3 * sd(s) ** 2) 
b = s[order(s)[6]]



color <- ifelse(1:args$p %in% order(s)[1:5], "red", "blue")
values = w[upper.tri(w)]
plot(s, find_param_sig(s, 5, args$p - 5)$s, type='p', col = color, ylim=c(0,1))

plot(s, get_sig_tr(s), type="p")


#heatmap_plot <- plot_ly(
#  z = sig_prop(args$p, find_param_sig(s, 5, 45)$s),
#  type = "heatmap",
#  colorscale = "Viridis"
#)
w = sig_prop(args$p, find_param_sig(s, 5, 45)$s)
# Display the heatmap
#heatmap_plot

plot(w[1, 1:args$p], main = 'Weights for IPCHD sigmoid', xlab = 'Variable', ylab = 'Weight of of entry 1,j', 
     ylim = c(min(w[50, 1:args$p]), quantile(w[50, 1:args$p], 0.95)))
points(w[50, 1:args$p], col = 'blue')
legend("topright", legend = c("Hub node - 1", "Non-hub node 50"), col = c("black", "blue"), pch = 19)


