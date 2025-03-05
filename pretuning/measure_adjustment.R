library(mclust)


find_param_sig = function(measure, n_hub, n_nh, tol = 0.05){
  a = pi / sqrt(3 * sd(measure) ** 2)
  pos = n_hub + 1 
  b = (measure[order(measure)[pos]] + measure[order(measure)[n_hub]]) /2
  factor = 1
  n_nh_hat = 0
  n_hub_hat = 0
  while(factor < 10 && n_nh >= n_nh_hat && n_hub >= n_hub_hat){
    transf = sigmoid(measure, a * factor, b)
    n_nh_hat = sum(transf >= 1 - tol)
    n_hub_hat = sum(transf <= tol)
    factor = factor + 0.01
  }
  return(list(a = a, b = b, s = transf))
}

linear_trans = function(measure, n_hub, n_nh){
  measure = measure + 1e-4
  measure_order = order(measure)
  pos_hub = measure_order[n_hub]
  pos_nh = tail(measure_order, n_nh)[1]
  
  new_measure = (measure - measure[pos_hub]) / measure[pos_nh]
  if(any(is.na(new_measure)) | any(is.infinite(measure))){
    new_measure = measure
  }
  new_measure[measure_order[1:n_hub]] = 0
  p = length(measure)
  new_measure[measure_order[(p - n_nh + 1):p]] = 1
  
  return(new_measure)
}


# get_linear_tr = function(measure){
#   GMM = suppressWarnings(invisible(Mclust(measure)))
#   optimal.size = min(GMM$G,3)
#   if(optimal.size == 1){
#     optimal.size = 2
#   }
#   fit = suppressWarnings(invisible(Mclust(measure, G=optimal.size)))
#   labels = fit$classification
#   first_group = measure[labels == 1]
#   second_group = measure[labels == 2]
#   groups = list("1" = first_group, "2" = second_group)
#   if(optimal.size == 2){
#     means_order = order(c(mean(first_group), mean(second_group)))
#     n_hub = length(groups[[means_order[1]]])
#     n_nh = length(groups[[means_order[2]]])
#   }else{
#     third_group = measure[labels == 3]
#     groups[[3]] = third_group
#     means_order = order(c(mean(first_group), mean(second_group), mean(third_group)))
#     n_hub = length(groups[[means_order[1]]])
#     n_nh = length(groups[[means_order[3]]])
#   }
#   return(linear_trans(measure, n_hub, n_nh))
# }


get_linear_tr = function(measure){
  mean_s = mean(measure)
  sd_s = sd(measure)
  n_hub = sum((measure <= (mean_s - 2 * sd_s)))
  n_nh = sum((measure >= mean_s - sd_s))
  return(linear_trans(measure, n_hub, n_nh))
}

get_sig_tr = function(measure){
  GMM = suppressWarnings(invisible(Mclust(measure)))
  optimal.size = min(GMM$G,3)
  if(optimal.size == 1){
    optimal.size = 2
  }
  fit = suppressWarnings(invisible(Mclust(measure, G=optimal.size)))
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