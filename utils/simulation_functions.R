simulation_BIC_ipchd = function(Theta, args, w_func ,lambda, threshold, folder, name, save = FALSE){
  #simulate data
  set.seed(args$seed)
  sim_data = sample_from_precision(args$n, args$p, Theta)
  X = sim_data$X
  
  #get correlation and covariance matrix
  CoVR_X = get_cov_cor(X)
  
  #get influence measures
  icphd_inf = sta_ipchd(X, mat_type = 'cor', mat = CoVR_X$Corr, var_inds = 1:args$p, overest_type = 'frac')
  w_icphd = w_func(args$p, icphd_inf)
  
  #get optimal lambda

  max_rho = tp.weightedGL(X = X, .W = w_icphd, p = args$p, n = args$n)
  lambda = lambda * max_rho
  
  #Fit GLASSO
  Theta_icphd_BIC = BICweighted_glasso(mat = CoVR_X$Cov, rho = lambda, .W = w_icphd, p = args$p, n = args$n,
                                       threshold = 0, maxit = 200,
                                       penalize.diagonal = FALSE)
  #Save results to analize
  if(save){
    folder_path = paste0(folder, '/IPCHD_BIC_',name)
    
    
    
    if (!dir.exists(folder_path)) {
      dir.create(folder_path)
    }
    #############################------------ Saving Results ------------------------------------------#
    # Save Influence Measure
    
    png(paste0(folder_path,"/IPCHD_influence.png"), width = 800, height = 600)  # Set file name and dimensions
    plot(icphd_inf, main = 'Influence measure of correlation', ylab = 'ICPHD Measure', xlab = 'Variable')
    dev.off()  # Save and close the file
    
    #Save weights for a hub and non hub
    png(paste0(folder_path,"/Weight_hub_node_IPCHD.png"), width = 800, height = 600)
    plot(w_icphd[1, 1:args$p], main = paste('Weights for IPCHD', name), xlab = 'Variable', ylab = 'Weight of of entry 1,j', 
         ylim = c(min(w_icphd[args$p, 1:args$p]), quantile(w_icphd[args$p, 1:args$p], 0.95)))
    points(w_icphd[args$p, 1:args$p], col = 'blue')
    legend("topright", legend = c("Hub node - 1", "Non-hub last node"), col = c("black", "blue"), pch = 19)
    dev.off()
    
    p <- plot_ly(z = w_icphd, type = "heatmap", colorscale = "Viridis")
    
    # Save the Plotly plot to an HTML file
    html_file <- paste0(folder_path,"/Weight_matrix.html")
    htmlwidgets::saveWidget(p, html_file)
    
    # Use webshot to capture the HTML file and save as PNG
    webshot(html_file, file = paste0(folder_path,"/Weight_matrix.png"))
    #Save Sanity check for fitting method
    
    png(paste0(folder_path,"/BIC_curve_IPCHD.png"), width = 800, height = 600)
    plot(x = lambda, y = Theta_icphd_BIC$BIC, log = 'x', main ='BIC for IPCHD', ylab = 'BIC', xlab = 'Lambda')
    dev.off()
    
    #Save matrices visualizations
    
    png(paste0(folder_path,"/Absolute_diff_IPCHD.png"), width = 800, height = 600)
    plot.pdm(abs(Theta - Theta_icphd_BIC$optimal.model$wi), main = 'Absolute Difference Estimated IPCHD')
    dev.off()
    
    png(paste0(folder_path,'/Precision_estimated_IPCHD.png'), width= 800, height = 600)
    plot.pdm(Theta_icphd_BIC$optimal.model$wi, main = 'Estimated IPCHD')
    dev.off()
    
  }  
  # Save Comparisons
  raw_results = all_comparisons(Theta, Theta_icphd_BIC$optimal.model$wi,  args$r, threshold = threshold)
  return(list(results = raw_results, fit = Theta_icphd_BIC))
}

#HW simulation
simulation_BIC_hw = function(Theta, args,lambda, threshold, folder, save = FALSE){
  set.seed(args$seed)
  #simulate data
  sim_data = sample_from_precision(args$n, args$p, Theta)
  X = sim_data$X
  
  #get correlation and covariance matrix
  CoVR_X = get_cov_cor(X)
  
  #get influence measures
  hw = hw_weights(CoVR_X$Cov)
  hw_inf = hw$influence
  w_hw = hw$w
  
  #get optimal lambda
  
  max_rho = tp.weightedGL(X = X, .W = w_hw, p = args$p, n = args$n)
  lambda = lambda * max_rho
  
  #Fit GLASSO
  #should be the same as the one below
  #Theta_hw_BIC = BICweighted_glasso(mat = CoVR_X$Cov, rho = lambda, .W = w_hw, p = args$p, n = args$n,
  #                                  threshold = 0, maxit = 200,
  #                                  penalize.diagonal = TRUE)
  
  Theta_hw_BIC = BIChwglasso(mat = CoVR_X$Cov, rho = lambda, p = args$p, n = args$n,
                             threshold = 0, maxit = 200,
                             penalize.diagonal = FALSE)
  
  if(save){
    #Save results to analize
    folder_path = paste0(folder, '/HW_BIC')
    
    
    if (!dir.exists(folder_path)) {
      dir.create(folder_path)
    }
    #############################------------ Saving Results ------------------------------------------#
    # Save Influence Measure
    
    png(paste0(folder_path,"/HW_influence.png"), width = 800, height = 600)  # Set file name and dimensions
    plot(hw_inf, main = 'L1 norm of rows', ylab = 'L1', xlab = 'Variable')
    dev.off()  # Save and close the file
    
    #Save weights for a hub and non hub
    png(paste0(folder_path,"/Weight_hub_node_hw.png"), width = 800, height = 600)
    plot(w_hw[1, 1:args$p], main = 'Weights for HWGLASSO', xlab = 'Variable', ylab = 'Weight of of entry 1,j', 
         ylim = c(min(w_hw[args$p, 1:args$p]), quantile(w_hw[args$p, 1:args$p], 0.95)))
    points(w_hw[args$p, 1:args$p], col = 'blue')
    legend("topright", legend = c("Hub node - 1", "Non-hub last node"), col = c("black", "blue"), pch = 19)
    dev.off()
    
    #Save Sanity check for fitting method
    
    png(paste0(folder_path,"/BIC_curve_hw.png"), width = 800, height = 600)
    plot(x = lambda, y = Theta_hw_BIC$BIC, log = 'x', main ='BIC for HWgLASSO', ylab = 'BIC', xlab = 'Lambda')
    dev.off()
    
    #Save matrices visualizations
    
    png(paste0(folder_path,"/Absolute_diff_hw.png"), width = 800, height = 600)
    plot.pdm(abs(Theta - Theta_hw_BIC$optimal.model$wi), main = 'Absolute Difference Estimated hw')
    dev.off()
    
    png(paste0(folder_path,'/Precision_estimated_hw.png'), width= 800, height = 600)
    plot.pdm(Theta_hw_BIC$optimal.model$wi, main = 'Estimated hw')
    dev.off()
    
    p <- plot_ly(z = w_hw, type = "heatmap", colorscale = "Viridis")
    
    # Save the Plotly plot to an HTML file
    html_file <- paste0(folder_path,"/Weight_matrix.html")
    htmlwidgets::saveWidget(p, html_file)
    
    # Use webshot to capture the HTML file and save as PNG
    webshot(html_file, file = paste0(folder_path,"/Weight_matrix.png"))
  }
    # Save Comparisons
  raw_results = all_comparisons(Theta, Theta_hw_BIC$optimal.model$wi, args$r, threshold = threshold)
  return(list(results = raw_results, fit = Theta_hw_BIC))
}


simulation_Lik_ipchd = function(Theta, args, w_func ,lambda, threshold, folder, save = FALSE, name=""){
  set.seed(args$seed)
  #simulate data
  sim_data = sample_from_precision(args$n, args$p, Theta)
  X_train = sim_data$X_train
  X_val = sim_data$X_val
  X = X_train
  
  #get correlation and covariance matrix
  CoVR_X_train = get_cov_cor(X_train)
  CoVR_X_val = get_cov_cor(X_val)
  
  #get influence measures
  icphd_inf = sta_ipchd(X_train, mat_type = 'cor', mat = CoVR_X_train$Corr, var_inds = 1:args$p, overest_type = 'frac')
  w_icphd = w_func(args$p, icphd_inf)
  
  #get optimal lambda
  
  max_rho = tp.weightedGL(X = X, .W = w_icphd, p = args$p, n = args$n)
  lambda = lambda * max_rho
  
  #Fit GLASSO
  Theta_icphd_Lik = Lik_weighted_glasso(mat_train = CoVR_X_train$Cov, mat_val = CoVR_X_val$Cov, rho = lambda, 
                                        .W = w_icphd , p = args$p, n = args$n,
                                        threshold = 0, maxit = 200,
                                        penalize.diagonal = FALSE)
  #Save results to analize
  
  folder_path = paste0(folder, '/IPCHD_Lik_',name)
  
  if(save){
  
    if (!dir.exists(folder_path)) {
      dir.create(folder_path)
    }
    #############################------------ Saving Results ------------------------------------------#
    # Save Influence Measure
    
    png(paste0(folder_path,"/IPCHD_influence.png"), width = 800, height = 600)  # Set file name and dimensions
    plot(icphd_inf, main = 'Influence measure of correlation', ylab = 'ICPHD Measure', xlab = 'Variable')
    dev.off()  # Save and close the file
    
    #Save weights for a hub and non hub
    png(paste0(folder_path,"/Weight_hub_node_IPCHD.png"), width = 800, height = 600)
    plot(w_icphd[1, 1:args$p], main = paste('Weights for IPCHD', name), xlab = 'Variable', ylab = 'Weight of of entry 1,j', 
         ylim = c(min(w_icphd[args$p, 1:args$p]), quantile(w_icphd[args$p, 1:args$p], 0.95)))
    points(w_icphd[args$p, 1:args$p], col = 'blue')
    legend("topright", legend = c("Hub node - 1", "Non-hub last node"), col = c("black", "blue"), pch = 19)
    dev.off()
    
    p <- plot_ly(z = w_icphd, type = "heatmap", colorscale = "Viridis")
    
    # Save the Plotly plot to an HTML file
    html_file <- paste0(folder_path,"/Weight_matrix.html")
    htmlwidgets::saveWidget(p, html_file)
    
    # Use webshot to capture the HTML file and save as PNG
    webshot(html_file, file = paste0(folder_path,"/Weight_matrix.png"))
    
    #Save Sanity check for fitting method
    
    png(paste0(folder_path,"/Likelihood_curve_IPCHD.png"), width = 800, height = 600)
    plot(x = lambda, y = Theta_icphd_Lik$BIC, log = 'x', main ='-Likelihood for IPCHD', ylab = '-Likelihood', xlab = 'Lambda')
    dev.off()
    
    #Save matrices visualizations
    
    png(paste0(folder_path,"/Absolute_diff_IPCHD.png"), width = 800, height = 600)
    plot.pdm(abs(Theta - Theta_icphd_Lik$optimal.model$wi), main = 'Absolute Difference Estimated IPCHD')
    dev.off()
    
    png(paste0(folder_path,'/Precision_estimated_IPCHD.png'), width= 800, height = 600)
    plot.pdm(Theta_icphd_Lik$optimal.model$wi, main = 'Estimated IPCHD')
    dev.off()
  }  
  
  # Save Comparisons
  raw_results = all_comparisons(Theta, Theta_icphd_Lik$optimal.model$wi, args$r, threshold = threshold)
  return(list(results = raw_results, fit = Theta_icphd_Lik))
}


#HW simulation
simulation_Lik_hw = function(Theta, args,lambda, threshold, folder, save = FALSE){
  set.seed(args$seed)
  #simulate data
  sim_data = sample_from_precision(args$n, args$p, Theta)
  X_train = sim_data$X_train
  X_val = sim_data$X_val
  X = X_train
  
  #get correlation and covariance matrix
  CoVR_X_train = get_cov_cor(X_train)
  CoVR_X_val = get_cov_cor(X_val)
  
  #get influence measures
  hw = hw_weights(CoVR_X_train$Cov)
  hw_inf = hw$influence
  w_hw = hw$w
  
  #get optimal lambda
  
  max_rho = tp.weightedGL(X = X, .W = w_hw, p = args$p, n = args$n)
  lambda = lambda * max_rho
  
  #Fit GLASSO
  Theta_hw_Lik = Lik_weighted_glasso(mat_train = CoVR_X_train$Cov, mat_val = CoVR_X_val$Cov, rho = lambda, 
                                     .W = w_hw , p = args$p, n = args$n,
                                     threshold = 0, maxit = 200,
                                     penalize.diagonal = FALSE)
  
  if(save){
    #Save results to analize
    folder_path = paste0(folder, '/HW_Lik')
    
    
    if (!dir.exists(folder_path)) {
      dir.create(folder_path)
    }
    #############################------------ Saving Results ------------------------------------------#
    # Save Influence Measure
    
    png(paste0(folder_path,"/HW_influence.png"), width = 800, height = 600)  # Set file name and dimensions
    plot(hw_inf, main = 'L1 norm of rows', ylab = 'L1', xlab = 'Variable')
    dev.off()  # Save and close the file
    
    #Save weights for a hub and non hub
    png(paste0(folder_path,"/Weight_hub_node_hw.png"), width = 800, height = 600)
    plot(w_hw[1, 1:args$p], main = 'Weights for HWGLASSO', xlab = 'Variable', ylab = 'Weight of of entry 1,j', 
         ylim = c(min(w_hw[args$p, 1:args$p]), quantile(w_hw[args$p, 1:args$p], 0.95)))
    points(w_hw[args$p, 1:args$p], col = 'blue')
    legend("topright", legend = c("Hub node - 1", "Non-hub last node"), col = c("black", "blue"), pch = 19)
    dev.off()
    
    p <- plot_ly(z = w_hw, type = "heatmap", colorscale = "Viridis")
    
    # Save the Plotly plot to an HTML file
    html_file <- paste0(folder_path,"/Weight_matrix.html")
    htmlwidgets::saveWidget(p, html_file)
    
    # Use webshot to capture the HTML file and save as PNG
    webshot(html_file, file = paste0(folder_path,"/Weight_matrix.png"))
    #Save Sanity check for fitting method
    
    png(paste0(folder_path,"/Likelihood_curve_hw.png"), width = 800, height = 600)
    plot(x = lambda, y = Theta_hw_Lik$BIC, log = 'x', main ='BIC for HWgLASSO', ylab = 'BIC', xlab = 'Lambda')
    dev.off()
    
    #Save matrices visualizations
    
    png(paste0(folder_path,"/Absolute_diff_hw.png"), width = 800, height = 600)
    plot.pdm(abs(Theta - Theta_hw_Lik$optimal.model$wi), main = 'Absolute Difference Estimated hw')
    dev.off()
    
    png(paste0(folder_path,'/Precision_estimated_hw.png'), width= 800, height = 600)
    plot.pdm(Theta_hw_Lik$optimal.model$wi, main = 'Estimated hw')
    dev.off()
  }  
  
  # Save Comparisons
  raw_results = all_comparisons(Theta, Theta_hw_Lik$optimal.model$wi, args$r, threshold = threshold)
  return(list(results = raw_results, fit = Theta_hw_Lik))
}

#adaptive Lasso simulation
simulation_BIC_adaGL = function(Theta, args,lambda, threshold, folder, save = FALSE){
  set.seed(args$seed)
  #simulate data
  sim_data = sample_from_precision(args$n, args$p, Theta)
  X = sim_data$X
  
  #get correlation and covariance matrix
  CoVR_X = get_cov_cor(X)
  
  #get influence measures
  W = adaGL_weights(CoVR_X$Cov)
  
  #get optimal lambda
  
  max_rho = tp.weightedGL(X = X, .W = W, p = args$p, n = args$n)
  lambda = lambda * max_rho
  
  #Fit GLASSO
  #should be the same as the one below
  #Theta_hw_BIC = BICweighted_glasso(mat = CoVR_X$Cov, rho = lambda, .W = w_hw, p = args$p, n = args$n,
  #                                  threshold = 0, maxit = 200,
  #                                  penalize.diagonal = TRUE)
  
  Theta_ada_BIC = BICweighted_glasso(mat = CoVR_X$Cov, rho = lambda, .W = W, p = args$p, n = args$n,
                                    threshold = 0, maxit = 200,
                                    penalize.diagonal = FALSE)
  
  if(save){
  
    #Save results to analize
    folder_path = paste0(folder, '/Ada_BIC')
    
    
    if (!dir.exists(folder_path)) {
      dir.create(folder_path)
    }
    #############################------------ Saving Results ------------------------------------------#
    # Save Influence Measure
  
    #Save weights for a hub and non hub
    png(paste0(folder_path,"/Weight_hub_node_ada.png"), width = 800, height = 600)
    plot(W[1, 1:args$p], main = 'Weights for ada', xlab = 'Variable', ylab = 'Weight of of entry 1,j', 
         ylim = c(min(W[args$p, 1:args$p]), quantile(W[args$p, 1:args$p], 0.95)))
    points(W[args$p, 1:args$p], col = 'blue')
    legend("topright", legend = c("Hub node - 1", "Non-hub last node"), col = c("black", "blue"), pch = 19)
    dev.off()
    
    p <- plot_ly(z = W, type = "heatmap", colorscale = "Viridis")
    
    # Save the Plotly plot to an HTML file
    html_file <- paste0(folder_path,"/Weight_matrix.html")
    htmlwidgets::saveWidget(p, html_file)
    
    # Use webshot to capture the HTML file and save as PNG
    webshot(html_file, file = paste0(folder_path,"/Weight_matrix.png"))
    
    
    #Save Sanity check for fitting method
    
    png(paste0(folder_path,"/BIC_curve_ada.png"), width = 800, height = 600)
    plot(x = lambda, y = Theta_ada_BIC$BIC, log = 'x', main ='BIC for Ada', ylab = 'BIC', xlab = 'Lambda')
    dev.off()
    
    #Save matrices visualizations
    
    png(paste0(folder_path,"/Absolute_diff_ada.png"), width = 800, height = 600)
    plot.pdm(abs(Theta - Theta_ada_BIC$optimal.model$wi), main = 'Absolute Difference Estimated ada BIC')
    dev.off()
    
    png(paste0(folder_path,'/Precision_estimated_ada.png'), width= 800, height = 600)
    plot.pdm(Theta_ada_BIC$optimal.model$wi, main = 'Estimated ada')
    dev.off()
  }  
  
  # Save Comparisons
  raw_results = all_comparisons(Theta, Theta_ada_BIC$optimal.model$wi, args$r, threshold = threshold)
  return(list(results = raw_results, fit = Theta_ada_BIC))
}

#GLasso simulation
simulation_BIC_GL = function(Theta, args,lambda, threshold, folder, save = FALSE){
  set.seed(args$seed)
  #simulate data
  sim_data = sample_from_precision(args$n, args$p, Theta)
  X = sim_data$X
  
  #get correlation and covariance matrix
  CoVR_X = get_cov_cor(X)
  
  #get influence measures
  W = (matrix(1, args$p, args$p) - diag(1, args$p, args$p)) * 0.5
  
  #get optimal lambda
  
  max_rho = tp.weightedGL(X = X, .W = W, p = args$p, n = args$n)
  lambda = lambda * max_rho
  
  #Fit GLASSO
  #should be the same as the one below
  #Theta_hw_BIC = BICweighted_glasso(mat = CoVR_X$Cov, rho = lambda, .W = w_hw, p = args$p, n = args$n,
  #                                  threshold = 0, maxit = 200,
  #                                  penalize.diagonal = TRUE)
  
  Theta_GL_BIC = BICweighted_glasso(mat = CoVR_X$Cov, rho = lambda, .W = W, p = args$p, n = args$n,
                                     threshold = 0, maxit = 200,
                                     penalize.diagonal = FALSE)
  
  
  if(save){
  
    #Save results to analize
    folder_path = paste0(folder, '/GL_BIC')
    
    
    if (!dir.exists(folder_path)) {
      dir.create(folder_path)
    }
    #############################------------ Saving Results ------------------------------------------#
    # Save Influence Measure
    
    #Save weights for a hub and non hub
    png(paste0(folder_path,"/Weight_hub_node_GL.png"), width = 800, height = 600)
    plot(W[1, 1:args$p], main = 'Weights for GL', xlab = 'Variable', ylab = 'Weight of of entry 1,j', 
         ylim = c(min(W[args$p, 1:args$p]), quantile(W[args$p, 1:args$p], 0.95)))
    points(W[args$p, 1:args$p], col = 'blue')
    legend("topright", legend = c("Hub node - 1", "Non-hub last node"), col = c("black", "blue"), pch = 19)
    dev.off()
    
    p <- plot_ly(z = W, type = "heatmap", colorscale = "Viridis")
    
    # Save the Plotly plot to an HTML file
    html_file <- paste0(folder_path,"/Weight_matrix.html")
    htmlwidgets::saveWidget(p, html_file)
    
    # Use webshot to capture the HTML file and save as PNG
    webshot(html_file, file = paste0(folder_path,"/Weight_matrix.png"))
    
    #Save Sanity check for fitting method
    
    png(paste0(folder_path,"/BIC_curve_GL.png"), width = 800, height = 600)
    plot(x = lambda, y = Theta_GL_BIC$BIC, log = 'x', main ='BIC for GL', ylab = 'BIC', xlab = 'Lambda')
    dev.off()
    
    #Save matrices visualizations
    
    png(paste0(folder_path,"/Absolute_diff_GL.png"), width = 800, height = 600)
    plot.pdm(abs(Theta - Theta_GL_BIC$optimal.model$wi), main = 'Absolute Difference Estimated GL BIC')
    dev.off()
    
    png(paste0(folder_path,'/Precision_estimated_gl.png'), width= 800, height = 600)
    plot.pdm(Theta_GL_BIC$optimal.model$wi, main = 'Estimated GL')
    dev.off()
  }  
  
  # Save Comparisons
  raw_results = all_comparisons(Theta, Theta_GL_BIC$optimal.model$wi, args$r, threshold = threshold)
  return(list(results = raw_results, fit = Theta_GL_BIC))
}


#HW simulation
simulation_Lik_adaGL = function(Theta, args,lambda, threshold, folder, save=FALSE){
  set.seed(args$seed)
  #simulate data
  sim_data = sample_from_precision(args$n, args$p, Theta)
  X_train = sim_data$X_train
  X_val = sim_data$X_val
  X = X_train
  
  #get correlation and covariance matrix
  CoVR_X_train = get_cov_cor(X_train)
  CoVR_X_val = get_cov_cor(X_val)
  
  #get influence measures
  W = adaGL_weights(CoVR_X_train$Cov)
  
  #get optimal lambda
  
  max_rho = tp.weightedGL(X = X, .W = W, p = args$p, n = args$n)
  lambda = lambda * max_rho
  
  #Fit GLASSO
  Theta_ada_Lik = Lik_weighted_glasso(mat_train = CoVR_X_train$Cov, mat_val = CoVR_X_val$Cov, rho = lambda, 
                                     .W = W , p = args$p, n = args$n,
                                     threshold = 0, maxit = 200,
                                     penalize.diagonal = FALSE)
  
  if(save){
    #Save results to analize
    folder_path = paste0(folder, '/Ada_Lik')
    
    
    if (!dir.exists(folder_path)) {
      dir.create(folder_path)
    }
    #############################------------ Saving Results ------------------------------------------#
    
    #Save weights for a hub and non hub
    png(paste0(folder_path,"/Weight_hub_node_ada.png"), width = 800, height = 600)
    plot(W[1, 1:args$p], main = 'Weights for Ada', xlab = 'Variable', ylab = 'Weight of of entry 1,j', 
         ylim = c(min(W[args$p, 1:args$p]), quantile(W[args$p, 1:args$p], 0.95)))
    points(W[args$p, 1:args$p], col = 'blue')
    legend("topright", legend = c("Hub node - 1", "Non-hub last node"), col = c("black", "blue"), pch = 19)
    dev.off()
    
    p <- plot_ly(z = W, type = "heatmap", colorscale = "Viridis")
    
    # Save the Plotly plot to an HTML file
    html_file <- paste0(folder_path,"/Weight_matrix.html")
    htmlwidgets::saveWidget(p, html_file)
    
    # Use webshot to capture the HTML file and save as PNG
    webshot(html_file, file = paste0(folder_path,"/Weight_matrix.png"))
    
    #Save Sanity check for fitting method
    
    png(paste0(folder_path,"/Likelihood_curve_ada.png"), width = 800, height = 600)
    plot(x = lambda, y = Theta_ada_Lik$BIC, log = 'x', main ='BIC for Ada', ylab = 'BIC', xlab = 'Lambda')
    dev.off()
    
    #Save matrices visualizations
    
    png(paste0(folder_path,"/Absolute_diff_ada.png"), width = 800, height = 600)
    plot.pdm(abs(Theta - Theta_ada_Lik$optimal.model$wi), main = 'Absolute Difference Estimated ada')
    dev.off()
    
    png(paste0(folder_path,'/Precision_estimated_ada.png'), width= 800, height = 600)
    plot.pdm(Theta_ada_Lik$optimal.model$wi, main = 'Estimated Ada')
    dev.off()
  }  
  
  # Save Comparisons
  raw_results = all_comparisons(Theta, Theta_ada_Lik$optimal.model$wi, args$r, threshold = threshold)
  return(list(results = raw_results, fit = Theta_ada_Lik))
}

#HW simulation
simulation_Lik_GL = function(Theta, args,lambda, threshold, folder, save = FALSE){
  set.seed(args$seed)
  #simulate data
  sim_data = sample_from_precision(args$n, args$p, Theta)
  X_train = sim_data$X_train
  X_val = sim_data$X_val
  X = X_train
  
  #get correlation and covariance matrix
  CoVR_X_train = get_cov_cor(X_train)
  CoVR_X_val = get_cov_cor(X_val)
  
  #get influence measures
  W = (matrix(1, args$p, args$p) - diag(1, args$p, args$p)) * 0.5
  
  #get optimal lambda
  
  max_rho = tp.weightedGL(X = X, .W = W, p = args$p, n = args$n)
  lambda = lambda * max_rho
  
  #Fit GLASSO
  Theta_GL_Lik = Lik_weighted_glasso(mat_train = CoVR_X_train$Cov, mat_val = CoVR_X_val$Cov, rho = lambda, 
                                      .W = W , p = args$p, n = args$n,
                                      threshold = 0, maxit = 200,
                                      penalize.diagonal = FALSE)
  
  if(save){
    #Save results to analize
    folder_path = paste0(folder, '/GL_Lik')
    
    
    if (!dir.exists(folder_path)) {
      dir.create(folder_path)
    }
    #############################------------ Saving Results ------------------------------------------#
    
    #Save weights for a hub and non hub
    png(paste0(folder_path,"/Weight_hub_node_GL.png"), width = 800, height = 600)
    plot(W[1, 1:args$p], main = 'Weights for GL', xlab = 'Variable', ylab = 'Weight of of entry 1,j', 
         ylim = c(min(W[args$p, 1:args$p]), quantile(W[args$p, 1:args$p], 0.95)))
    points(W[args$p, 1:args$p], col = 'blue')
    legend("topright", legend = c("Hub node - 1", "Non-hub last node"), col = c("black", "blue"), pch = 19)
    dev.off()
    
    p <- plot_ly(z = W, type = "heatmap", colorscale = "Viridis")
    
    # Save the Plotly plot to an HTML file
    html_file <- paste0(folder_path,"/Weight_matrix.html")
    htmlwidgets::saveWidget(p, html_file)
    
    # Use webshot to capture the HTML file and save as PNG
    webshot(html_file, file = paste0(folder_path,"/Weight_matrix.png"))
    
    #Save Sanity check for fitting method
    
    png(paste0(folder_path,"/Likelihood_curve_GL.png"), width = 800, height = 600)
    plot(x = lambda, y = Theta_GL_Lik$BIC, log = 'x', main ='BIC for GL', ylab = 'BIC', xlab = 'Lambda')
    dev.off()
    
    #Save matrices visualizations
    
    png(paste0(folder_path,"/Absolute_diff_GL.png"), width = 800, height = 600)
    plot.pdm(abs(Theta - Theta_GL_Lik$optimal.model$wi), main = 'Absolute Difference Estimated GL')
    dev.off()
    
    png(paste0(folder_path,'/Precision_estimated_GL.png'), width= 800, height = 600)
    plot.pdm(Theta_GL_Lik$optimal.model$wi, main = 'Estimated GL')
    dev.off()
  }  
  
  # Save Comparisons
  raw_results = all_comparisons(Theta, Theta_GL_Lik$optimal.model$wi, args$r, threshold = threshold)
  return(list(results = raw_results, fit = Theta_GL_Lik))
}


#GLasso simulation
simulation_BIC_oracle = function(Theta, args,lambda, threshold, folder, save = FALSE){
  set.seed(args$seed)
  #simulate data
  sim_data = sample_from_precision(args$n, args$p, Theta)
  X = sim_data$X
  
  #get correlation and covariance matrix
  CoVR_X = get_cov_cor(X)
  
  #get influence measures
  W = matrix(1, args$p, args$p)
  W[1:args$r, ] = 0
  W[, 1:args$r] = 0
  diag(W) = 0
  
  #get optimal lambda
  
  max_rho = tp.weightedGL(X = X, .W = W, p = args$p, n = args$n)
  lambda = lambda * max_rho
  
  #Fit GLASSO
  #should be the same as the one below
  #Theta_hw_BIC = BICweighted_glasso(mat = CoVR_X$Cov, rho = lambda, .W = w_hw, p = args$p, n = args$n,
  #                                  threshold = 0, maxit = 200,
  #                                  penalize.diagonal = TRUE)
  
  Theta_GL_BIC = BICweighted_glasso(mat = CoVR_X$Cov, rho = lambda, .W = W, p = args$p, n = args$n,
                                    threshold = 0, maxit = 200,
                                    penalize.diagonal = FALSE)
  
  if(save){
    #Save results to analize
    folder_path = paste0(folder, '/Oracle_BIC')
    
    
    if (!dir.exists(folder_path)) {
      dir.create(folder_path)
    }
    #############################------------ Saving Results ------------------------------------------#
    # Save Influence Measure
    
    #Save weights for a hub and non hub
    png(paste0(folder_path,"/Weight_hub_node_oracle.png"), width = 800, height = 600)
    plot(W[1, 1:args$p], main = 'Weights for oracle', xlab = 'Variable', ylab = 'Weight of of entry 1,j', 
         ylim = c(min(W[args$p, 1:args$p]), quantile(W[args$p, 1:args$p], 0.95)))
    points(W[args$p, 1:args$p], col = 'blue')
    legend("topright", legend = c("Hub node - 1", "Non-hub last node"), col = c("black", "blue"), pch = 19)
    dev.off()
    
    p <- plot_ly(z = W, type = "heatmap", colorscale = "Viridis")
    
    # Save the Plotly plot to an HTML file
    html_file <- paste0(folder_path,"/Weight_matrix.html")
    htmlwidgets::saveWidget(p, html_file)
    
    # Use webshot to capture the HTML file and save as PNG
    webshot(html_file, file = paste0(folder_path,"/Weight_matrix.png"))
    
    #Save Sanity check for fitting method
    
    png(paste0(folder_path,"/BIC_curve_oracle.png"), width = 800, height = 600)
    plot(x = lambda, y = Theta_GL_BIC$BIC, log = 'x', main ='BIC for GL', ylab = 'BIC', xlab = 'Lambda')
    dev.off()
    
    #Save matrices visualizations
    
    png(paste0(folder_path,"/Absolute_diff_oracle.png"), width = 800, height = 600)
    plot.pdm(abs(Theta - Theta_GL_BIC$optimal.model$wi), main = 'Absolute Difference Estimated oracle BIC')
    dev.off()
    
    png(paste0(folder_path,'/Precision_estimated_oracle.png'), width= 800, height = 600)
    plot.pdm(Theta_GL_BIC$optimal.model$wi, main = 'Estimated oracle')
    dev.off()
  }  
  
  # Save Comparisons
  raw_results = all_comparisons(Theta, Theta_GL_BIC$optimal.model$wi, args$r, threshold = threshold)
  return(list(results = raw_results, fit = Theta_GL_BIC))
}


simulation_Lik_oracle = function(Theta, args,lambda, threshold, folder, save = FALSE){
  set.seed(args$seed)
  #simulate data
  sim_data = sample_from_precision(args$n, args$p, Theta)
  X_train = sim_data$X_train
  X_val = sim_data$X_val
  X = X_train
  
  #get correlation and covariance matrix
  CoVR_X_train = get_cov_cor(X_train)
  CoVR_X_val = get_cov_cor(X_val)
  
  #get influence measures
  #get influence measures
  W = matrix(1, args$p, args$p)
  W[1:args$r, ] = 0
  W[, 1:args$r] = 0
  diag(W) = 0
  #get optimal lambda
  
  max_rho = tp.weightedGL(X = X, .W = W, p = args$p, n = args$n)
  lambda = lambda * max_rho
  
  #Fit GLASSO
  Theta_GL_Lik = Lik_weighted_glasso(mat_train = CoVR_X_train$Cov, mat_val = CoVR_X_val$Cov, rho = lambda, 
                                     .W = W , p = args$p, n = args$n,
                                     threshold = 0, maxit = 200,
                                     penalize.diagonal = FALSE)
  
  if(save){
    #Save results to analize
    folder_path = paste0(folder, '/oracle_Lik')
    
    
    if (!dir.exists(folder_path)) {
      dir.create(folder_path)
    }
    #############################------------ Saving Results ------------------------------------------#
    
    #Save weights for a hub and non hub
    png(paste0(folder_path,"/Weight_hub_node_oracle.png"), width = 800, height = 600)
    plot(W[1, 1:args$p], main = 'Weights for GL', xlab = 'Variable', ylab = 'Weight of of entry 1,j', 
         ylim = c(min(W[args$p, 1:args$p]), quantile(W[args$p, 1:args$p], 0.95)))
    points(W[args$p, 1:args$p], col = 'blue')
    legend("topright", legend = c("Hub node - 1", "Non-hub last node"), col = c("black", "blue"), pch = 19)
    dev.off()
    
    p <- plot_ly(z = W, type = "heatmap", colorscale = "Viridis")
    
    # Save the Plotly plot to an HTML file
    html_file <- paste0(folder_path,"/Weight_matrix.html")
    htmlwidgets::saveWidget(p, html_file)
    
    # Use webshot to capture the HTML file and save as PNG
    webshot(html_file, file = paste0(folder_path,"/Weight_matrix.png"))
    
    
    #Save Sanity check for fitting method
    
    png(paste0(folder_path,"/Likelihood_curve_oracle.png"), width = 800, height = 600)
    plot(x = lambda, y = Theta_GL_Lik$BIC, log = 'x', main ='BIC for Ada', ylab = 'BIC', xlab = 'Lambda')
    dev.off()
    
    #Save matrices visualizations
    
    png(paste0(folder_path,"/Absolute_diff_oracle.png"), width = 800, height = 600)
    plot.pdm(abs(Theta - Theta_GL_Lik$optimal.model$wi), main = 'Absolute Difference Estimated oracle')
    dev.off()
    
    png(paste0(folder_path,'/Precision_estimated_oracle.png'), width= 800, height = 600)
    plot.pdm(Theta_GL_Lik$optimal.model$wi, main = 'Estimated oracle')
    dev.off()
  }  
  
  # Save Comparisons
  raw_results = all_comparisons(Theta, Theta_GL_Lik$optimal.model$wi, args$r, threshold = threshold)
  return(list(results = raw_results, fit = Theta_GL_Lik))
}


#GLasso simulation
simulation_BIC_2step = function(Theta, args,lambda, threshold, folder, save=FALSE){
  set.seed(args$seed)
  #simulate data
  sim_data = sample_from_precision(args$n, args$p, Theta)
  X = sim_data$X
  
  #get correlation and covariance matrix
  CoVR_X = get_cov_cor(X)
  
  #get influence measures
  W = matrix(1, args$p, args$p)
  diag(W) = 0
  
  #get optimal lambda
  
  max_rho = tp.weightedGL(X = X, .W = W, p = args$p, n = args$n)
  lambda = lambda * max_rho * 2
  
  #Fit GLASSO
  #should be the same as the one below
  #Theta_hw_BIC = BICweighted_glasso(mat = CoVR_X$Cov, rho = lambda, .W = w_hw, p = args$p, n = args$n,
  #                                  threshold = 0, maxit = 200,
  #                                  penalize.diagonal = TRUE)
  
  Theta_GL_BIC = BIC_2step_glasso(mat = CoVR_X$Cov, rho = lambda, p = args$p, n = args$n,
                                  threshold = 0, maxit = 200,
                                  penalize.diagonal = FALSE)
  
  if(save){
    #Save results to analize
    folder_path = paste0(folder, '/2step_BIC')
    
    
    if (!dir.exists(folder_path)) {
      dir.create(folder_path)
    }
    #############################------------ Saving Results ------------------------------------------#
    # Save Influence Measure
  
    #Save Sanity check for fitting method
    
    png(paste0(folder_path,"/BIC_curve_2step.png"), width = 800, height = 600)
    persp(log(lambda), log(lambda), Theta_GL_BIC$BIC, 
                 xlab = "log Lambda 1", ylab = "log Lambda 2", zlab = "BIC",
                 theta = 30, phi = 20, col = "lightblue", main = "3D BIC Plot")
    dev.off()
    
    png(paste0(folder_path,"/BIC_contour.png"), width = 800, height = 600)
    filled.contour(log(lambda), log(lambda), Theta_GL_BIC$BIC,
                                   xlab = "log lambda 1 ", ylab = "log lambda 2",
                                   main = "Contour Plot of BIC", color.palette = terrain.colors)
    dev.off()
    #Save matrices visualizations
    
    png(paste0(folder_path,"/Absolute_diff.png"), width = 800, height = 600)
    plot.pdm(abs(Theta - Theta_GL_BIC$optimal.model$wi), main = 'Absolute Difference Estimated 2step')
    dev.off()
    
    png(paste0(folder_path,'/Precision_estimated.png'), width= 800, height = 600)
    plot.pdm(Theta_GL_BIC$optimal.model$wi, main = 'Estimated 2step')
    dev.off()
  }  
  
  # Save Comparisons
  raw_results = all_comparisons(Theta, Theta_GL_BIC$optimal.model$wi, args$r, threshold = threshold)
  return(list(results = raw_results, fit = Theta_GL_BIC))
}

# simulation_BIC_ipchd_frob = function(Theta, args, w_func ,lambda, threshold, folder, name){
#   cat("BIC")
#   #simulate data
#   set.seed(args$seed)
#   sim_data = sample_from_precision(args$n, args$p, Theta)
#   X = sim_data$X
#   
#   #get correlation and covariance matrix
#   CoVR_X = get_cov_cor(X)
#   
#   #get influence measures
#   icphd_inf = sta_ipchd(X, mat_type = 'cor', mat = CoVR_X$Corr, var_inds = 1:args$p, overest_type = 'frac')
#   w_icphd = w_func(args$p, icphd_inf)
#   
#   #get optimal lambda
#   
#   max_rho = tp.weightedGL(X = X, .W = w_icphd, p = args$p, n = args$n)
#   lambda = lambda * max_rho
#   
#   #Fit GLASSO
#   Theta_icphd_BIC = BICweighted_glasso_frob(Theta = Theta, mat = CoVR_X$Cov, rho = lambda, .W = w_icphd, p = args$p, n = args$n,
#                                        threshold = 0, maxit = 200,
#                                        penalize.diagonal = FALSE)
#   #Save results to analize
#   
#   folder_path = paste0(folder, '/IPCHD_frob_',name)
#   
#   
#   
#   if (!dir.exists(folder_path)) {
#     dir.create(folder_path)
#   }
#   #############################------------ Saving Results ------------------------------------------#
#   # Save Influence Measure
#   
#   png(paste0(folder_path,"/IPCHD_influence.png"), width = 800, height = 600)  # Set file name and dimensions
#   plot(icphd_inf, main = 'Influence measure of correlation', ylab = 'ICPHD Measure', xlab = 'Variable')
#   dev.off()  # Save and close the file
#   
#   #Save weights for a hub and non hub
#   png(paste0(folder_path,"/Weight_hub_node_IPCHD.png"), width = 800, height = 600)
#   plot(w_icphd[1, 1:args$p], main = paste('Weights for IPCHD', name), xlab = 'Variable', ylab = 'Weight of of entry 1,j', 
#        ylim = c(min(w_icphd[args$p, 1:args$p]), quantile(w_icphd[args$p, 1:args$p], 0.95)))
#   points(w_icphd[args$p, 1:args$p], col = 'blue')
#   legend("topright", legend = c("Hub node - 1", "Non-hub last node"), col = c("black", "blue"), pch = 19)
#   dev.off()
#   
#   #Save Sanity check for fitting method
#   
#   png(paste0(folder_path,"/frob_curve_IPCHD.png"), width = 800, height = 600)
#   plot(x = lambda, y = Theta_icphd_BIC$BIC, log = 'x', main ='BIC for IPCHD', ylab = 'BIC', xlab = 'Lambda')
#   dev.off()
#   
#   #Save matrices visualizations
#   
#   png(paste0(folder_path,"/Absolute_diff_IPCHD.png"), width = 800, height = 600)
#   plot.pdm(abs(Theta - Theta_icphd_BIC$optimal.model$wi), main = 'Absolute Difference Estimated IPCHD')
#   dev.off()
#   
#   png(paste0(folder_path,'/Precision_estimated_IPCHD.png'), width= 800, height = 600)
#   plot.pdm(Theta_icphd_BIC$optimal.model$wi, main = 'Estimated IPCHD')
#   dev.off()
#   
#   
#   # Save Comparisons
#   raw_results = all_comparisons(Theta, Theta_icphd_BIC$optimal.model$wi, threshold = threshold)
#   return(list(results = raw_results, fit = Theta_icphd_BIC))
# }
############################################################# NEW CHECK

simulation_BIC_sig_adj = function(Theta, args, lambda, threshold, folder, save=FALSE, name =""){
  cat("BIC")
  #simulate data
  set.seed(args$seed)
  sim_data = sample_from_precision(args$n, args$p, Theta)
  X = sim_data$X
  
  #get correlation and covariance matrix
  CoVR_X = get_cov_cor(X)
  
  #get influence measures
  icphd_inf = sta_ipchd(X, mat_type = 'cor', mat = CoVR_X$Corr, var_inds = 1:args$p, overest_type = 'frac')
  s = 1 - icphd_inf
  measure = get_sig_tr(s)
  w_icphd = sig_prop(args$p, measure)

  #get optimal lambda
  
  max_rho = 3
  lambda = lambda * max_rho
  
  #Fit GLASSO
  Theta_icphd_BIC = BICweighted_glasso(mat = CoVR_X$Cov, rho = lambda, .W = w_icphd, p = args$p, n = args$n,
                                       threshold = 0, maxit = 200,
                                       penalize.diagonal = FALSE)
  #Save results to analize
  
  if(save){
  
    folder_path = paste0(folder, '/SIG_ADJ_BIC_',name)
    
    
    
    if (!dir.exists(folder_path)) {
      dir.create(folder_path)
    }
    #############################------------ Saving Results ------------------------------------------#
    # Save Influence Measure
    
    png(paste0(folder_path,"/SIG_ADJ_influence.png"), width = 800, height = 600)  # Set file name and dimensions
    plot(measure, main = 'Influence adjusted measure of correlation', ylab = 'adjusted ICPHD Measure', xlab = 'Variable')
    dev.off()  # Save and close the file
    
    png(paste0(folder_path,"/IPCHD_influence.png"), width = 800, height = 600)  # Set file name and dimensions
    plot(icphd_inf, main = 'Influence measure of correlation', ylab = 'ICPHD Measure', xlab = 'Variable')
    dev.off()  # Save and close the file
    
    png(paste0(folder_path,"/IPCHD_adj_curv influence.png"), width = 800, height = 600)  # Set file name and dimensions
    plot(s, measure, main = 'adjusted curve', ylab = 'ICPHD Measure', xlab = 'Variable')
    dev.off()  # Save and close the file
    
    #Save weights for a hub and non hub
    png(paste0(folder_path,"/Weight_hub_node_IPCHD.png"), width = 800, height = 600)
    plot(w_icphd[1, 1:args$p], main = paste('Weights for IPCHD', name), xlab = 'Variable', ylab = 'Weight of of entry 1,j', 
         ylim = c(min(w_icphd[args$p, 1:args$p]), quantile(w_icphd[args$p, 1:args$p], 0.95)))
    points(w_icphd[args$p, 1:args$p], col = 'blue')
    legend("topright", legend = c("Hub node - 1", "Non-hub last node"), col = c("black", "blue"), pch = 19)
    dev.off()
    
    p <- plot_ly(z = w_icphd, type = "heatmap", colorscale = "Viridis")
    
    # Save the Plotly plot to an HTML file
    html_file <- paste0(folder_path,"/Weight_matrix.html")
    htmlwidgets::saveWidget(p, html_file)
    
    # Use webshot to capture the HTML file and save as PNG
    webshot(html_file, file = paste0(folder_path,"/Weight_matrix.png"))
    
  
    
    #Save Sanity check for fitting method
    
    png(paste0(folder_path,"/BIC_curve_IPCHD.png"), width = 800, height = 600)
    plot(x = lambda, y = Theta_icphd_BIC$BIC, log = 'x', main ='BIC for IPCHD', ylab = 'BIC', xlab = 'Lambda')
    dev.off()
    
    #Save matrices visualizations
    
    png(paste0(folder_path,"/Absolute_diff_IPCHD.png"), width = 800, height = 600)
    plot.pdm(abs(Theta - Theta_icphd_BIC$optimal.model$wi), main = 'Absolute Difference Estimated IPCHD')
    dev.off()
    
    png(paste0(folder_path,'/Precision_estimated_IPCHD.png'), width= 800, height = 600)
    plot.pdm(Theta_icphd_BIC$optimal.model$wi, main = 'Estimated IPCHD')
    dev.off()
  }  
  
  # Save Comparisons
  raw_results = all_comparisons(Theta, Theta_icphd_BIC$optimal.model$wi, args$r, threshold = threshold)
  return(list(results = raw_results, fit = Theta_icphd_BIC))
}

simulation_BIC_lin_adj = function(Theta, args, lambda, threshold, folder, save=FALSE, name=""){
  cat("BIC")
  #simulate data
  set.seed(args$seed)
  sim_data = sample_from_precision(args$n, args$p, Theta)
  X = sim_data$X
  
  #get correlation and covariance matrix
  CoVR_X = get_cov_cor(X)
  
  #get influence measures
  icphd_inf = sta_ipchd(X, mat_type = 'cor', mat = CoVR_X$Corr, var_inds = 1:args$p, overest_type = 'frac')
  s = 1 - icphd_inf
  measure = get_linear_tr(s)
  w_icphd = sig_prop(args$p, measure)
  
  #get optimal lambda
  
  max_rho = tp.weightedGL(X = X, .W = w_icphd, p = args$p, n = args$n)
  lambda = lambda * max_rho
  
  #Fit GLASSO
  Theta_icphd_BIC = BICweighted_glasso(mat = CoVR_X$Cov, rho = lambda, .W = w_icphd, p = args$p, n = args$n,
                                       threshold = 0, maxit = 200,
                                       penalize.diagonal = FALSE)
  
  if(save){
    
    #Save results to analize
    
    folder_path = paste0(folder, '/LINEAR_ADJ_BIC_',name)
    
    
    
    if (!dir.exists(folder_path)) {
      dir.create(folder_path)
    }
    #############################------------ Saving Results ------------------------------------------#
    # Save Influence Measure
    
    png(paste0(folder_path,"/LINEAR_ADJ_influence.png"), width = 800, height = 600)  # Set file name and dimensions
    plot(measure, main = 'Influence adjusted measure of correlation', ylab = 'adjusted ICPHD Measure', xlab = 'Variable')
    dev.off()  # Save and close the file
    
    png(paste0(folder_path,"/IPCHD_adj_curv influence.png"), width = 800, height = 600)  # Set file name and dimensions
    plot(s, measure, main = 'adjusted curve', ylab = 'ICPHD Measure', xlab = 'Variable')
    dev.off()  # Save and close the file
    
    
    png(paste0(folder_path,"/IPCHD_influence.png"), width = 800, height = 600)  # Set file name and dimensions
    plot(icphd_inf, main = 'Influence measure of correlation', ylab = 'ICPHD Measure', xlab = 'Variable')
    dev.off()  # Save and close the file
    
    #Save weights for a hub and non hub
    png(paste0(folder_path,"/Weight_hub_node_IPCHD.png"), width = 800, height = 600)
    plot(w_icphd[1, 1:args$p], main = paste('Weights for IPCHD', name), xlab = 'Variable', ylab = 'Weight of of entry 1,j', 
         ylim = c(min(w_icphd[args$p, 1:args$p]), quantile(w_icphd[args$p, 1:args$p], 0.95)))
    points(w_icphd[args$p, 1:args$p], col = 'blue')
    legend("topright", legend = c("Hub node - 1", "Non-hub last node"), col = c("black", "blue"), pch = 19)
    dev.off()
    
    p <- plot_ly(z = w_icphd, type = "heatmap", colorscale = "Viridis")
    
    # Save the Plotly plot to an HTML file
    html_file <- paste0(folder_path,"/Weight_matrix.html")
    htmlwidgets::saveWidget(p, html_file)
    
    # Use webshot to capture the HTML file and save as PNG
    webshot(html_file, file = paste0(folder_path,"/Weight_matrix.png"))
    
    #Save Sanity check for fitting method
    
    png(paste0(folder_path,"/BIC_curve_IPCHD.png"), width = 800, height = 600)
    plot(x = lambda, y = Theta_icphd_BIC$BIC, log = 'x', main ='BIC for IPCHD', ylab = 'BIC', xlab = 'Lambda')
    dev.off()
    
    #Save matrices visualizations
    
    png(paste0(folder_path,"/Absolute_diff_IPCHD.png"), width = 800, height = 600)
    plot.pdm(abs(Theta - Theta_icphd_BIC$optimal.model$wi), main = 'Absolute Difference Estimated IPCHD')
    dev.off()
    
    png(paste0(folder_path,'/Precision_estimated_IPCHD.png'), width= 800, height = 600)
    plot.pdm(Theta_icphd_BIC$optimal.model$wi, main = 'Estimated IPCHD')
    dev.off()
  }  
  
  # Save Comparisons
  raw_results = all_comparisons(Theta, Theta_icphd_BIC$optimal.model$wi, args$r, threshold = threshold)
  return(list(results = raw_results, fit = Theta_icphd_BIC))
}


simulation_BIC_lin_oracle = function(Theta, args, lambda, threshold, folder, name = "", save=FALSE){
  cat("BIC")
  #simulate data
  set.seed(args$seed)
  sim_data = sample_from_precision(args$n, args$p, Theta)
  X = sim_data$X
  
  #get correlation and covariance matrix
  CoVR_X = get_cov_cor(X)
  
  #get influence measures
  icphd_inf = sta_ipchd(X, mat_type = 'cor', mat = CoVR_X$Corr, var_inds = 1:args$p, overest_type = 'frac')
  s = 1 - icphd_inf
  n_nh = 0
  if(args$p == args$T0){
    n_nh = args$p - args$r
  }else{
    n_nh = args$p - args$T0
  }
  measure = linear_trans(s, n_hub = args$r, n_nh = n_nh)
  w_icphd = sig_prop(args$p, measure)
  
  #get optimal lambda
  
  max_rho = tp.weightedGL(X = X, .W = w_icphd, p = args$p, n = args$n)
  lambda = lambda * max_rho
  
  #Fit GLASSO
  Theta_icphd_BIC = BICweighted_glasso(mat = CoVR_X$Cov, rho = lambda, .W = w_icphd, p = args$p, n = args$n,
                                       threshold = 0, maxit = 200,
                                       penalize.diagonal = FALSE)
  
  
  if(save){
    #Save results to analize
    
    folder_path = paste0(folder, '/LINEAR_ORACLE_BIC_',name)
    
    
    
    if (!dir.exists(folder_path)) {
      dir.create(folder_path)
    }
    #############################------------ Saving Results ------------------------------------------#
    # Save Influence Measure
    
    png(paste0(folder_path,"/LINEAR_ORACLE_influence.png"), width = 800, height = 600)  # Set file name and dimensions
    plot(measure, main = 'Influence adjusted measure of correlation', ylab = 'adjusted ICPHD Measure', xlab = 'Variable')
    dev.off()  # Save and close the file
    
    png(paste0(folder_path,"/IPCHD_influence.png"), width = 800, height = 600)  # Set file name and dimensions
    plot(icphd_inf, main = 'Influence measure of correlation', ylab = 'ICPHD Measure', xlab = 'Variable')
    dev.off()  # Save and close the file
    
    png(paste0(folder_path,"/IPCHD_adj_curv influence.png"), width = 800, height = 600)  # Set file name and dimensions
    plot(s, measure, main = 'adjusted curve', ylab = 'ICPHD Measure', xlab = 'Variable')
    dev.off()  # Save and close the file
    
    
    #Save weights for a hub and non hub
    png(paste0(folder_path,"/Weight_hub_node_IPCHD.png"), width = 800, height = 600)
    plot(w_icphd[1, 1:args$p], main = paste('Weights for IPCHD', name), xlab = 'Variable', ylab = 'Weight of of entry 1,j', 
         ylim = c(min(w_icphd[args$p, 1:args$p]), quantile(w_icphd[args$p, 1:args$p], 0.95)))
    points(w_icphd[args$p, 1:args$p], col = 'blue')
    legend("topright", legend = c("Hub node - 1", "Non-hub last node"), col = c("black", "blue"), pch = 19)
    dev.off()
    
    p <- plot_ly(z = w_icphd, type = "heatmap", colorscale = "Viridis")
    
    # Save the Plotly plot to an HTML file
    html_file <- paste0(folder_path,"/Weight_matrix.html")
    htmlwidgets::saveWidget(p, html_file)
    
    # Use webshot to capture the HTML file and save as PNG
    webshot(html_file, file = paste0(folder_path,"/Weight_matrix.png"))
    
    #Save Sanity check for fitting method
    
    png(paste0(folder_path,"/BIC_curve_IPCHD.png"), width = 800, height = 600)
    plot(x = lambda, y = Theta_icphd_BIC$BIC, log = 'x', main ='BIC for IPCHD', ylab = 'BIC', xlab = 'Lambda')
    dev.off()
    
    #Save matrices visualizations
    
    png(paste0(folder_path,"/Absolute_diff_IPCHD.png"), width = 800, height = 600)
    plot.pdm(abs(Theta - Theta_icphd_BIC$optimal.model$wi), main = 'Absolute Difference Estimated IPCHD')
    dev.off()
    
    png(paste0(folder_path,'/Precision_estimated_IPCHD.png'), width= 800, height = 600)
    plot.pdm(Theta_icphd_BIC$optimal.model$wi, main = 'Estimated IPCHD')
    dev.off()
  }  
  
  # Save Comparisons
  raw_results = all_comparisons(Theta, Theta_icphd_BIC$optimal.model$wi, args$r, threshold = threshold)
  return(list(results = raw_results, fit = Theta_icphd_BIC))
}


#GLasso simulation
simulation_Lik_2step = function(Theta, args,lambda, threshold, folder, save=FALSE){
  set.seed(args$seed)
  #simulate data
  sim_data = sample_from_precision(args$n, args$p, Theta)
  X_train = sim_data$X_train
  X_val = sim_data$X_val
  X = X_train
  
  #get correlation and covariance matrix
  CoVR_X_train = get_cov_cor(X_train)
  CoVR_X_val = get_cov_cor(X_val)
  
  #get influence measures
  W = matrix(1, args$p, args$p)
  diag(W) = 0
  
  #get optimal lambda
  
  max_rho = tp.weightedGL(X = X_train, .W = W, p = args$p, n = args$n)
  lambda = lambda * max_rho * 2
  
  #Fit GLASSO
  #should be the same as the one below
  #Theta_hw_BIC = BICweighted_glasso(mat = CoVR_X$Cov, rho = lambda, .W = w_hw, p = args$p, n = args$n,
  #                                  threshold = 0, maxit = 200,
  #                                  penalize.diagonal = TRUE)
  
  Theta_GL_BIC = Lik_2step_glasso(mat = CoVR_X_train$Cov, mat_val = CoVR_X_val$Cov, 
                                  rho = lambda, p = args$p, n = args$n,
                                  threshold = 0, maxit = 200,
                                  penalize.diagonal = FALSE)
  
  if(save){
  
    #Save results to analize
    folder_path = paste0(folder, '/2step_Lik')
    
    
    if (!dir.exists(folder_path)) {
      dir.create(folder_path)
    }
    #############################------------ Saving Results ------------------------------------------#
    # Save Influence Measure
    
    #Save Sanity check for fitting method
    
    png(paste0(folder_path,"/Lik_curve_2step.png"), width = 800, height = 600)
    persp(log(lambda), log(lambda), Theta_GL_BIC$BIC, 
          xlab = "log Lambda 1", ylab = "log Lambda 2", zlab = "-Lik",
          theta = 30, phi = 20, col = "lightblue", main = "3D BIC Plot")
    dev.off()
    
    png(paste0(folder_path,"/Lik_contour.png"), width = 800, height = 600)
    filled.contour(log(lambda), log(lambda), Theta_GL_BIC$BIC,
                   xlab = "log lambda 1 ", ylab = "log lambda 2",
                   main = "Contour Plot of BIC", color.palette = terrain.colors)
    dev.off()
    #Save matrices visualizations
    
    png(paste0(folder_path,"/Absolute_diff.png"), width = 800, height = 600)
    plot.pdm(abs(Theta - Theta_GL_BIC$optimal.model$wi), main = 'Absolute Difference Estimated 2step')
    dev.off()
    
    png(paste0(folder_path,'/Precision_estimated.png'), width= 800, height = 600)
    plot.pdm(Theta_GL_BIC$optimal.model$wi, main = 'Estimated 2step')
    dev.off()
  }  
  
  # Save Comparisons
  raw_results = all_comparisons(Theta, Theta_GL_BIC$optimal.model$wi, args$r, threshold = threshold)
  return(list(results = raw_results, fit = Theta_GL_BIC))
}


simulation_Lik_sig_adj = function(Theta, args, lambda, threshold, folder, name="", save=""){
  set.seed(args$seed)
  #simulate data
  sim_data = sample_from_precision(args$n, args$p, Theta)
  X_train = sim_data$X_train
  X = X_train
  
  X_val = sim_data$X_val
  X = X_train
  
  #get correlation and covariance matrix
  CoVR_X_train = get_cov_cor(X_train)
  CoVR_X_val = get_cov_cor(X_val)
  CoVR_X = CoVR_X_train
  
  #get influence measures
  icphd_inf = sta_ipchd(X_train, mat_type = 'cor', mat = CoVR_X_train$Corr, var_inds = 1:args$p, overest_type = 'frac')
  s = 1 - icphd_inf
  measure = get_sig_tr(s)
  w_icphd = sig_prop(args$p, measure)
  
  #get optimal lambda
  
  max_rho = 3
  lambda = lambda * max_rho
  
  #Fit GLASSO
  Theta_icphd_BIC = Lik_weighted_glasso(mat = CoVR_X$Cov, mat_val = CoVR_X_val$Cov, rho = lambda, .W = w_icphd, p = args$p, n = args$n,
                                       threshold = 0, maxit = 200,
                                       penalize.diagonal = FALSE)
  
  if(save){
    #Save results to analize
    
    folder_path = paste0(folder, '/SIG_ADJ_Lik_',name)
    
    
    
    if (!dir.exists(folder_path)) {
      dir.create(folder_path)
    }
    #############################------------ Saving Results ------------------------------------------#
    # Save Influence Measure
    
    png(paste0(folder_path,"/Lik_ADJ_influence.png"), width = 800, height = 600)  # Set file name and dimensions
    plot(measure, main = 'Influence adjusted measure of correlation', ylab = 'adjusted ICPHD Measure', xlab = 'Variable')
    dev.off()  # Save and close the file
    
    png(paste0(folder_path,"/IPCHD_influence.png"), width = 800, height = 600)  # Set file name and dimensions
    plot(icphd_inf, main = 'Influence measure of correlation', ylab = 'ICPHD Measure', xlab = 'Variable')
    dev.off()  # Save and close the file
    
    png(paste0(folder_path,"/IPCHD_adj_curv influence.png"), width = 800, height = 600)  # Set file name and dimensions
    plot(s, measure, main = 'adjusted curve', ylab = 'ICPHD Measure', xlab = 'Variable')
    dev.off()  # Save and close the file
    
    #Save weights for a hub and non hub
    png(paste0(folder_path,"/Weight_hub_node_IPCHD.png"), width = 800, height = 600)
    plot(w_icphd[1, 1:args$p], main = paste('Weights for IPCHD', name), xlab = 'Variable', ylab = 'Weight of of entry 1,j', 
         ylim = c(min(w_icphd[args$p, 1:args$p]), quantile(w_icphd[args$p, 1:args$p], 0.95)))
    points(w_icphd[args$p, 1:args$p], col = 'blue')
    legend("topright", legend = c("Hub node - 1", "Non-hub last node"), col = c("black", "blue"), pch = 19)
    dev.off()
    
    p <- plot_ly(z = w_icphd, type = "heatmap", colorscale = "Viridis")
    
    # Save the Plotly plot to an HTML file
    html_file <- paste0(folder_path,"/Weight_matrix.html")
    htmlwidgets::saveWidget(p, html_file)
    
    # Use webshot to capture the HTML file and save as PNG
    webshot(html_file, file = paste0(folder_path,"/Weight_matrix.png"))
    
    
    
    #Save Sanity check for fitting method
    
    png(paste0(folder_path,"/Lik_curve_IPCHD.png"), width = 800, height = 600)
    plot(x = lambda, y = Theta_icphd_BIC$BIC, log = 'x', main ='BIC for IPCHD', ylab = '-Lik', xlab = 'Lambda')
    dev.off()
    
    #Save matrices visualizations
    
    png(paste0(folder_path,"/Absolute_diff_IPCHD.png"), width = 800, height = 600)
    plot.pdm(abs(Theta - Theta_icphd_BIC$optimal.model$wi), main = 'Absolute Difference Estimated IPCHD')
    dev.off()
    
    png(paste0(folder_path,'/Precision_estimated_IPCHD.png'), width= 800, height = 600)
    plot.pdm(Theta_icphd_BIC$optimal.model$wi, main = 'Estimated IPCHD')
    dev.off()
  }  
  
  # Save Comparisons
  raw_results = all_comparisons(Theta, Theta_icphd_BIC$optimal.model$wi, args$r, threshold = threshold)
  return(list(results = raw_results, fit = Theta_icphd_BIC))
}


simulation_Lik_lin_adj = function(Theta, args, lambda, threshold, folder, name="", save=FALSE){
  set.seed(args$seed)
  #simulate data
  sim_data = sample_from_precision(args$n, args$p, Theta)
  X_train = sim_data$X_train
  X = X_train
  
  X_val = sim_data$X_val
  X = X_train
  
  #get correlation and covariance matrix
  CoVR_X_train = get_cov_cor(X_train)
  CoVR_X_val = get_cov_cor(X_val)
  CoVR_X = CoVR_X_train
  
  #get influence measures
  icphd_inf = sta_ipchd(X, mat_type = 'cor', mat = CoVR_X$Corr, var_inds = 1:args$p, overest_type = 'frac')
  s = 1 - icphd_inf
  measure = get_linear_tr(s)
  w_icphd = sig_prop(args$p, measure)
  
  #get optimal lambda
  
  max_rho = tp.weightedGL(X = X, .W = w_icphd, p = args$p, n = args$n)
  lambda = lambda * max_rho
  
  #Fit GLASSO
  Theta_icphd_BIC = Lik_weighted_glasso(mat = CoVR_X$Cov, mat_val = CoVR_X_val$Cov, rho = lambda, .W = w_icphd, p = args$p, n = args$n,
                                       threshold = 0, maxit = 200,
                                       penalize.diagonal = FALSE)
  #Save results to analize
  
  if(save){
    folder_path = paste0(folder, '/LINEAR_ADJ_Lik_',name)
    
    
    
    if (!dir.exists(folder_path)) {
      dir.create(folder_path)
    }
    #############################------------ Saving Results ------------------------------------------#
    # Save Influence Measure
    
    png(paste0(folder_path,"/LINEAR_ADJ_influence.png"), width = 800, height = 600)  # Set file name and dimensions
    plot(measure, main = 'Influence adjusted measure of correlation', ylab = 'adjusted ICPHD Measure', xlab = 'Variable')
    dev.off()  # Save and close the file
    
    png(paste0(folder_path,"/IPCHD_adj_curv influence.png"), width = 800, height = 600)  # Set file name and dimensions
    plot(s, measure, main = 'adjusted curve', ylab = 'ICPHD Measure', xlab = 'Variable')
    dev.off()  # Save and close the file
    
    
    png(paste0(folder_path,"/IPCHD_influence.png"), width = 800, height = 600)  # Set file name and dimensions
    plot(icphd_inf, main = 'Influence measure of correlation', ylab = 'ICPHD Measure', xlab = 'Variable')
    dev.off()  # Save and close the file
    
    #Save weights for a hub and non hub
    png(paste0(folder_path,"/Weight_hub_node_IPCHD.png"), width = 800, height = 600)
    plot(w_icphd[1, 1:args$p], main = paste('Weights for IPCHD', name), xlab = 'Variable', ylab = 'Weight of of entry 1,j', 
         ylim = c(min(w_icphd[args$p, 1:args$p]), quantile(w_icphd[args$p, 1:args$p], 0.95)))
    points(w_icphd[args$p, 1:args$p], col = 'blue')
    legend("topright", legend = c("Hub node - 1", "Non-hub last node"), col = c("black", "blue"), pch = 19)
    dev.off()
    
    p <- plot_ly(z = w_icphd, type = "heatmap", colorscale = "Viridis")
    
    # Save the Plotly plot to an HTML file
    html_file <- paste0(folder_path,"/Weight_matrix.html")
    htmlwidgets::saveWidget(p, html_file)
    
    # Use webshot to capture the HTML file and save as PNG
    webshot(html_file, file = paste0(folder_path,"/Weight_matrix.png"))
    
    #Save Sanity check for fitting method
    
    png(paste0(folder_path,"/Lik_curve_IPCHD.png"), width = 800, height = 600)
    plot(x = lambda, y = Theta_icphd_BIC$BIC, log = 'x', main ='BIC for IPCHD', ylab = 'Lik', xlab = 'Lambda')
    dev.off()
    
    #Save matrices visualizations
    
    png(paste0(folder_path,"/Absolute_diff_IPCHD.png"), width = 800, height = 600)
    plot.pdm(abs(Theta - Theta_icphd_BIC$optimal.model$wi), main = 'Absolute Difference Estimated IPCHD')
    dev.off()
    
    png(paste0(folder_path,'/Precision_estimated_IPCHD.png'), width= 800, height = 600)
    plot.pdm(Theta_icphd_BIC$optimal.model$wi, main = 'Estimated IPCHD')
    dev.off()
  }  
  
  # Save Comparisons
  raw_results = all_comparisons(Theta, Theta_icphd_BIC$optimal.model$wi, args$r, threshold = threshold)
  return(list(results = raw_results, fit = Theta_icphd_BIC))
}



simulation_Lik_lin_oracle = function(Theta, args, lambda, threshold, folder, name = "", save=FALSE){
  set.seed(args$seed)
  #simulate data
  sim_data = sample_from_precision(args$n, args$p, Theta)
  X = sim_data$X_train
  X_val = sim_data$X_val
  
  #get correlation and covariance matrix
  CoVR_X = get_cov_cor(X)
  CoVR_X_val = get_cov_cor(X_val)
  
  
  #get influence measures
  icphd_inf = sta_ipchd(X, mat_type = 'cor', mat = CoVR_X$Corr, var_inds = 1:args$p, overest_type = 'frac')
  s = 1 - icphd_inf
  n_nh = 0
  if(args$p == args$T0){
    n_nh = args$p - args$r
  }else{
    n_nh = args$p - args$T0
  }
  measure = linear_trans(s, n_hub = args$r, n_nh = n_nh)
  w_icphd = sig_prop(args$p, measure)
  
  #get optimal lambda
  
  max_rho = tp.weightedGL(X = X, .W = w_icphd, p = args$p, n = args$n)
  lambda = lambda * max_rho
  
  #Fit GLASSO
  Theta_icphd_BIC = Lik_weighted_glasso(mat = CoVR_X$Cov, mat_val = CoVR_X_val$Cov, 
                                       rho = lambda, .W = w_icphd, p = args$p, n = args$n,
                                       threshold = 0, maxit = 200,
                                       penalize.diagonal = FALSE)
  if(save){
    #Save results to analize
    
    folder_path = paste0(folder, '/LINEAR_ORACLE_Lik_',name)
    
    
    
    if (!dir.exists(folder_path)) {
      dir.create(folder_path)
    }
    #############################------------ Saving Results ------------------------------------------#
    # Save Influence Measure
    
    png(paste0(folder_path,"/LINEAR_ORACLE_influence.png"), width = 800, height = 600)  # Set file name and dimensions
    plot(measure, main = 'Influence adjusted measure of correlation', ylab = 'adjusted ICPHD Measure', xlab = 'Variable')
    dev.off()  # Save and close the file
    
    png(paste0(folder_path,"/IPCHD_influence.png"), width = 800, height = 600)  # Set file name and dimensions
    plot(icphd_inf, main = 'Influence measure of correlation', ylab = 'ICPHD Measure', xlab = 'Variable')
    dev.off()  # Save and close the file
    
    png(paste0(folder_path,"/IPCHD_adj_curv influence.png"), width = 800, height = 600)  # Set file name and dimensions
    plot(s, measure, main = 'adjusted curve', ylab = 'ICPHD Measure', xlab = 'Variable')
    dev.off()  # Save and close the file
    
    
    #Save weights for a hub and non hub
    png(paste0(folder_path,"/Weight_hub_node_IPCHD.png"), width = 800, height = 600)
    plot(w_icphd[1, 1:args$p], main = paste('Weights for IPCHD', name), xlab = 'Variable', ylab = 'Weight of of entry 1,j', 
         ylim = c(min(w_icphd[args$p, 1:args$p]), quantile(w_icphd[args$p, 1:args$p], 0.95)))
    points(w_icphd[args$p, 1:args$p], col = 'blue')
    legend("topright", legend = c("Hub node - 1", "Non-hub last node"), col = c("black", "blue"), pch = 19)
    dev.off()
    
    p <- plot_ly(z = w_icphd, type = "heatmap", colorscale = "Viridis")
    
    # Save the Plotly plot to an HTML file
    html_file <- paste0(folder_path,"/Weight_matrix.html")
    htmlwidgets::saveWidget(p, html_file)
    
    # Use webshot to capture the HTML file and save as PNG
    webshot(html_file, file = paste0(folder_path,"/Weight_matrix.png"))
    
    #Save Sanity check for fitting method
    
    png(paste0(folder_path,"/Lik_curve_IPCHD.png"), width = 800, height = 600)
    plot(x = lambda, y = Theta_icphd_BIC$BIC, log = 'x', main ='BIC for IPCHD', ylab = '-Lik', xlab = 'Lambda')
    dev.off()
    
    #Save matrices visualizations
    
    png(paste0(folder_path,"/Absolute_diff_IPCHD.png"), width = 800, height = 600)
    plot.pdm(abs(Theta - Theta_icphd_BIC$optimal.model$wi), main = 'Absolute Difference Estimated IPCHD')
    dev.off()
    
    png(paste0(folder_path,'/Precision_estimated_IPCHD.png'), width= 800, height = 600)
    plot.pdm(Theta_icphd_BIC$optimal.model$wi, main = 'Estimated IPCHD')
    dev.off()
  }  
  
  # Save Comparisons
  raw_results = all_comparisons(Theta, Theta_icphd_BIC$optimal.model$wi, args$r, threshold = threshold)
  return(list(results = raw_results, fit = Theta_icphd_BIC))
}


simulation_Lik_Full_oracle = function(Theta, args,lambda, threshold, folder, save=FALSE){
  set.seed(args$seed)
  #simulate data
  sim_data = sample_from_precision(args$n, args$p, Theta)
  X_train = sim_data$X_train
  X_val = sim_data$X_val
  X = X_train
  
  #get correlation and covariance matrix
  CoVR_X = get_cov_cor(X_train)
  CoVR_X_val = get_cov_cor(X_val)
  
  #get influence measures
  #get influence measures
  W = 1 - (abs(Theta) >= 1e-5) * 1
  
  #get optimal lambda
  
  max_rho = tp.weightedGL(X = X, .W = W, p = args$p, n = args$n)
  lambda = lambda * max_rho
  
  #Fit GLASSO
  Theta_GL_Lik = Lik_weighted_glasso(mat = CoVR_X$Cov, mat_val = CoVR_X_val$Cov, rho = lambda, 
                                     .W = W , p = args$p, n = args$n,
                                     threshold = 0, maxit = 200,
                                     penalize.diagonal = FALSE)
  if(save){
    #Save results to analize
    folder_path = paste0(folder, '/Full_oracle_Lik')
    
    
    if (!dir.exists(folder_path)) {
      dir.create(folder_path)
    }
    #############################------------ Saving Results ------------------------------------------#
    
    #Save weights for a hub and non hub
    png(paste0(folder_path,"/Weight_hub_node_oracle.png"), width = 800, height = 600)
    plot(W[1, 1:args$p], main = 'Weights for GL', xlab = 'Variable', ylab = 'Weight of of entry 1,j', 
         ylim = c(min(W[args$p, 1:args$p]), quantile(W[args$p, 1:args$p], 0.95)))
    points(W[args$p, 1:args$p], col = 'blue')
    legend("topright", legend = c("Hub node - 1", "Non-hub last node"), col = c("black", "blue"), pch = 19)
    dev.off()
    
    p <- plot_ly(z = W, type = "heatmap", colorscale = "Viridis")
    
    # Save the Plotly plot to an HTML file
    html_file <- paste0(folder_path,"/Weight_matrix.html")
    htmlwidgets::saveWidget(p, html_file)
    
    # Use webshot to capture the HTML file and save as PNG
    webshot(html_file, file = paste0(folder_path,"/Weight_matrix.png"))
    
    
    #Save Sanity check for fitting method
    
    png(paste0(folder_path,"/Likelihood_curve_oracle.png"), width = 800, height = 600)
    plot(x = lambda, y = Theta_GL_Lik$BIC, log = 'x', main ='Lik', ylab = 'BIC', xlab = 'Lambda')
    dev.off()
    
    #Save matrices visualizations
    
    png(paste0(folder_path,"/Absolute_diff_oracle.png"), width = 800, height = 600)
    plot.pdm(abs(Theta - Theta_GL_Lik$optimal.model$wi), main = 'Absolute Difference Estimated oracle')
    dev.off()
    
    png(paste0(folder_path,'/Precision_estimated_oracle.png'), width= 800, height = 600)
    plot.pdm(Theta_GL_Lik$optimal.model$wi, main = 'Estimated oracle')
    dev.off()
  }  
  
  # Save Comparisons
  raw_results = all_comparisons(Theta, Theta_GL_Lik$optimal.model$wi, args$r, threshold = threshold)
  return(list(results = raw_results, fit = Theta_GL_Lik))
}

simulation_BIC_Full_oracle = function(Theta, args,lambda, threshold, folder, save=FALSE){
  set.seed(args$seed)
  sim_data = sample_from_precision(args$n, args$p, Theta)
  X = sim_data$X
  
  #get correlation and covariance matrix
  CoVR_X = get_cov_cor(X)
  
  #get influence measures
  #get influence measures
  W = 1 - (abs(Theta) >= 1e-5) * 1
  
  #get optimal lambda
  
  max_rho = tp.weightedGL(X = X, .W = W, p = args$p, n = args$n)
  lambda = lambda * max_rho
  
  #Fit GLASSO
  Theta_GL_Lik = BICweighted_glasso(mat = CoVR_X$Cov, rho = lambda, 
                                     .W = W , p = args$p, n = args$n,
                                     threshold = 0, maxit = 200,
                                     penalize.diagonal = FALSE)
  
  if(save){
    #Save results to analize
    folder_path = paste0(folder, '/Full_oracle_BIC')
    
    
    if (!dir.exists(folder_path)) {
      dir.create(folder_path)
    }
    #############################------------ Saving Results ------------------------------------------#
    
    #Save weights for a hub and non hub
    png(paste0(folder_path,"/Weight_hub_node_oracle.png"), width = 800, height = 600)
    plot(W[1, 1:args$p], main = 'Weights for GL', xlab = 'Variable', ylab = 'Weight of of entry 1,j', 
         ylim = c(min(W[args$p, 1:args$p]), quantile(W[args$p, 1:args$p], 0.95)))
    points(W[args$p, 1:args$p], col = 'blue')
    legend("topright", legend = c("Hub node - 1", "Non-hub last node"), col = c("black", "blue"), pch = 19)
    dev.off()
    
    p <- plot_ly(z = W, type = "heatmap", colorscale = "Viridis")
    
    # Save the Plotly plot to an HTML file
    html_file <- paste0(folder_path,"/Weight_matrix.html")
    htmlwidgets::saveWidget(p, html_file)
    
    # Use webshot to capture the HTML file and save as PNG
    webshot(html_file, file = paste0(folder_path,"/Weight_matrix.png"))
    
    
    #Save Sanity check for fitting method
    
    png(paste0(folder_path,"/Likelihood_curve_oracle.png"), width = 800, height = 600)
    plot(x = lambda, y = Theta_GL_Lik$BIC, log = 'x', main ='BIC', ylab = 'BIC', xlab = 'Lambda')
    dev.off()
    
    #Save matrices visualizations
    
    png(paste0(folder_path,"/Absolute_diff_oracle.png"), width = 800, height = 600)
    plot.pdm(abs(Theta - Theta_GL_Lik$optimal.model$wi), main = 'Absolute Difference Estimated oracle')
    dev.off()
    
    png(paste0(folder_path,'/Precision_estimated_oracle.png'), width= 800, height = 600)
    plot.pdm(Theta_GL_Lik$optimal.model$wi, main = 'Estimated oracle')
    dev.off()
  }  
  
  # Save Comparisons
  raw_results = all_comparisons(Theta, Theta_GL_Lik$optimal.model$wi, args$r, threshold = threshold)
  return(list(results = raw_results, fit = Theta_GL_Lik))
}
