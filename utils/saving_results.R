save_comparisons = function(experiment_name, n_sample, full_results, type = ''){
  png(paste0(experiment_name, '/comparison_', "Frobenius_", type ,'.png'), width = 800, height = 600)
  plot_metric_comparison(full_results, n_sample, "Frobenius")  
  dev.off()
  
  png(paste0(experiment_name, '/comparison_',  "spectral_norm_", type,'.png'), width = 800, height = 600)
  plot_metric_comparison(full_results, n_sample,  "spectral_norm")  
  dev.off()
  
  png(paste0(experiment_name, '/comparison_', 'TPR_', type,'.png'), width = 800, height = 600)
  plot_metric_comparison(full_results, n_sample, 'TPR')  
  dev.off()
  
  png(paste0(experiment_name, '/comparison_', "kl_divergence_", type,'.png'), width = 800, height = 600)
  plot_metric_comparison(full_results, n_sample, "kl_divergence")  
  dev.off()
  
  png(paste0(experiment_name, '/comparison_', "jaccard_", type ,'.png'), width = 800, height = 600)
  plot_metric_comparison(full_results, n_sample, "jaccard" )  
  dev.off()
  
  png(paste0(experiment_name, '/comparison_', "hamming_", type, '.png'), width = 800, height = 600)
  plot_metric_comparison(full_results, n_sample, "hamming")  
  dev.off()
  
  png(paste0(experiment_name, '/comparison_', 'TPR_', type,'.png'), width = 800, height = 600)
  plot_metric_comparison(full_results, n_sample, 'TPR')  
  dev.off()
  
  png(paste0(experiment_name, '/comparison_', 'TNR_', type,'.png'), width = 800, height = 600)
  plot_metric_comparison(full_results, n_sample, 'TNR')  
  dev.off()
  
  save(full_results, file = paste0(experiment_name,"/full_results_", type, ".RData"))
  
}