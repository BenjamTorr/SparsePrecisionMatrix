plot_metric_comparison = function(full_results, n_sample, metric){
  n_sim = length(n_sample)
  chosen_metric = c()
  for(method in names(full_results)){
    for(j in 1:n_sim){
      if(is.null(full_results[[method]][[j]][[metric]])){
        print(method)
      }
      chosen_metric = c(chosen_metric, full_results[[method]][[j]][[metric]])
    }
  }
  # Sample dataset
  df <- data.frame(
    Method = rep(names(full_results), each = n_sim),
    Sample_Size = rep(1:length(n_sample), times = length(full_results)),
    Metric = chosen_metric
  )

  # Compute IQR-based limits to remove outliers
  q1 <- quantile(df$Metric, 0.99, na.rm = TRUE)
  upper_bound <- q1
  
  # Filter out extreme outliers
  df <- df[df$Metric <= upper_bound, ]
  
  # Create the plot with distinct linetypes and shapes
  ggplot(df, aes(x = Sample_Size, y = Metric, color = Method, group = Method)) +
    geom_line(aes(linetype = Method), size = 1) +  # Different linetypes
    geom_point(aes(shape = Method), size = 3) +   # Different point shapes
    scale_x_continuous(breaks = unique(df$Sample_Size), labels= n_sample) +
    labs(
      title = paste0("Comparison of Methods using ", metric),
      x = "sample size",
      y = metric,
      color = "Method",
      linetype = "Method",
      shape = "Method"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
      axis.text = element_text(size = 12),
      axis.title = element_text(size = 12),
      legend.text = element_text(size = 10)
    )
  
}
