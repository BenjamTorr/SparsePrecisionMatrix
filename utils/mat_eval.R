TPR_TNR = function(Theta_true, Theta_est, threshold = 1e-6){
  p = nrow(Theta_true)
  Theta0_true = 1 * (abs(Theta_true) > threshold)
  Theta0_est = 1 * (abs(Theta_est) > threshold)
  TN = 0
  TNFP = 0
  TP = 0
  TPFN = 0
  for(i in 1:p){
    for(j in i:p){
      if(i < j){
        TN = TN +  (Theta0_true[i,j] == 0) * (Theta0_est[i,j] == 0)
        TNFP = TNFP + (Theta0_true[i,j] == 0)
      }
      TP = TP + (Theta0_true[i,j] != 0) * (Theta0_est[i,j] != 0)
      TPFN = TPFN + (Theta0_true[i,j] != 0)
    }
  }
  return(list(TNR = (TN / TNFP), TPR = (TP/TPFN)))
}


Frob_norm = function(mat, mat_est){
  return(sqrt(sum((mat - mat_est) ** 2)))
}

spectral_norm <- function(mat, mat_est) {
  # Calculate the difference matrix
  diff_matrix <- mat - mat_est
  
  # Compute the largest singular value (spectral norm)
  norm <- max(svd(diff_matrix)$d)  # $d gives the singular values
  
  return(norm)
}


log_det <- function(matrix) {
  eigenvalues <- eigen(matrix, symmetric = TRUE, only.values = TRUE)$values
  return(sum(log(eigenvalues)))
}

# Function to compute KL divergence between two multivariate Gaussians
kl_divergence <- function(Theta1, Theta2) {
  
  # Calculate matrix dimensions
  d <- nrow(Theta1)
  
  # Compute trace of (Theta2^-1 * Theta1)
  Theta2_inv <- solve(Theta2)
  trace_term <- sum(diag(Theta2_inv %*% Theta1))
  
  # Compute log(det(Theta2) / det(Theta1)) using eigenvalues
  log_det_ratio <- log_det(Theta2) - log_det(Theta1)
  
  # Compute KL divergence
  kl <- 0.5 * (trace_term - d + log_det_ratio)
  
  return(kl)
}

# Function to compute sparsity pattern comparison metrics
sparsity_comparison <- function(Theta1, Theta2, threshold = 1e-5) {
  # Binarize the precision matrices based on threshold
  sparsity1 <- abs(Theta1) > threshold
  sparsity2 <- abs(Theta2) > threshold
  
  # Ensure the matrices are of the same size
  if (!all(dim(Theta1) == dim(Theta2))) {
    stop("Precision matrices must have the same dimensions.")
  }
  
  # Compute Jaccard similarity: |A ∩ B| / |A ∪ B|
  intersection <- sum(sparsity1 & sparsity2)
  union <- sum(sparsity1 | sparsity2)
  jaccard_similarity <- intersection / union
  
  # Compute Hamming distance: |A XOR B| / Total Elements
  xor <- sparsity1 != sparsity2
  hamming_distance <- sum(xor) / length(xor)
  
  return(list(
    Jaccard_Similarity = jaccard_similarity,
    Hamming_Distance = hamming_distance
  ))
}


# Function to compute true positives for non-zero and zero entries
compare_TPR_TNR <- function(inv_mat, inv_mat_est, threshold = 1e-5) {
  # Ensure matrices have the same dimensions
  if (!all(dim(inv_mat) == dim(inv_mat_est))) {
    stop("True and estimated matrices must have the same dimensions.")
  }
  
  sparsity1 <- abs(inv_mat) > threshold
  sparsity2 <- abs(inv_mat_est) > threshold
  
  true_zero_est = (1 - sparsity1) * (1 - sparsity2)
  all_real_zero = (1 - sparsity1)
  
  TNR = sum(true_zero_est[upper.tri(true_zero_est)]) / sum(all_real_zero[upper.tri(all_real_zero)])
  
  true_nzero_est = sparsity1 * sparsity2
  all_real_nzero = sparsity1
  
  TPR = (sum(true_nzero_est[upper.tri(true_nzero_est, diag = TRUE)]) / 
           sum(all_real_nzero[upper.tri(all_real_nzero, diag = TRUE)])
  )
  
  
  return(list(
    True_Positive_Rate = TPR * 100,
    True_Negative_Rate = TNR * 100
  ))
}

TPR_TNR_hub = function(Theta_true, Theta_est, r, threshold = 1e-6){
  p = nrow(Theta_true)
  Theta0_true = 1 * (abs(Theta_true) > threshold)
  Theta0_est = 1 * (abs(Theta_est) > threshold)
  TN = 0
  TNFP = 0
  TP = 0
  TPFN = 0
  for(i in 1:r){
    for(j in i:p){
      if(i < j){
        TN = TN +  (Theta0_true[i,j] == 0) * (Theta0_est[i,j] == 0)
        TNFP = TNFP + (Theta0_true[i,j] == 0)
      }
      TP = TP + (Theta0_true[i,j] != 0) * (Theta0_est[i,j] != 0)
      TPFN = TPFN + (Theta0_true[i,j] != 0)
    }
  }
  if(TNFP != 0){
    TNR = TN / TNFP  
  }else{
    TNR = 1
  }
  if(TPFN != 0){
    TPR = TP/TPFN  
  }else{
    TPR = 1
  }
  
  return(list(TNR = TNR, TPR = TPR))
}


TPR_TNR_nonhub = function(Theta_true, Theta_est, r, threshold = 1e-6){
  p = nrow(Theta_true)
  Theta0_true = 1 * (abs(Theta_true) > threshold)
  Theta0_est = 1 * (abs(Theta_est) > threshold)
  TN = 0
  TNFP = 0
  TP = 0
  TPFN = 0
  for(i in (r + 1):p){
    for(j in i:p){
      if(i < j){
        TN = TN +  (Theta0_true[i,j] == 0) * (Theta0_est[i,j] == 0)
        TNFP = TNFP + (Theta0_true[i,j] == 0)
      }
      TP = TP + (Theta0_true[i,j] != 0) * (Theta0_est[i,j] != 0)
      TPFN = TPFN + (Theta0_true[i,j] != 0)
    }
  }
  if(TNFP != 0){
    TNR = TN / TNFP  
  }else{
    TNR = 1
  }
  if(TPFN != 0){
    TPR = TP/TPFN  
  }else{
    TPR = 1
  }
  
  return(list(TNR = TNR, TPR = TPR))
}

hub_nonhub_accuracy <- function(precision_matrix, r) {
  # Number of nodes
  n <- nrow(precision_matrix)
  threshold <- 0.1 * (n - 1)  # More than 10% of nodes
  
  # Identify hub nodes based on nonzero connections (ignore diagonal)
  connectivity <- rowSums(abs(precision_matrix) >= 1e-5) - 1  # Exclude self-connections
  estimated_hub_nodes <- which(connectivity > threshold)
  
  # Define true hub and non-hub nodes
  true_hub_nodes <- 1:r
  true_nonhub_nodes <- setdiff(1:n, true_hub_nodes)
  
  # Estimated non-hub nodes (complement of estimated hub nodes)
  estimated_nonhub_nodes <- setdiff(1:n, estimated_hub_nodes)
  
  # Correctly identified hub nodes
  correctly_identified_hub <- sum(estimated_hub_nodes %in% true_hub_nodes)
  hub_accuracy <- if (r == 0) ifelse(length(estimated_hub_nodes) == 0, 1.0, 0.0) else correctly_identified_hub / r
  
  # Correctly identified non-hub nodes
  correctly_identified_nonhub <- sum(estimated_nonhub_nodes %in% true_nonhub_nodes)
  nonhub_accuracy <- if ((n - r) == 0) ifelse(length(estimated_nonhub_nodes) == 0, 1.0, 0.0) else correctly_identified_nonhub / (n - r)
  
  return(list(hub_accuracy = hub_accuracy, nonhub_accuracy = nonhub_accuracy))
}



all_comparisons = function(inv_mat, inv_mat_est,   r, threshold = 1e-5){
  
  Frobenius = Frob_norm(inv_mat, inv_mat_est)
  spectral = spectral_norm(inv_mat, inv_mat_est)
  kl = kl_divergence(inv_mat, inv_mat_est)
  jaccard_hamming = sparsity_comparison(inv_mat, inv_mat_est, threshold = threshold)
  TPR_TNR = compare_TPR_TNR(inv_mat, inv_mat_est, threshold = threshold)
  TPR_TNR_hub_res = TPR_TNR_hub(inv_mat, inv_mat_est, r = r, threshold = threshold)
  TPR_TNR_nonhub_res = TPR_TNR_nonhub(inv_mat, inv_mat_est, r = r, threshold = threshold)
  #hub_nonhub = hub_nonhub_accuracy(inv_mat, r)
  
  
  
  return(list(Frobenius = Frobenius,
              spectral_norm = spectral,
              kl_divergence = kl,
              jaccard = jaccard_hamming$Jaccard_Similarity,
              hamming = jaccard_hamming$Hamming_Distance,
              TPR = TPR_TNR$True_Positive_Rate,
              TNR = TPR_TNR$True_Negative_Rate,
              TPR_hub = TPR_TNR_hub_res$TPR,
              TNR_hub = TPR_TNR_hub_res$TNR,
              TPR_nhub = TPR_TNR_nonhub_res$TPR,
              TNR_nhub = TPR_TNR_nonhub_res$TNR))
}
