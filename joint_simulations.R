################################################################################
rm(list=ls())
setwd("D:/UNC/PhD UNC/Research/Sparse Precision Matrix/code")

invisible(lapply(list.files("JoseMethods", pattern = "\\.R$", full.names = TRUE), source))
invisible(lapply(list.files("utils", pattern = "\\.R$", full.names = TRUE), source))
invisible(lapply(list.files("HWGLASSO", pattern = "\\.R$", full.names = TRUE), source))
invisible(lapply(list.files('pretuning',pattern = "\\.R$", full.names = TRUE), source))
source('SimulationParameters.R')
source('simulation_functions.R')
library(ggplot2)
library(reshape2)
library(plotly)
library(matrixcalc)
library(webshot)
library(dplyr)
library(tidyr)
#---------------------------Hyperparameters-------------------------------------
# p_val = c(30, 30, 30, 50, 100)
# p_ph = c(0.9, 1, 0.85, 0.8, 0.9)
# p_pnh = c(0.05, 0.0, 0.05, 0.1, 0.07)
# p_r = c(5, 5, 5, 7, 10)
# case_exp_name = c("OverNight1_", "OverNight2_", "OverNight3_", "OverNight4_", "OverNight5_")
# 
# for(case.i in 1:5){
#   if(case.i > 1){
#     args$p = p_val[case.i]
#     args$ph = p_ph[case.i]
#     args$pnh = p_pnh[case.i]
#     args$r = p_r[case.i]
#     
#   }
  



test_name = paste0("Adjusted2_",args$p)

experiment_name = paste0('Simulations_new/', test_name)

if (!dir.exists(experiment_name)) {
  dir.create(experiment_name)
}

save(args, file = paste0(experiment_name,"/args.RData"))

weights_fun = list(inverse_prop = inverse_prop)

lambda = 10^seq(-5, 3, length = 50)
lambda_small = 10^seq(-5, 1, length = 50)
lambda_2step = 10^seq(-5, 1, length = 50)



threshold = 1e-5

n_props = c(1.1, 1.5, 2, 3, 5, 10 ,100, 500)

correlation = FALSE

rep_num = 20

#### generate and save simulation settings

set.seed(args$seed)

Theta = r.sparse.pdhubmat(
  args$p, T0 = args$T0, r = args$r, ph = args$ph, pnh = args$pnh, pneff = args$pneff, diagonal_shift = args$diagonal_shift, 
  shuffle = args$shuffle, type = "unif", hmin = args$hmin, hmax = args$hmax, nhmin = args$nhmin, nhmax = args$nhmax, 
  neffmin = args$neffmin, neffmax = args$neffmax)

cat("The condition number of the Precision Matrix is:", kappa(Theta))

png(paste0(experiment_name, '/Real_precision_matrix.png'), width = 800, height = 600)
plot.pdm(Theta, main = 'Real Precision Matrix')
dev.off()

png(paste0(experiment_name, '/Eigenvalues_precision.png'), width = 800, height = 600)
plot(eigen(Theta)$values, main = 'Eigenvalues of Precision matrix', xlab = 'Variable', ylab = 'Eigenvalue')
dev.off()

png(paste0(experiment_name, '/Eigenvalues_inverse_corr.png'), width = 800, height = 600)
plot(eigen(.PMtoIC(Theta))$values, main = 'Eigenvalues of inverse correlation matrix', xlab = 'Variable', ylab = 'Eigenvalue')
dev.off()

if(correlation){
  Theta = .PMtoIC(Theta)
}


#simulation

n_sample = floor(n_props * args$p)
n_sim = length(n_sample)

# full_results_BIC = list(inverse_prop = list(),
#                         hw = list(),
#                         adaGL = list(),
#                         GL = list(),
#                         hub_oracle = list(),
#                         lin_adj = list(),
#                         Full_oracle = list())
# 
# full_results_Lik = list(inverse_prop = list(),
#                         hw = list(),
#                         adaGL = list(),
#                         GL = list(),
#                         hub_oracle = list(),
#                         lin_adj = list(),
#                         Full_oracle = list())


full_simulations = run_simulations(experiment_name, n_sample, n_sim, rep_num = rep_num, Theta, args, lambda, lambda_2step, threshold, folder, save)
full_results_BIC = full_simulations$BIC
full_results_Lik = full_simulations$Lik



# full_results_BIC = list()
# full_results_Lik = list()
# 
# 
# 
# # save simulations results
# for(j in 1:n_sim){
#   args$seed = args$seed + 1
#   save = TRUE
#   args$n = n_sample[j]
#   cat('----------- ', args$n,' ----------------\n' )
#   folder = paste0(experiment_name,'/sample_', n_sample[j])
#   if (!dir.exists(folder)) {
#     dir.create(folder)
#   }
#   ##Real oracle
#   full_results_BIC$Full_oracle[[j]] = simulation_BIC_Full_oracle(Theta, args, lambda, threshold, folder, save=save)$results
#   full_results_Lik$Full_oracle[[j]] = simulation_Lik_Full_oracle(Theta, args, lambda, threshold, folder, save=save)$results
#   #adjusted
#   cat('\n Sig Adj\n')
# #  full_results_BIC$sig_adj[[j]] = simulation_BIC_sig_adj(Theta, args, lambda, threshold, folder, save)$results
# #  full_results_Lik$sig_adj[[j]] = simulation_Lik_sig_adj(Theta, args, lambda, threshold, folder)$results
#   
#   cat('\n Lin Adj \n')
#   full_results_BIC$lin_adj[[j]] = simulation_BIC_lin_adj(Theta, args, lambda, threshold, folder)$results
#   full_results_Lik$lin_adj[[j]] = simulation_Lik_lin_adj(Theta, args, lambda, threshold, folder)$results
#   
#   cat("\n Lin Oracle \n")
#   #full_results_BIC$lin_hub_oracle[[j]] = simulation_BIC_lin_oracle(Theta, args, lambda, threshold, folder, save=save)$results
#   #full_results_Lik$lin_hub_oracle[[j]] = simulation_Lik_lin_oracle(Theta, args, lambda, threshold, folder, save=save)$results
#   
#   #HWGLASSO
#   cat("\n Working on HWGLASSO\n")
#   full_results_BIC$hw[[j]] = simulation_BIC_hw(Theta, args, lambda, threshold, folder, save=save)$results
#   full_results_Lik$hw[[j]] = simulation_Lik_hw(Theta, args, lambda, threshold, folder, save=save)$results
#   
#   #Adaptive GLASSO
#   cat("\n working on Adapative GLASSO\n")
#   full_results_BIC$adaGL[[j]] = simulation_BIC_adaGL(Theta, args, lambda, threshold, folder, save=save)$results
#   full_results_Lik$adaGL[[j]] = simulation_Lik_adaGL(Theta, args, lambda, threshold, folder, save=save)$results
#   
#   #GLASSO no weight
#   cat("\n Working on pure GLASSO\n")
#   full_results_BIC$GL[[j]] = simulation_BIC_GL(Theta, args, lambda, threshold, folder, save=save)$results
#   full_results_Lik$GL[[j]] = simulation_Lik_GL(Theta, args, lambda, threshold, folder, save=save)$results
#   
#   #Glasso oracle
#   cat("\n Working on oracle\n")
#   full_results_BIC$hub_oracle[[j]] = simulation_BIC_oracle(Theta, args, lambda, threshold, folder, save=save)$results
#   full_results_Lik$hub_oracle[[j]] = simulation_Lik_oracle(Theta, args, lambda, threshold, folder, save=save)$results
#   
#   #two step
#   cat("Working on 2 step\n")
# #  full_results_BIC$two_step[[j]] = simulation_BIC_2step(Theta, args, lambda_2step, threshold, folder)$results
# #  full_results_Lik$two_step[[j]] = simulation_Lik_2step(Theta, args, lambda_2step, threshold, folder)$results  
#   
#   cat("Working on New weights\n")
#   for(func in names(weights_fun)){
#     full_results_BIC[[func]][[j]] = simulation_BIC_ipchd(Theta, args,  weights_fun[[func]] ,lambda, threshold, folder, name = func, save=save)$results
#     full_results_Lik[[func]][[j]] = simulation_Lik_ipchd(Theta, args, weights_fun[[func]] ,lambda, threshold, folder, name = func, save=save)$results
#     #full_results_BIC[[paste0(func,"_frob")]][[j]] = simulation_BIC_ipchd_frob(Theta, args, weights_fun[[func]] ,lambda, threshold, folder, name = func)$results
#   }
# }

#plot final results

type = 'BIC'

BIC_results_folder = paste0(experiment_name, '/BIC_results')

if (!dir.exists(BIC_results_folder)) {
  dir.create(BIC_results_folder)
}

png(paste0(BIC_results_folder, '/comparison_', "Frobenius_", type ,'.png'), width = 800, height = 600)
plot_metric_comparison(full_results_BIC, n_sample, "Frobenius")  
dev.off()

png(paste0(BIC_results_folder, '/comparison_',  "spectral_norm_", type,'.png'), width = 800, height = 600)
plot_metric_comparison(full_results_BIC, n_sample,  "spectral_norm")  
dev.off()

png(paste0(BIC_results_folder, '/comparison_', 'TPR_', type,'.png'), width = 800, height = 600)
plot_metric_comparison(full_results_BIC, n_sample, 'TPR')  
dev.off()

png(paste0(BIC_results_folder, '/comparison_', "kl_divergence_", type,'.png'), width = 800, height = 600)
plot_metric_comparison(full_results_BIC, n_sample, "kl_divergence")  
dev.off()

png(paste0(BIC_results_folder, '/comparison_', "jaccard_", type ,'.png'), width = 800, height = 600)
plot_metric_comparison(full_results_BIC, n_sample, "jaccard" )  
dev.off()

png(paste0(BIC_results_folder, '/comparison_', "hamming_", type, '.png'), width = 800, height = 600)
plot_metric_comparison(full_results_BIC, n_sample, "hamming")  
dev.off()

png(paste0(BIC_results_folder, '/comparison_', 'TPR_', type,'.png'), width = 800, height = 600)
plot_metric_comparison(full_results_BIC, n_sample, 'TPR')  
dev.off()

png(paste0(BIC_results_folder, '/comparison_', 'TNR_', type,'.png'), width = 800, height = 600)
plot_metric_comparison(full_results_BIC, n_sample, 'TNR')  
dev.off()

png(paste0(BIC_results_folder, '/comparison_', 'TPR_hub_', type,'.png'), width = 800, height = 600)
plot_metric_comparison(full_results_BIC, n_sample, 'TPR_hub')  
dev.off()

png(paste0(BIC_results_folder, '/comparison_', 'TNR_hub_', type,'.png'), width = 800, height = 600)
plot_metric_comparison(full_results_BIC, n_sample, 'TNR_hub')  
dev.off()

png(paste0(BIC_results_folder, '/comparison_', 'TPR_nonhub_', type,'.png'), width = 800, height = 600)
plot_metric_comparison(full_results_BIC, n_sample, 'TPR_nhub')  
dev.off()

png(paste0(BIC_results_folder, '/comparison_', 'TNR_nonhub_', type,'.png'), width = 800, height = 600)
plot_metric_comparison(full_results_BIC, n_sample, 'TNR_nhub')  
dev.off()


# png(paste0(BIC_results_folder, '/comparison_', 'hub_accuracy_', type,'.png'), width = 800, height = 600)
# plot_metric_comparison(full_results_BIC, n_sample, 'hub_accuracy')  
# dev.off()
# 
# png(paste0(BIC_results_folder, '/comparison_', 'nonhub_accuracy_', type,'.png'), width = 800, height = 600)
# plot_metric_comparison(full_results_BIC, n_sample, 'nonhub_accuracy')  
# dev.off()


save(full_results_BIC, file = paste0(BIC_results_folder,"/full_results_", type, ".RData"))

 
type = 'Lik'

Lik_results_folder = paste0(experiment_name, '/Lik_results')

if (!dir.exists(Lik_results_folder)) {
  dir.create(Lik_results_folder)
}

png(paste0(Lik_results_folder, '/comparison_', "Frobenius_", type ,'.png'), width = 800, height = 600)
plot_metric_comparison(full_results_Lik, n_sample, "Frobenius")
dev.off()

png(paste0(Lik_results_folder, '/comparison_',  "spectral_norm_", type,'.png'), width = 800, height = 600)
plot_metric_comparison(full_results_Lik, n_sample,  "spectral_norm")
dev.off()

png(paste0(Lik_results_folder, '/comparison_', 'TPR_', type,'.png'), width = 800, height = 600)
plot_metric_comparison(full_results_Lik, n_sample, 'TPR')
dev.off()

png(paste0(Lik_results_folder, '/comparison_', "kl_divergence_", type,'.png'), width = 800, height = 600)
plot_metric_comparison(full_results_Lik, n_sample, "kl_divergence")
dev.off()

png(paste0(Lik_results_folder, '/comparison_', "jaccard_", type ,'.png'), width = 800, height = 600)
plot_metric_comparison(full_results_Lik, n_sample, "jaccard" )
dev.off()

png(paste0(Lik_results_folder, '/comparison_', "hamming_", type, '.png'), width = 800, height = 600)
plot_metric_comparison(full_results_Lik, n_sample, "hamming")
dev.off()

png(paste0(Lik_results_folder, '/comparison_', 'TPR_', type,'.png'), width = 800, height = 600)
plot_metric_comparison(full_results_Lik, n_sample, 'TPR')
dev.off()

png(paste0(Lik_results_folder, '/comparison_', 'TNR_', type,'.png'), width = 800, height = 600)
plot_metric_comparison(full_results_Lik, n_sample, 'TNR')
dev.off()

png(paste0(Lik_results_folder, '/comparison_', 'TPR_hub_', type,'.png'), width = 800, height = 600)
plot_metric_comparison(full_results_BIC, n_sample, 'TPR_hub')  
dev.off()

png(paste0(Lik_results_folder, '/comparison_', 'TNR_hub_', type,'.png'), width = 800, height = 600)
plot_metric_comparison(full_results_BIC, n_sample, 'TNR_hub')  
dev.off()

png(paste0(Lik_results_folder, '/comparison_', 'TPR_nonhub_', type,'.png'), width = 800, height = 600)
plot_metric_comparison(full_results_BIC, n_sample, 'TPR_nhub')  
dev.off()

png(paste0(Lik_results_folder, '/comparison_', 'TNR_nonhub_', type,'.png'), width = 800, height = 600)
plot_metric_comparison(full_results_BIC, n_sample, 'TNR_nhub')  
dev.off()


# png(paste0(Lik_results_folder, '/comparison_', 'hub_accuracy_', type,'.png'), width = 800, height = 600)
# plot_metric_comparison(full_results_BIC, n_sample, 'hub_accuracy')  
# dev.off()
# 
# png(paste0(Lik_results_folder, '/comparison_', 'nonhub_accuracy_', type,'.png'), width = 800, height = 600)
# plot_metric_comparison(full_results_BIC, n_sample, 'nonhub_accuracy')  
# dev.off()



save(full_results_Lik, file = paste0(Lik_results_folder,"/full_results_", type, ".RData"))


# }

