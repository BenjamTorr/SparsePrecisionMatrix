################################################################################

rm(list=ls())
setwd("D:/UNC/PhD UNC/Research/Sparse Precision Matrix/code")

invisible(lapply(list.files("JoseMethods", pattern = "\\.R$", full.names = TRUE), source))
invisible(lapply(list.files("utils", pattern = "\\.R$", full.names = TRUE), source))
invisible(lapply(list.files("HWGLASSO", pattern = "\\.R$", full.names = TRUE), source))
invisible(lapply(list.files('pretuning',pattern = "\\.R$", full.names = TRUE), source))

source('SimulationParameters.R')
##############################################################################s##
library(ggplot2)
library(reshape2)
#---------------------------Generate data---------------------------------------#

set.seed(21)

Theta = r.sparse.pdhubmat(
    args$p, T0 = args$T0, r = args$r, ph = args$ph, pnh = args$pnh, pneff = args$pneff, diagonal_shift = args$diagonal_shift, 
    shuffle = args$shuffle, type = "unif", hmin = args$hmin, hmax = args$hmax, nhmin = args$nhmin, nhmax = args$nhmax, 
    neffmin = args$neffmin, neffmax = args$neffmax)

sim_data = sample_from_precision(args$n, args$p, Theta)
X = sim_data$X
X_train = sim_data$X_train
X_val = sim_data$X_val

CoVR_X = get_cov_cor(X)
CoVR_X_train = get_cov_cor(X_train)
CoVR_X_val = get_cov_cor(X_val)
#--------------------------- Weight generation----------------------------------#

icphd_inf = sta_ipchd(X, mat_type = 'cor', mat = CoVR_X$Corr, var_inds = 1:args$p, overest_type = 'frac')
w_icphd = inverse_prop(args$p, icphd_inf)

hw = hw_weights(CoVR_X$Cov)
hw_inf = hw$influence
w_hw = hw$w

#training and validation for likelihood method
icphd_inf_tr = sta_ipchd(X, mat_type = 'cor', mat = CoVR_X_train$Corr, var_inds = 1:args$p, overest_type = 'frac')
w_icphd_tr = inverse_prop(args$p, icphd_inf_tr)

hw_train = hw_weights(CoVR_X_train$Cov)
hw_inf_train = hw_train$influence
w_hw_train = hw_train$w

#---------------------------Precision Matrix recovery---------------------------------#

lambda = 10^seq(-5, -2, length = 500)

Theta_icphd_BIC = BICweighted_glasso(mat = CoVR_X$Cov, rho = lambda, .W = w_icphd, p = args$p, n = args$n,
                                     threshold = 0, maxit = 200,
                                     penalize.diagonal = TRUE)

Theta_hw_BIC = BICweighted_glasso(mat = CoVR_X$Cov, rho = lambda, .W = w_hw, p = args$p, n = args$n,
                                  threshold = 0, maxit = 200,
                                  penalize.diagonal = TRUE)
## Check what happens with likelihood
Theta_icphd_Lik = Lik_weighted_glasso(mat_train = CoVR_X_train$Cov, mat_val = CoVR_X_val$Cov, rho = lambda, 
                                      .W = w_icphd_tr , p = args$p, n = args$n,
                                       threshold = 0, maxit = 200,
                                       penalize.diagonal = TRUE)

Theta_hw_Lik = Lik_weighted_glasso(mat_train = CoVR_X_train$Cov, mat_val = CoVR_X_val$Cov, rho = lambda, 
                                      .W = w_hw_train , p = args$p, n = args$n,
                                      threshold = 0, maxit = 200,
                                      penalize.diagonal = TRUE)

#-----------------------------Plots-----------------------------------------------#
plot(eigen(Theta)$values, main = 'Eigenvalues of Precision matrix', xlab = 'Variable', ylab = 'Eigenvalue')
plot(eigen(.PMtoIC(Theta))$values, main = 'Eigenvalues of inverse correlation matrix', xlab = 'Variable', ylab = 'Eigenvalue')

plot(icphd_inf, main = 'Influence measure of correlation', ylab = 'ICPHD Measure', xlab = 'Variable')
plot(hw_inf, main = 'L1 norm of rows', ylab = 'L1', xlab = 'Variable')

plot(w_icphd[1, 1:args$p], main = 'Weights for IPCHD', xlab = 'Variable', ylab = 'Weight of of entry 1,j', 
     ylim = c(min(w_icphd[50, 1:args$p]), quantile(w_icphd[50, 1:args$p], 0.95)))
points(w_icphd[50, 1:args$p], col = 'blue')
legend("topright", legend = c("Hub node - 1", "Non-hub node 50"), col = c("black", "blue"), pch = 19)

plot(w_hw[1, 1:args$p], main = 'Weights for HWGLASSO', xlab = 'Variable', ylab = 'Weight of of entry 1,j', 
     ylim = c(min(w_hw[50, 1:args$p]), quantile(w_hw[50, 1:args$p], 0.95)))
points(w_hw[50, 1:args$p], col = 'blue')
legend("topright", legend = c("Hub node - 1", "Non-hub node 50"), col = c("black", "blue"), pch = 19)

## Sanity Check for GLasso fitthing methods
plot(x = lambda, y = Theta_icphd_BIC$BIC, log = 'x', main ='BIC for IPCHD', ylab = 'BIC', xlab = 'Lambda')
plot(x = lambda, y = Theta_hw_BIC$BIC, log = 'x', main ='BIC for HWgLASSO', ylab = 'BIC', xlab = 'Lambda')
plot(x = lambda, y = Theta_icphd_Lik$BIC, log = 'x', main ='Negative Likelihood for IPCHD', ylab = 'Negative Likelihood', xlab = 'Lambda')
plot(x = lambda, y = Theta_hw_Lik$BIC, log = 'x', main ='Negative Likelihood for HWGLASSO', ylab = 'Negative Likelihood', xlab = 'Lambda')

raw_results_icphd_BIC = all_comparisons(Theta, Theta_icphd_BIC$optimal.model$wi, threshold = 1e-05)

plot.pdm(Theta, main = 'Real Precision Matrix')

plot.pdm(Theta - Theta_icphd_BIC$optimal.model$wi, main = 'Difference')

plot.pdm(Theta_icphd_BIC$optimal.model$wi, main = 'Estimated IPCHD refined')

plot.pdm(Theta_hw_BIC$optimal.model$wi, main = 'Estimated HWGLASSO')
