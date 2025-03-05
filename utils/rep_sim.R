run_simulations = function(experiment_name, n_sample, n_sim, rep_num, Theta, args, lambda, lambda_2step, threshold, folder, save){
  full_results_BIC = list()
  full_results_Lik = list()
  lambda_small = 10^seq(-8, 1, length = 50)
for(j in 1:n_sim){
  folder = paste0(experiment_name,'/sample_', n_sample[j])
  if (!dir.exists(folder)) {
    dir.create(folder)
  }
  for(rep in 1:rep_num){
    cat("\n\n\n----------------------Simulation ", rep,"-------sample---", n_sample[j], "-------------------------------\n\n\n")
    if(rep == 1){
      save = TRUE
    }else{
      save = FALSE
    }
    # save simulations results
      args$seed = args$seed + 1
      args$n = n_sample[j]
      if(rep == 1){
        #penHubs  
        cat("\n Working on penHubs\n")
        full_results_BIC$penHubs[[j]] = simulation_BIC_penHubs(Theta, args, lambda_small, threshold, folder, save=save)$results
        full_results_Lik$penHubs[[j]] = simulation_Lik_penHubs(Theta, args, lambda_small, threshold, folder, save=save)$results
        
        
        ##Real oracle
        full_results_BIC$Full_oracle[[j]] = simulation_BIC_Full_oracle(Theta, args, lambda, threshold, folder, save=save)$results
        full_results_Lik$Full_oracle[[j]] = simulation_Lik_Full_oracle(Theta, args, lambda, threshold, folder, save=save)$results
  
        cat('\n Lin Adj \n')
        full_results_BIC$lin_adj[[j]] = simulation_BIC_lin_adj(Theta, args, lambda, threshold, folder, save=save)$results
        full_results_Lik$lin_adj[[j]] = simulation_Lik_lin_adj(Theta, args, lambda, threshold, folder, save = save)$results
  
        #HWGLASSO
        cat("\n Working on HWGLASSO\n")
        full_results_BIC$hw[[j]] = simulation_BIC_hw(Theta, args, lambda, threshold, folder, save=save)$results
        full_results_Lik$hw[[j]] = simulation_Lik_hw(Theta, args, lambda, threshold, folder, save=save)$results
        
        #Adaptive GLASSO
        cat("\n working on Adapative GLASSO\n")
        full_results_BIC$adaGL[[j]] = simulation_BIC_adaGL(Theta, args, lambda, threshold, folder, save=save)$results
        full_results_Lik$adaGL[[j]] = simulation_Lik_adaGL(Theta, args, lambda, threshold, folder, save=save)$results
        
        #GLASSO no weight
        cat("\n Working on pure GLASSO\n")
        full_results_BIC$GL[[j]] = simulation_BIC_GL(Theta, args, lambda, threshold, folder, save=save)$results
        full_results_Lik$GL[[j]] = simulation_Lik_GL(Theta, args, lambda, threshold, folder, save=save)$results
        
        #Glasso oracle
        cat("\n Working on oracle\n")
        full_results_BIC$hub_oracle[[j]] = simulation_BIC_oracle(Theta, args, lambda, threshold, folder, save=save)$results
        full_results_Lik$hub_oracle[[j]] = simulation_Lik_oracle(Theta, args, lambda, threshold, folder, save=save)$results
        
        #two step
        cat("Working on 2 step\n")
        full_results_BIC$two_step[[j]] = simulation_BIC_2step(Theta, args, lambda_2step, threshold, folder, save=save)$results
        full_results_Lik$two_step[[j]] = simulation_Lik_2step(Theta, args, lambda_2step, threshold, folder,save=save)$results  
        
        cat("Working on New weights\n")
        for(func in names(weights_fun)){
          full_results_BIC[[func]][[j]] = simulation_BIC_ipchd(Theta, args,  weights_fun[[func]] ,lambda, threshold, folder, name = func, save=save)$results
          full_results_Lik[[func]][[j]] = simulation_Lik_ipchd(Theta, args, weights_fun[[func]] ,lambda, threshold, folder, name = func, save=save)$results
          #full_results_BIC[[paste0(func,"_frob")]][[j]] = simulation_BIC_ipchd_frob(Theta, args, weights_fun[[func]] ,lambda, threshold, folder, name = func)$results
        }
      }else{
        ##Real oracle
        full_results_BIC$Full_oracle[[j]] = sum_list(simulation_BIC_Full_oracle(Theta, args, lambda, threshold, folder, save=save)$results, full_results_BIC$Full_oracle[[j]])
        full_results_Lik$Full_oracle[[j]] = sum_list(simulation_Lik_Full_oracle(Theta, args, lambda, threshold, folder, save=save)$results, full_results_Lik$Full_oracle[[j]])
        
        cat('\n Lin Adj \n')
        full_results_BIC$lin_adj[[j]] = sum_list(simulation_BIC_lin_adj(Theta, args, lambda, threshold, folder)$results, full_results_BIC$lin_adj[[j]])
        full_results_Lik$lin_adj[[j]] = sum_list(simulation_Lik_lin_adj(Theta, args, lambda, threshold, folder)$results, full_results_Lik$lin_adj[[j]])
        
        #HWGLASSO
        cat("\n Working on HWGLASSO\n")
        full_results_BIC$hw[[j]] = sum_list(simulation_BIC_hw(Theta, args, lambda, threshold, folder, save=save)$results, full_results_BIC$hw[[j]])
        full_results_Lik$hw[[j]] = sum_list(simulation_Lik_hw(Theta, args, lambda, threshold, folder, save=save)$results, full_results_Lik$hw[[j]])
        
        #Adaptive GLASSO
        cat("\n working on Adapative GLASSO\n")
        full_results_BIC$adaGL[[j]] = sum_list(simulation_BIC_adaGL(Theta, args, lambda, threshold, folder, save=save)$results, full_results_BIC$adaGL[[j]])
        full_results_Lik$adaGL[[j]] = sum_list(simulation_Lik_adaGL(Theta, args, lambda, threshold, folder, save=save)$results, full_results_Lik$adaGL[[j]])
        
        #GLASSO no weight
        cat("\n Working on pure GLASSO\n")
        full_results_BIC$GL[[j]] = sum_list(simulation_BIC_GL(Theta, args, lambda, threshold, folder, save=save)$results, full_results_BIC$GL[[j]])
        full_results_Lik$GL[[j]] = sum_list(simulation_Lik_GL(Theta, args, lambda, threshold, folder, save=save)$results, full_results_Lik$GL[[j]])
        
        #Glasso oracle
        cat("\n Working on oracle\n")
        full_results_BIC$hub_oracle[[j]] = sum_list(simulation_BIC_oracle(Theta, args, lambda, threshold, folder, save=save)$results, full_results_BIC$hub_oracle[[j]])
        full_results_Lik$hub_oracle[[j]] = sum_list(simulation_Lik_oracle(Theta, args, lambda, threshold, folder, save=save)$results, full_results_Lik$hub_oracle[[j]])
        
        full_results_BIC$two_step[[j]] = sum_list(simulation_BIC_2step(Theta, args, lambda, threshold, folder, save=save)$results, full_results_BIC$two_step[[j]]) 
        full_results_Lik$two_step[[j]] = sum_list(simulation_Lik_2step(Theta, args, lambda, threshold, folder, save=save)$results, full_results_Lik$two_step[[j]])
        
        full_results_BIC$penHubs[[j]] = sum_list(simulation_BIC_penHubs(Theta, args, lambda_small, threshold, folder, save=save)$results, full_results_BIC$penHubs[[j]])
        full_results_Lik$penHubs[[j]] = sum_list(simulation_Lik_penHubs(Theta, args, lambda_small, threshold, folder, save=save)$results, full_results_Lik$penHubs[[j]])
        
        cat("Working on New weights\n")
        for(func in names(weights_fun)){
          full_results_BIC[[func]][[j]] = sum_list(simulation_BIC_ipchd(Theta, args,  weights_fun[[func]] ,lambda, threshold, folder, name = func, save=save)$results, full_results_BIC[[func]][[j]])
          full_results_Lik[[func]][[j]] = sum_list(simulation_Lik_ipchd(Theta, args, weights_fun[[func]] ,lambda, threshold, folder, name = func, save=save)$results, full_results_Lik[[func]][[j]])
          #full_results_BIC[[paste0(func,"_frob")]][[j]] = simulation_BIC_ipchd_frob(Theta, args, weights_fun[[func]] ,lambda, threshold, folder, name = func)$results
        }
      }
  }
  ##Real oracle
  full_results_BIC$Full_oracle[[j]] = sum_div(full_results_BIC$Full_oracle[[j]], rep_num)
  full_results_Lik$Full_oracle[[j]] = sum_div(full_results_Lik$Full_oracle[[j]], rep_num)
  
  cat('\n Lin Adj \n')
  full_results_BIC$lin_adj[[j]] = sum_div(full_results_BIC$lin_adj[[j]], rep_num)
  full_results_Lik$lin_adj[[j]] = sum_div(full_results_Lik$lin_adj[[j]], rep_num)
  
  #HWGLASSO
  cat("\n Working on HWGLASSO\n")
  full_results_BIC$hw[[j]] = sum_div(full_results_BIC$hw[[j]], rep_num)
  full_results_Lik$hw[[j]] = sum_div(full_results_Lik$hw[[j]], rep_num)
  
  #Adaptive GLASSO
  cat("\n working on Adapative GLASSO\n")
  full_results_BIC$adaGL[[j]] = sum_div(full_results_BIC$adaGL[[j]], rep_num)
  full_results_Lik$adaGL[[j]] = sum_div(full_results_Lik$adaGL[[j]], rep_num)
  
  #GLASSO no weight
  cat("\n Working on pure GLASSO\n")
  full_results_BIC$GL[[j]] = sum_div(full_results_BIC$GL[[j]], rep_num)
  full_results_Lik$GL[[j]] = sum_div( full_results_Lik$GL[[j]], rep_num)
  
  #Glasso oracle
  cat("\n Working on oracle\n")
  full_results_BIC$hub_oracle[[j]] = sum_div(full_results_BIC$hub_oracle[[j]], rep_num)
  full_results_Lik$hub_oracle[[j]] = sum_div(full_results_Lik$hub_oracle[[j]], rep_num)
  
  #two step
  cat("Working on 2 step\n")
  full_results_BIC$two_step[[j]] = sum_div(full_results_BIC$two_step[[j]], rep_num)
  full_results_Lik$two_step[[j]] = sum_div(full_results_Lik$two_step[[j]], rep_num)
    
  cat("\nWorking on penHubs\n")  
  full_results_BIC$penHubs[[j]] = sum_div(full_results_BIC$penHubs[[j]], rep_num)
  full_results_Lik$penHubs[[j]] = sum_div(full_results_Lik$penHubs[[j]], rep_num)
  
  
  cat("Working on New weights\n")
  for(func in names(weights_fun)){
    full_results_BIC[[func]][[j]] = sum_div(full_results_BIC[[func]][[j]], rep_num)
    full_results_Lik[[func]][[j]] = sum_div(full_results_Lik[[func]][[j]], rep_num)
    #full_results_BIC[[paste0(func,"_frob")]][[j]] = simulation_BIC_ipchd_frob(Theta, args, weights_fun[[func]] ,lambda, threshold, folder, name = func)$results
  }
}
  
  return(list(BIC = full_results_BIC, Lik = full_results_Lik))
}