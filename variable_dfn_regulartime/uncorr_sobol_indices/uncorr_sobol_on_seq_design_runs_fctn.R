###############################################################################
### Analysis of emulator built on sequential design draws from dfnWorks
### Using Global SA: A Primer by A. Saltelli & Saltelli 2010
### This script calculates the 'uncorrelated' sobol indices for the functional 
### variables.
## Author: Alexander C. Murph
## Date: January 2024
library(hetGP)
library(ggplot2)
library(reshape2)
library(gridExtra)
library(latex2exp)
setwd("~/GitLab/sa_for_chemistry/variable_dfn_regulartime/uncorr_sobol_indices")
source("MC_sobol_indices.R")

# Saltelli suggests at least 500 samples:
num_samples     = 750
n_bootstraps    = 300

sim_num = commandArgs()
print(sim_num)
print(sim_num[6])
sim_num = as.integer(sim_num[6])
set.seed(sim_num)
boot_idx       = sim_num
input_var_idx  = (sim_num-1) %% 7 + 1
sim_num        = floor((sim_num-1)/7) + 1
output_var_idx = sim_num + 20

# Grab the final hetGP model:
models_loc = "~/GitLab/sa_for_chemistry/variable_dfn/joint_emulation/ftcn_models"

#########################
## Load in all the data and information on which vars we are analyzing.
# Load in the full training data:
training_data   = read.csv(file="../joint_emulation/training_data.csv") 
training_data$X = NULL
# Load in the testing data:
testing_data   = read.csv(file="../joint_emulation/testing_data.csv") 
testing_data$X = NULL
# Load in all data:
full_data      = read.csv(file="../joint_emulation/full_data.csv") 
full_data$X    = NULL

# Variable info:
new_names = c( "dfn_seed",                   "curr_inflow_rate",           "curr_diffusion_coef",        
               "curr_porosity",              "curr_gypsum_rate_constant",
               "curr_calcite_rate_constant", "curr_gypsum_surface_area",   "curr_calcite_surface_area",  
               "dfn_p32",                    "dfn_volume",                
               "backbone_volume",            "backbone_p32",               "final_gypsum",               
               "final_calcite",              "Pe",                        
               "Da_1_gypsum",                "Da_2_gypsum",                "Da_1_calcite",               
               "Da_2_calcite",               "tau",                       
               "calcite_flush_0_nondim",     "gypsum_flush_0_nondim",      "calcite_flush_1_nondim",     
               "gypsum_flush_1_nondim",      "calcite_flush_2_nondim",    
               "gypsum_flush_2_nondim",      "calcite_flush_3_nondim",     "gypsum_flush_3_nondim",      
               "calcite_flush_4_nondim",     "gypsum_flush_4_nondim",     
               "calcite_flush_5_nondim",     "gypsum_flush_5_nondim",      "calcite_flush_6_nondim",     
               "gypsum_flush_6_nondim",      "calcite_flush_7_nondim",    
               "gypsum_flush_7_nondim",      "calcite_flush_8_nondim",     "gypsum_flush_8_nondim",      
               "calcite_flush_9_nondim",     "gypsum_flush_9_nondim",     
               "calcite_flush_10_nondim",    "gypsum_flush_10_nondim",     
               "calcite_flush_0",            "gypsum_flush_0",             "calcite_flush_1",           
               "gypsum_flush_1",             "calcite_flush_2",            "gypsum_flush_2",             
               "calcite_flush_3",            "gypsum_flush_3",            
               "calcite_flush_4",            "gypsum_flush_4",             "calcite_flush_5",            
               "gypsum_flush_5",             "calcite_flush_6",           
               "gypsum_flush_6",             "max_calcite",                "min_gypsum" )

# Variable updates after the EDA.
taopoints_to_remove = c("calcite_flush_0_nondim",     "gypsum_flush_0_nondim",      "calcite_flush_1_nondim",     
                        "gypsum_flush_1_nondim",      "calcite_flush_2_nondim",    
                        "gypsum_flush_2_nondim",      "calcite_flush_3_nondim",     "gypsum_flush_3_nondim",      
                        "calcite_flush_4_nondim",     "gypsum_flush_4_nondim",     
                        "calcite_flush_5_nondim",     "gypsum_flush_5_nondim",      "calcite_flush_6_nondim",     
                        "gypsum_flush_6_nondim",      "calcite_flush_7_nondim",    
                        "gypsum_flush_7_nondim",      "calcite_flush_8_nondim",     "gypsum_flush_8_nondim",      
                        "calcite_flush_9_nondim",     "gypsum_flush_9_nondim",     
                        "calcite_flush_10_nondim",    "gypsum_flush_10_nondim",
                        "gypsum_flush_0", 
                        "gypsum_flush_1", "gypsum_flush_2",             
                        "gypsum_flush_3",            
                        "gypsum_flush_4",
                        "gypsum_flush_5",
                        "gypsum_flush_6",             
                        "max_calcite",                
                        "min_gypsum")

new_names           = new_names[which(!(new_names%in%taopoints_to_remove))]


cols_of_log_vars                 = c(2,3,5,6,8)
cols_of_net_vars                 = c(9,10,11,12)
cols_of_non_fctn_vars            = c(2:14)
cols_of_fctn_vars                = c(15,16,17,18,19,20)
output_vars                      = c(21:length(new_names))
input_vars_for_emulation_fctn    = c(cols_of_net_vars[2], cols_of_fctn_vars)
input_var_idx_raw                = input_var_idx
input_var_idx                    = input_vars_for_emulation_fctn[input_var_idx]
input_vars_for_emulation_nonfctn = c(2:8, 10)

########################
## Lets bootstrap this boi!
tuning_2 = c("calcite_flush_2")
tuning_3 = c("calcite_flush_0", "calcite_flush_1")
tuning_4 = c('gypsum_flush_10_nondim')

input_vars_names = new_names[input_vars_for_emulation_fctn]
n_obs            = nrow(full_data)
input_dimension  = length(input_vars_for_emulation_fctn)

model_name  = paste("../joint_emulation/ftcn_models/het_gp_output_",
                    new_names[output_var_idx], ".Rdata", sep = "")
model       = NULL
load(model_name)

orig_model          = model
bootstrapped_FOmean = NULL
bootstrapped_TEmean = NULL
bootstrapped_FOsd2  = NULL
bootstrapped_TEsd2  = NULL
bootstrapped_TEseed = c()

########################
## For this script, I need to uncorrelated the variables with a regression-
## based graham-schmidt process. 
other_input_vars = input_vars_for_emulation_fctn[which(input_vars_for_emulation_fctn!=input_var_idx)]
temp_full_data   = full_data
for(input_idx in 1:length(other_input_vars)){
  temp_data      = as.data.frame(temp_full_data[,new_names[input_var_idx]])
  temp_data$out  = temp_full_data[,other_input_vars[input_idx]]
  temp_lm        = lm(out ~ ., data = temp_data)
  cond_mean      = predict(temp_lm, x = as.data.frame(temp_full_data[,new_names[input_var_idx]]))
  
  # Remove correlated effect of variables so far from this input var.
  temp_full_data[,other_input_vars[input_idx]] = temp_full_data[,other_input_vars[input_idx]] - cond_mean
}
# Finish by removing the effect for the correlated input of interest:
temp_data     = as.data.frame(temp_full_data[,other_input_vars])
temp_data$out = temp_full_data[,input_var_idx]
temp_lm       = lm(out ~ ., data = temp_data)
cond_mean     = predict(temp_lm, x = as.data.frame(temp_full_data[,other_input_vars]))

# Remove correlated effect of variables so far from this input var.
temp_full_data[,input_var_idx] = temp_full_data[,input_var_idx] - cond_mean

# Set full data variables to appropriate data for this experiment.
full_data_X = temp_full_data[,input_vars_for_emulation_fctn]
full_data_Z = temp_full_data[,output_var_idx]

model               = NULL
bootstrapped_FOmean = NULL
bootstrapped_TEmean = NULL
bootstrapped_FOsd2  = NULL
bootstrapped_TEsd2  = NULL
bootstrapped_TEseed = c()

print(paste("Starting bootstrap iteration:", boot_idx))
resample_idxs = sample(1:n_obs, size = n_obs, replace = TRUE)
resampled_X   = full_data_X[resample_idxs,]
resampled_Z   = full_data_Z[resample_idxs]

complete_resampled = which(!is.na(resampled_Z))
var_of_response    = var(resampled_Z[complete_resampled])

if(new_names[output_var_idx] %in% tuning_2){
  model = mleHetGP( X = list(X0 = as.matrix(resampled_X[complete_resampled,]), 
                             Z0 = as.vector(resampled_Z[complete_resampled]), 
                             mult = rep(1,nrow(as.matrix(resampled_X[complete_resampled,])))),
                    noiseControl = list(g_min = 3),
                    lower = 1,
                    upper = 200,
                    Z = as.vector(resampled_Z[complete_resampled]),
                    settings = list(checkHom = FALSE ) )
} else if(new_names[output_var_idx] %in% tuning_3){
  model = mleHetGP( X = list(X0 = as.matrix(resampled_X[complete_resampled,]), 
                             Z0 = as.vector(resampled_Z[complete_resampled]), 
                             mult = rep(1,nrow(as.matrix(resampled_X[complete_resampled,])))),
                    noiseControl = list(g_min = 3),
                    lower = 0.5,
                    upper = 200,
                    Z = as.vector(resampled_Z[complete_resampled]),
                    settings = list(checkHom = FALSE ) )
} else if(new_names[output_var_idx] %in% tuning_4){
  model = mleHetGP( X = list(X0 = as.matrix(resampled_X[complete_resampled,]), 
                             Z0 = as.vector(resampled_Z[complete_resampled]), 
                             mult = rep(1,nrow(as.matrix(resampled_X[complete_resampled,])))),
                    noiseControl = list(g_min = 3),
                    Z = as.vector(resampled_Z[complete_resampled]), covtype = "Matern5_2",
                    settings = list(checkHom = FALSE ) )
} else {
  model = mleHetGP( X = list(X0 = as.matrix(resampled_X[complete_resampled,]), 
                             Z0 = as.vector(resampled_Z[complete_resampled]), 
                             mult = rep(1,nrow(as.matrix(resampled_X[complete_resampled,])))),
                    noiseControl = list(g_min = 3),
                    Z = as.vector(resampled_Z[complete_resampled]), 
                    covtype = "Matern3_2",
                    settings = list(checkHom = FALSE ) )
}
model               = rebuild(model, robust = TRUE)
temp_SAs            = estimate_sobol_index(model, num_samples, input_dimension, input_var_idx_raw,
                                           var_of_response, var_names = input_vars_names)

bootstrapped_FOmean = as.data.frame(t(as.matrix(unlist(unlist(temp_SAs[['first_order_mean']])))))
bootstrapped_TEmean = as.data.frame(t(as.matrix(unlist(unlist(temp_SAs[['total_effect_mean']])))))
bootstrapped_FOsd2  = as.data.frame(t(as.matrix(unlist(unlist(temp_SAs[['first_order_sd2']])))))
bootstrapped_TEsd2  = as.data.frame(t(as.matrix(unlist(unlist(temp_SAs[['total_effect_sd2']])))))

############################
# Get the predictions from the dispersion GP
predict_values     = predict(object=model, x = as.matrix(resampled_X))
dispersion_GP_vals = predict_values$nugs
mean_GP_vals       = predict_values$mean
avg_dis_GP_vals    = mean(dispersion_GP_vals)
var_mean_GP_vals   = var(mean_GP_vals)

# Predict the total effect:
# bootstrapped_TEseed               = c(bootstrapped_TEseed, avg_dis_GP_vals / (var_of_response))
bootstrapped_TEseed_df            = data.frame(TEseed = avg_dis_GP_vals / (var_of_response))
bootstrapped_TEseed_df$output_var = rep(new_names[output_var_idx], times = 1)
bootstrapped_TEseed_df$input_var  = rep(new_names[input_var_idx],  times = 1)
bootstrapped_TEseed_df$sobol_type = rep('TEseed', times = 1)
bootstrapped_TEseed_df$boot_num   = boot_idx

write.csv(bootstrapped_TEseed_df,  file = paste("bootstrap_data_fctn/", 
                                                "seed_variable__", new_names[output_var_idx], 
                                                "__", new_names[input_var_idx], 
                                                "__", boot_idx, ".csv", sep = "") )

# Now put together the full data for this output variable.
bootstrapped_FOmean$output_var = rep(new_names[output_var_idx], times = 1)
bootstrapped_FOmean$input_var  = rep(new_names[input_var_idx], times = 1)
bootstrapped_FOmean$sobol_type = rep('FOmean', times = 1)
bootstrapped_FOmean$boot_num   = boot_idx
bootstrapped_TEmean$output_var = rep(new_names[output_var_idx], times = 1)
bootstrapped_TEmean$input_var  = rep(new_names[input_var_idx], times = 1)
bootstrapped_TEmean$sobol_type = rep('TEmean', times = 1)
bootstrapped_TEmean$boot_num   = boot_idx
bootstrapped_FOsd2$output_var  = rep(new_names[output_var_idx], times = 1)
bootstrapped_FOsd2$input_var   = rep(new_names[input_var_idx], times = 1)
bootstrapped_FOsd2$sobol_type  = rep('FOsd2', times = 1)
bootstrapped_FOsd2$boot_num    = boot_idx
bootstrapped_TEsd2$output_var  = rep(new_names[output_var_idx], times = 1)
bootstrapped_TEsd2$input_var   = rep(new_names[input_var_idx], times = 1)
bootstrapped_TEsd2$sobol_type  = rep('TEsd2', times = 1)
bootstrapped_TEsd2$boot_num    = boot_idx

full_bootstrap_data = rbind(bootstrapped_FOmean, bootstrapped_TEmean,
                            bootstrapped_FOsd2,  bootstrapped_TEsd2)
write.csv(full_bootstrap_data, file = paste('bootstrap_data_fctn/', new_names[output_var_idx], 
                                            "__", new_names[input_var_idx], 
                                            "__", boot_idx, ".csv", sep = ''))


