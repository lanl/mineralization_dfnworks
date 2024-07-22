###############################################################################
### Analysis of emulator built on sequential design draws from dfnWorks
### Using Global SA: A Primer by A. Saltelli & Saltelli 2010
## Author: Alexander C. Murph
## Date: September 2023
library(hetGP)
library(ggplot2)
library(reshape2)
library(gridExtra)
library(latex2exp)
setwd("~/GitLab/sa_for_chemistry/fixed_dfn/sobol_indices_fixedDFN")
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
sim_num        = (sim_num-1) %% 6 + 1
# output_vars    = c(20,22,24,26,28,30,32,34) + 1


# Grab the final hetGP model:
models_loc = "~/GitLab/sa_for_chemistry/fixed_dfn/joint_emulation/non_ftcn_models"

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
new_names = c("dfn_seed",                   "curr_inflow_rate",           "curr_diffusion_coef",        "curr_porosity",              "curr_gypsum_rate_constant", 
              "curr_calcite_rate_constant", "curr_gypsum_surface_area",   "curr_calcite_surface_area",  "dfn_p32",                    "dfn_volume",                
              "backbone_volume",            "backbone_p32",               "final_gypsum",               "final_calcite",              "Pe",                        
              "Da_1_gypsum",                "Da_2_gypsum",                "Da_1_calcite",               "Da_2_calcite",               "tau",                       
              "calcite_flush_0_nondim",     "gypsum_flush_0_nondim",      "calcite_flush_1_nondim",     "gypsum_flush_1_nondim",      "calcite_flush_2_nondim",    
              "gypsum_flush_2_nondim",      "calcite_flush_3_nondim",     "gypsum_flush_3_nondim",      "calcite_flush_4_nondim",     "gypsum_flush_4_nondim",     
              "calcite_flush_5_nondim",     "gypsum_flush_5_nondim",      "calcite_flush_6_nondim",     "gypsum_flush_6_nondim",      "calcite_flush_7_nondim",    
              "gypsum_flush_7_nondim",      "calcite_flush_8_nondim",     "gypsum_flush_8_nondim",      "calcite_flush_9_nondim",     "gypsum_flush_9_nondim",     
              "calcite_flush_10_nondim",    "gypsum_flush_10_nondim",     "calcite_flush_11_nondim",    "gypsum_flush_11_nondim",     "calcite_flush_12_nondim",   
              "gypsum_flush_12_nondim")

gypsum_outputs  = c("gypsum_flush_1_nondim", "gypsum_flush_1_nondim", "gypsum_flush_2_nondim",
                    "gypsum_flush_3_nondim", "gypsum_flush_4_nondim",     
                    "gypsum_flush_5_nondim", "gypsum_flush_6_nondim", 
                    "gypsum_flush_7_nondim", "gypsum_flush_8_nondim", "gypsum_flush_9_nondim",     
                    "gypsum_flush_10_nondim", "gypsum_flush_11_nondim",   
                    "gypsum_flush_12_nondim")
outputs_abrv    = c("flush_0", "flush_1",
                    "flush_2", "flush_3", "flush_4",     
                    "flush_5", "flush_6", 
                    "flush_7", "flush_8", "flush_9",     
                    "flush_10", "flush_11",   
                    "flush_12")
calcite_outputs = c("calcite_flush_0_nondim",      "calcite_flush_1_nondim", 
                    "calcite_flush_2_nondim",    
                    "calcite_flush_3_nondim",     "calcite_flush_4_nondim",
                    "calcite_flush_5_nondim",     "calcite_flush_6_nondim",
                    "calcite_flush_7_nondim",    
                    "calcite_flush_8_nondim",     "calcite_flush_9_nondim",      
                    "calcite_flush_10_nondim",    "calcite_flush_11_nondim",  "calcite_flush_12_nondim")



# For a fixed DFN, this turned out to be slightly different:
# taopoints_to_remove = c("calcite_flush_0_nondim", "calcite_flush_1_nondim", 
#                         "calcite_flush_2_nondim", "calcite_flush_3_nondim",
#                         "gypsum_flush_0_nondim", "gypsum_flush_1_nondim",
#                         "gypsum_flush_2_nondim", "gypsum_flush_3_nondim")
taopoints_to_remove = c("calcite_flush_0_nondim", "calcite_flush_1_nondim",
                        "calcite_flush_2_nondim", "calcite_flush_3_nondim",
                        "gypsum_flush_0_nondim", "gypsum_flush_1_nondim",
                        "gypsum_flush_2_nondim", "gypsum_flush_3_nondim",
                        "gypsum_flush_4_nondim", "gypsum_flush_5_nondim",
                        "gypsum_flush_6_nondim", "gypsum_flush_7_nondim",
                        "gypsum_flush_8_nondim", "gypsum_flush_9_nondim",
                        "calcite_flush_10_nondim",    "gypsum_flush_10_nondim",     "calcite_flush_11_nondim",    "gypsum_flush_11_nondim",     "calcite_flush_12_nondim",   
                        "gypsum_flush_12_nondim")


# # From here, we rescale.  Then split and save the data.
# # Let's be more principled with this rescaling:
# calcite_data      = as.data.frame(norm_data[,calcite_outputs])
# final_calcite     = norm_data[,"final_calcite"]
# # This is often zero....why might calcite go down?
# 
# f = function(x){max(x,na.rm=T)}
# max_values = as.numeric(apply(calcite_data, 1, f))
# # max_values = final_calcite
# ff = function(x){
#   x / max_values
# }
# for(col_idx in 1:ncol(calcite_data)){
#   calcite_data[,col_idx] = calcite_data[,col_idx] / max_values
# }
# norm_data[,calcite_outputs] = calcite_data


## Based on the EDA, I am going to drop the initial calcite points and do some fixes to the gypsum data.
# Drop the first few calcite:
new_names = new_names[which(!(new_names%in%taopoints_to_remove))]
#norm_data = norm_data[new_names]

# Breaking the input vars into different sets, separating out the output vars.
cols_of_log_vars              = c(2,3,5,6,8)
cols_of_net_vars              = c(9,10,11,12)
cols_of_non_fctn_vars         = c(2:14)
# Pe and Tau are perfectly correlated.  I'm removing tau
cols_of_fctn_vars             = c(15,16,17,18,19)
# I'm going to focus on calcite here for a bit:
output_vars                   = c(21:length(new_names))
output_var_idx                = output_vars[sim_num]
# output_vars                   = c(21,22,24,26,28,30,32,34)
# For the fixed DFN, I no longer use any network variables
# input_vars_for_emulation_fctn = c(cols_of_net_vars[2], cols_of_fctn_vars)
input_vars_for_emulation_fctn = c(cols_of_fctn_vars)
input_vars_for_emulation_nonfctn = c(2:8)

########################
## Lets bootstrap this boi!
tuning_2 = c('gypsum_flush_0_nondim', "gypsum_flush_1_nondim",'calcite_flush_5_nondim',
             'calcite_flush_10_nondim',
             'calcite_flush_11_nondim',
             'calcite_flush_12_nondim')
tuning_3 = c('calcite_flush_1_nondim',
             'gypsum_flush_12_nondim', 
             'gypsum_flush_11_nondim')
tuning_4 = c('calcite_flush_2_nondim','calcite_flush_4_nondim', 'calcite_flush_8_nondim',
             'calcite_flush_9_nondim',  
             'gypsum_flush_10_nondim')
tuning_5 = c('calcite_flush_3_nondim',  
             'gypsum_flush_9_nondim',
             'gypsum_flush_4_nondim', "gypsum_flush_8_nondim", 
             'gypsum_flush_7_nondim','gypsum_flush_2_nondim', 'gypsum_flush_3_nondim', 
             'gypsum_flush_5_nondim', 'gypsum_flush_6_nondim')

input_vars_names = new_names[input_vars_for_emulation_nonfctn]
n_obs            = nrow(full_data)
input_dimension  = length(input_vars_for_emulation_nonfctn)

# Set full data variables to appropriate data for this experiment.
full_data_X = full_data[,input_vars_for_emulation_nonfctn]
full_data_Z = full_data[,output_var_idx]

model_name = paste("../joint_emulation/non_ftcn_models/het_gp_output_",
                   new_names[output_var_idx], ".Rdata", sep = "")
model      = NULL
load(model_name)

orig_model          = model
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
hetGP_inputs       = find_reps(X = as.matrix(resampled_X[complete_resampled,]), 
                               Z = as.vector(resampled_Z[complete_resampled]))

var_of_response = var(resampled_Z[complete_resampled])

if(new_names[output_var_idx] %in% tuning_2){
  model = mleHomGP( X = list(X0 = as.matrix(resampled_X[complete_resampled,]), 
                             Z0 = as.vector(resampled_Z[complete_resampled]), 
                             mult = rep(1,nrow(as.matrix(resampled_X[complete_resampled,])))),
                    # noiseControl = list(g_min = 3),
                    Z = as.vector(resampled_Z[complete_resampled]),
                    lower = 0.5,
                    upper = 100 )
} else if(new_names[output_var_idx] %in% tuning_3){
  model = mleHomGP( X = list(X0 = as.matrix(resampled_X[complete_resampled,]), 
                             Z0 = as.vector(resampled_Z[complete_resampled]), 
                             mult = rep(1,nrow(as.matrix(resampled_X[complete_resampled,])))),
                    # noiseControl = list(g_min = 3),
                    Z = as.vector(resampled_Z[complete_resampled]), covtype = "Matern3_2",
                    lower = 0.5,
                    upper = 100 )
} else if(new_names[output_var_idx] %in% tuning_4){
  model = mleHomGP( X = list(X0 = as.matrix(resampled_X[complete_resampled,]), 
                             Z0 = as.vector(resampled_Z[complete_resampled]), 
                             mult = rep(1,nrow(as.matrix(resampled_X[complete_resampled,])))),
                    # noiseControl = list(g_min = 3),
                    Z = as.vector(resampled_Z[complete_resampled]),
                    lower = 5,
                    upper = 110 )
}else if(new_names[output_var_idx] %in% tuning_5){
  model = mleHomGP( X = list(X0 = as.matrix(resampled_X[complete_resampled,]), 
                             Z0 = as.vector(resampled_Z[complete_resampled]), 
                             mult = rep(1,nrow(as.matrix(resampled_X[complete_resampled,])))),
                    lower = 15,
                    Z = as.vector(resampled_Z[complete_resampled]),
                    upper = 500 )
} else {
  
  model = mleHomGP( X = list(X0 = as.matrix(resampled_X[complete_resampled,]), 
                             Z0 = as.vector(resampled_Z[complete_resampled]), 
                             mult = rep(1,nrow(as.matrix(resampled_X[complete_resampled,])))),
                    lower = 0.5,
                    upper = 100,
                    Z = as.vector(resampled_Z[complete_resampled]), covtype = "Matern3_2" )
}

model               = rebuild(model, robust = TRUE)
temp_SAs            = estimate_sobol_indices(model, num_samples, input_dimension, var_of_response,
                                              var_names = input_vars_names)

bootstrapped_FOmean = as.data.frame(t(as.matrix(unlist(unlist(temp_SAs[['first_order_mean']])))))
# bootstrapped_FOmean = rbind(bootstrapped_FOmean, temp_row)

bootstrapped_TEmean = as.data.frame(t(as.matrix(unlist(unlist(temp_SAs[['total_effect_mean']])))))
# bootstrapped_TEmean = rbind(bootstrapped_TEmean, temp_row)

bootstrapped_FOsd2  = as.data.frame(t(as.matrix(unlist(unlist(temp_SAs[['first_order_sd2']])))))
# bootstrapped_FOsd2  = rbind(bootstrapped_FOsd2, temp_row)

bootstrapped_TEsd2  = as.data.frame(t(as.matrix(unlist(unlist(temp_SAs[['total_effect_sd2']])))))
# bootstrapped_TEsd2  = rbind(bootstrapped_TEsd2, temp_row)

############################
# Get the predictions from the dispersion GP
predict_values     = predict(object=model, x = as.matrix(resampled_X))
dispersion_GP_vals = predict_values$nugs
mean_GP_vals       = predict_values$mean
avg_dis_GP_vals    = mean(dispersion_GP_vals)
var_mean_GP_vals   = var(mean_GP_vals)

# Predict the total effect:
# bootstrapped_TEseed = c(bootstrapped_TEseed, avg_dis_GP_vals / (var_of_response))
  
bootstrapped_TEseed_df            = data.frame(TEseed = avg_dis_GP_vals / (var_of_response))
bootstrapped_TEseed_df$output_var = rep(new_names[output_var_idx], times = 1)
bootstrapped_TEseed_df$sobol_type = rep('TEseed', times = 1)
bootstrapped_TEseed_df$boot_num   = boot_idx

write.csv(bootstrapped_TEseed_df,  file = paste("bootstrap_data_nonfctn/", "seed_variable__", 
                                                new_names[output_var_idx], "__", boot_idx, ".csv", sep = "") )

# Now put together the full data for this output variable.
bootstrapped_FOmean$output_var = rep(new_names[output_var_idx], times = 1)
bootstrapped_FOmean$sobol_type = rep('FOmean', times = 1)
bootstrapped_FOmean$boot_num   = boot_idx
bootstrapped_TEmean$output_var = rep(new_names[output_var_idx], times = 1)
bootstrapped_TEmean$sobol_type = rep('TEmean', times = 1)
bootstrapped_TEmean$boot_num   = boot_idx
bootstrapped_FOsd2$output_var  = rep(new_names[output_var_idx], times = 1)
bootstrapped_FOsd2$sobol_type  = rep('FOsd2', times = 1)
bootstrapped_FOsd2$boot_num    = boot_idx
bootstrapped_TEsd2$output_var  = rep(new_names[output_var_idx], times = 1)
bootstrapped_TEsd2$sobol_type  = rep('TEsd2', times = 1)
bootstrapped_TEsd2$boot_num    = boot_idx

full_bootstrap_data = rbind(bootstrapped_FOmean,bootstrapped_TEmean,
                            bootstrapped_FOsd2, bootstrapped_TEsd2)
write.csv(full_bootstrap_data, file = paste('bootstrap_data_nonfctn/', 
                                            new_names[output_var_idx], "_", 
                                            boot_idx, '.csv', sep = ''))


