# Sequential design on input locations for dfnWorks study.
# This is an updated code suite where I update some issues I had created
# by using bkde unnecessarily.  For this file, I just use the quantile
# function with 'percentile.'
## Author: Alexander C. Murph
## Date: August 2023
library(hetGP)
library(parallel)
setwd("~/GitLab/sa_for_chemistry/fixed_dfn/joint_emulation")
source("sequential_helpers.R")

# This is not a simulation study, so I'll set the seed.
set.seed(13)

# Data location:
data_path = "~/GitLab/sa_for_chemistry/fixed_dfn/data"

# Load data:
norm_data    = read.csv(paste(data_path, "/combined_data.csv", sep = ""))
norm_data$X  = NULL

# Updating names with which are on the log-scale:
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
norm_data = norm_data[new_names]

# # Fix the early gypsum issues:
# early_gypsum_data = norm_data[early_gypsum_points]
# for(gyp_row in 1:nrow(early_gypsum_data)){
#   temp_row                         = early_gypsum_data[gyp_row,]
#   temp_row[which(is.na(temp_row))] = max(temp_row, na.rm = T)
#   temp_row[which(temp_row < 250)]  = max(temp_row, na.rm = T)
#   early_gypsum_data[gyp_row,]      = temp_row
# }
# norm_data[early_gypsum_points] = early_gypsum_data



# Breaking the input vars into different sets, separating out the output vars.
cols_of_log_vars              = c(2,3,5,6,8)
cols_of_net_vars              = c(9,10,11,12)
cols_of_non_fctn_vars         = c(2:14)
# Pe and Tau are perfectly correlated.  I'm removing tau
cols_of_fctn_vars             = c(15,16,17,18,19)
# I'm going to focus on calcite here for a bit:
output_vars                   = c(21:length(new_names))
# output_vars                   = c(21,22,24,26,28,30,32,34)
# output_vars                   = c(20,22,24,26,28,30,32,34) + 1
# For the fixed DFN, I no longer use any network variables
# input_vars_for_emulation_fctn = c(cols_of_net_vars[2], cols_of_fctn_vars)
input_vars_for_emulation_fctn = c(cols_of_fctn_vars)

# EDA on the functional variables shows that they should really be considered on the log-scale.
# This MIGHT lead to a change needed in the sobol calculation.  
norm_data[,input_vars_for_emulation_fctn] = log(norm_data[,input_vars_for_emulation_fctn])

# Make sure there is no order to these data:
norm_data = norm_data[sample(1:nrow(norm_data), nrow(norm_data)),]

# all_inputs            = c(input_vars_for_emulation_fctn,input_vars_for_emulation_nonfctn)
# zero_one              = function(x){(x-min(x, na.rm = T))/(max(x, na.rm = T)-min(x, na.rm = T))}
# norm_data[all_inputs] = as.data.frame(apply(norm_data[all_inputs], 2, zero_one))

# all_inputs            = c(input_vars_for_emulation_fctn,input_vars_for_emulation_nonfctn)
zero_one              = function(x){(x-min(x, na.rm = T))/(max(x, na.rm = T)-min(x, na.rm = T))}
norm_data             = as.data.frame(apply(norm_data, 2, zero_one))

testing_data     = norm_data[ceil(nrow(norm_data)*.75):nrow(norm_data),]
training_data    = norm_data[1:floor(nrow(norm_data)*.75),]
write.csv(training_data, "training_data.csv")
write.csv(testing_data,  "testing_data.csv")
write.csv(norm_data,     "full_data.csv")

# Double check that these are right:
new_names[cols_of_log_vars]
# should be:
# 'curr_inflow_rate', 'curr_diffusion_coef', 'curr_gypsum_rate_constant', 'curr_calcite_rate_constant',
# 'curr_calcite_surface_area'

new_names[cols_of_net_vars]
# should be:
# 'dfn_p32','dfn_volume', 'backbone_volume', 'backbone_p32', 

new_names[cols_of_fctn_vars]
# should be
# 'Pe', 'Da_1_gypsum', 'Da_2_gypsum', 'Da_1_calcite','Da_2_calcite', 'tau', 

##########
## Correlated Inputs (non-ftcn vars)
# The inputted variables are meant to be uncorrelated.  This means that I'll need to
# handle the functional inputs separately.  The input variables that we sample
# via a LHS are made to be uncorrelated.  The remaining question is are there
# (original) input variables that are correlated with the graph-measured vars.
# pairs(training_data[,cols_of_non_fctn_vars])
# Inflow rate seems to drive final gypsum and final calcite.  
# The graph-measured variables are all strongly positively correlated with
# one another.
# I am going to take only 1 graph-measured variable (dfn_p32).  I am going
# to remove final_calcite and final_gypsum.
# For the fixed DFN, I no longer use network variables here.
# input_vars_for_emulation_nonfctn = c(2:8, 10)
input_vars_for_emulation_nonfctn = c(2:8)

##########
## Correlated Inputs (w/ ftcn vars)
pairs(training_data[,input_vars_for_emulation_fctn])
# Hmmm...I think I need to get these relationships explicitly from Jeffrey.
# Also everything but dfn_p32 seems like it needs to be on the log-scale before 
# we emulate.
# training_data[,input_vars_for_emulation_fctn[2:length(input_vars_for_emulation_fctn)]] = log(training_data[,input_vars_for_emulation_fctn[2:length(input_vars_for_emulation_fctn)]])


#########
## Fitting the emulators.
# I'll fit an emulator on each output variable, for each var set (ftcn and non fctn).

## Non ftcn:
# Big note: if you change any of these calibrations, make sure that you also change the repeated
# tuning that happens in the boostrapping in MC_sobol_indices.R.
# tuning_2 = c('gypsum_flush_2_nondim', 'gypsum_flush_1_nondim','calcite_flush_3_nondim',
#              'calcite_flush_5_nondim', 'gypsum_flush_6_nondim', 'gypsum_flush_7_nondim',
#              'gypsum_flush_8_nondim', 'gypsum_flush_9_nondim', 'gypsum_flush_0_nondim', 
#              'gypsum_flush_10_nondim', 'gypsum_flush_11_nondim',
#              'gypsum_flush_12_nondim')
# tuning_3 = c('gypsum_flush_3_nondim',
#              'gypsum_flush_4_nondim', 'calcite_flush_2_nondim')

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
training_data_nonfctn = training_data[,input_vars_for_emulation_nonfctn]

f_nonfctn = function(var_num){
  # Organize the data in the proper way to input into mleHomGP.
  complete_cases    = which(!is.na(training_data[,var_num]))
  hetGP_inputs      = find_reps(X = as.matrix(training_data_nonfctn[complete_cases,]),
                                Z = training_data[complete_cases,var_num])
  # Create the model.
  model = NULL
  
  if(new_names[var_num] %in% tuning_2){
    model = mleHomGP( X = list(X0 = hetGP_inputs$X0, Z0 = hetGP_inputs$Z0, mult = hetGP_inputs$mult),
                      # noiseControl = list(g_min = 3),
                      Z = hetGP_inputs$Z,
                      lower = 0.5,
                      upper = 100 )
  } else if(new_names[var_num] %in% tuning_3){
    model = mleHomGP( X = list(X0 = hetGP_inputs$X0, Z0 = hetGP_inputs$Z0, mult = hetGP_inputs$mult),
                      # noiseControl = list(g_min = 3),
                      Z = hetGP_inputs$Z, covtype = "Matern3_2",
                      lower = 0.5,
                      upper = 100 )
  } else if(new_names[var_num] %in% tuning_4){
    model = mleHomGP( X = list(X0 = hetGP_inputs$X0, Z0 = hetGP_inputs$Z0, mult = hetGP_inputs$mult),
                      noiseControl = list(g_min = 3),
                      Z = hetGP_inputs$Z,
                      lower = 5,
                      upper = 110 )
  }else if(new_names[var_num] %in% tuning_5){
    model = mleHomGP( X = list(X0 = hetGP_inputs$X0, Z0 = hetGP_inputs$Z0, mult = hetGP_inputs$mult),
                      lower = 15,
                      Z = hetGP_inputs$Z,
                      upper = 500 )
  } else {
    
    model = mleHomGP( X = list(X0 = hetGP_inputs$X0, Z0 = hetGP_inputs$Z0, mult = hetGP_inputs$mult),
                      # noiseControl = list(g_min = 3),
                      lower = 0.5,
                      upper = 100,
                      Z = hetGP_inputs$Z, covtype = "Matern3_2" )
  }
  
  name_of_model = paste("non_ftcn_models/het_gp_output_", 
                        new_names[var_num], ".Rdata", sep = "")
  save(model, file=name_of_model)
}

mclapply(output_vars, f_nonfctn, mc.cores = 10)
# output_vars = c(20,22,24,26,28,30,32,34) + 2
# mclapply(output_vars, f_nonfctn, mc.cores = 10)

## ftcn:
tuning_2 = c('gypsum_flush_9_nondim',
             'gypsum_flush_10_nondim', 'gypsum_flush_11_nondim', 'gypsum_flush_12_nondim',
             'gypsum_flush_0_nondim', "gypsum_flush_1_nondim", 'calcite_flush_8_nondim',
             'calcite_flush_9_nondim','calcite_flush_10_nondim',
             'calcite_flush_11_nondim','calcite_flush_12_nondim',
             "gypsum_flush_8_nondim",
             "calcite_flush_2_nondim")
tuning_3 = c('calcite_flush_1_nondim',
             'gypsum_flush_12_nondim',
             'gypsum_flush_11_nondim')
tuning_4 = c("calcite_flush_3_nondim", "calcite_flush_4_nondim","calcite_flush_5_nondim", 'gypsum_flush_7_nondim')
tuning_5 = c('gypsum_flush_9_nondim')

training_data_fctn = training_data[,input_vars_for_emulation_fctn]

f_fctn = function(var_num){
  complete_cases    = which(!is.na(training_data[,var_num]))
  hetGP_inputs      = find_reps(X = as.matrix(training_data_fctn[complete_cases,]),
                                Z = training_data[complete_cases,var_num])
  if(new_names[var_num] %in% tuning_2){
    model = mleHomGP( X = list(X0 = hetGP_inputs$X0, Z0 = hetGP_inputs$Z0, mult = hetGP_inputs$mult),
                      noiseControl = list(g_min = 3),
                      Z = hetGP_inputs$Z,
                      lower = 0.5,
                      upper = 100 )
  }else if(new_names[var_num] %in% tuning_3){
    model = mleHomGP( X = list(X0 = hetGP_inputs$X0, Z0 = hetGP_inputs$Z0, mult = hetGP_inputs$mult),
                      noiseControl = list(g_min = 3),
                      Z = hetGP_inputs$Z, covtype = "Matern3_2",
                      lower = 0.5,
                      upper = 100 )
  } else if(new_names[var_num] %in% tuning_4){
    model = mleHomGP( X = list(X0 = hetGP_inputs$X0, Z0 = hetGP_inputs$Z0, mult = hetGP_inputs$mult),
                      noiseControl = list(g_min = 3),
                      Z = hetGP_inputs$Z,
                      lower = 15,
                      upper = 100 )
  } else if(new_names[var_num] %in% tuning_5){
    model = mleHomGP( X = list(X0 = hetGP_inputs$X0, Z0 = hetGP_inputs$Z0, mult = hetGP_inputs$mult),
                      noiseControl = list(g_min = 3),
                      Z = hetGP_inputs$Z,
                      lower = 100,
                      upper = 500, covtype = "Matern3_2" )
  } else {

    model = mleHomGP( X = list(X0 = hetGP_inputs$X0, Z0 = hetGP_inputs$Z0, mult = hetGP_inputs$mult),
                      noiseControl = list(g_min = 3),
                      Z = hetGP_inputs$Z, covtype = "Matern3_2" )
  }
  name_of_model = paste("ftcn_models/het_gp_output_",
                        new_names[var_num], ".Rdata", sep = "")
  save(model, file=name_of_model)
}

# output_vars = c(20,22,24,26,28,30,32,34) + 1
mclapply(output_vars, f_fctn, mc.cores = 10)
# output_vars = c(20,22,24,26,28,30,32,34) + 2
# mclapply(output_vars, f_fctn, mc.cores = 10)
