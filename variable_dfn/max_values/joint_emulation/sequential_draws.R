# Sequential design on input locations for dfnWorks study.
# This is an updated code suite where I update some issues I had created
# by using bkde unnecessarily.  For this file, I just use the quantile
# function with 'percentile.'
## Author: Alexander C. Murph
## Date: August 2023
library(hetGP)
setwd("~/GitLab/sa_for_chemistry/variable_dfn/max_values/joint_emulation")
source("sequential_helpers.R")

# This is not a simulation study, so I'll set the seed.
set.seed(13)

# Data location:
data_path = "~/GitLab/sa_for_chemistry/variable_dfn/max_values/data"

# Load data:
norm_data    = read.csv(paste(data_path, "/combined_data.csv", sep = ""))
norm_data$X  = NULL

# Updating names with which are on the log-scale:
new_names = c("dfn_seed",                   "curr_inflow_rate",           "curr_diffusion_coef",        "curr_porosity",              "curr_gypsum_rate_constant", 
              "curr_calcite_rate_constant", "curr_gypsum_surface_area",   "curr_calcite_surface_area",  "dfn_p32",                    "dfn_volume",                
              "backbone_volume",            "backbone_p32",               "final_gypsum",               "final_calcite",              "Pe",                        
              "Da_1_gypsum",                "Da_2_gypsum",                "Da_1_calcite",               "Da_2_calcite",               "tau",                       
              "max_calcite",                "min_gypsum")

## Based on the EDA, I am going to drop the initial calcite points and do some fixes to the gypsum data.
# Drop the first few calcite:
# new_names = new_names[which(!(new_names%in%taopoints_to_remove))]
norm_data             = norm_data[new_names]
norm_data$max_calcite = log(norm_data$max_calcite)

# Breaking the input vars into different sets, separating out the output vars.
cols_of_log_vars              = c(2,3,5,6,8)
cols_of_net_vars              = c(9,10,11,12)
cols_of_non_fctn_vars         = c(2:14)
cols_of_fctn_vars             = c(15,16,17,18,19,20)
output_vars                   = c(21:length(new_names))
input_vars_for_emulation_fctn = c(cols_of_net_vars[2], cols_of_fctn_vars)

# EDA on the functional variables shows that they should really be considered on the log-scale.
# This MIGHT lead to a change needed in the sobol calculation.  
norm_data[,input_vars_for_emulation_fctn[2:length(input_vars_for_emulation_fctn)]] = log(norm_data[,input_vars_for_emulation_fctn[2:length(input_vars_for_emulation_fctn)]])

# Make sure there is no order to these data:
norm_data = norm_data[sample(1:nrow(norm_data), nrow(norm_data)),]

# From here, we rescale.  Then split and save the data.
zero_one         = function(x){(x-min(x, na.rm = T))/(max(x, na.rm = T)-min(x, na.rm = T))}
norm_data        = apply(norm_data, 2, zero_one)
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
pairs(training_data[,cols_of_non_fctn_vars])
# Inflow rate seems to drive final gypsum and final calcite.  
# The graph-measured variables are all strongly positively correlated with
# one another.
# I am going to take only 1 graph-measured variable (dfn_p32).  I am going
# to remove final_calcite and final_gypsum.
input_vars_for_emulation_nonfctn = c(2:8, 10)

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
tuning_2 = c('gypsum_flush_2_nondim', 'gypsum_flush_1_nondim',
             'calcite_flush_3_nondim', 'calcite_flush_5_nondim', 
             'gypsum_flush_6_nondim', 'gypsum_flush_7_nondim',
             'gypsum_flush_8_nondim', 'gypsum_flush_9_nondim', 
             'gypsum_flush_0_nondim', 'gypsum_flush_10_nondim', 
             'gypsum_flush_11_nondim', 'gypsum_flush_12_nondim')
tuning_3 = c('gypsum_flush_3_nondim','gypsum_flush_4_nondim', 
             'calcite_flush_2_nondim')
tuning_4 = c()
training_data_nonfctn = training_data[,input_vars_for_emulation_nonfctn]
for(var_num in output_vars){
  # Organize the data in the proper way to input into mleHetGP.
  complete_cases    = which(!is.na(training_data[,var_num]))
  hetGP_inputs      = find_reps(X = as.matrix(training_data_nonfctn[complete_cases,]),
                                Z = training_data[complete_cases,var_num])
  # Create the model.
  model = NULL
  
  if(new_names[var_num] %in% tuning_2){
    model = mleHetGP( X = list(X0 = hetGP_inputs$X0, Z0 = hetGP_inputs$Z0, mult = hetGP_inputs$mult),
                      noiseControl = list(g_min = 3),
                      Z = hetGP_inputs$Z,
                      settings = list(checkHom = FALSE ) )
  } else if(new_names[var_num] %in% tuning_3){
    model = mleHetGP( X = list(X0 = hetGP_inputs$X0, Z0 = hetGP_inputs$Z0, mult = hetGP_inputs$mult),
                      noiseControl = list(g_min = 3),
                      Z = hetGP_inputs$Z, covtype = "Matern3_2",
                      settings = list(checkHom = FALSE ) )
  } else if(new_names[var_num] %in% tuning_4){
    model = mleHetGP( X = list(X0 = hetGP_inputs$X0, Z0 = hetGP_inputs$Z0, mult = hetGP_inputs$mult),
                      Z = hetGP_inputs$Z,
                      noiseControl = list(g_min = 3),
                      settings = list(checkHom = FALSE ) )
  } else {
    
    model = mleHetGP( X = list(X0 = hetGP_inputs$X0, Z0 = hetGP_inputs$Z0, mult = hetGP_inputs$mult),
                      noiseControl = list(g_min = 3),
                      lower = 1.5,
                      upper = 110,
                      Z = hetGP_inputs$Z,
                      settings = list(checkHom = FALSE ) )
  }
  
  name_of_model = paste("non_ftcn_models/het_gp_output_", 
                        new_names[var_num], ".Rdata", sep = "")
  save(model, file=name_of_model)
}

## ftcn:
tuning_2 = c('calcite_flush_3_nondim', 'gypsum_flush_6_nondim',
             'gypsum_flush_7_nondim', 'gypsum_flush_9_nondim',
             'gypsum_flush_10_nondim', 'gypsum_flush_11_nondim', 
             'gypsum_flush_12_nondim', 'gypsum_flush_0_nondim', 
             "gypsum_flush_1_nondim", 'calcite_flush_8_nondim',
             'calcite_flush_9_nondim', 'calcite_flush_10_nondim',
             'calcite_flush_11_nondim', 'calcite_flush_12_nondim', 
             "gypsum_flush_8_nondim", "calcite_flush_2_nondim")
tuning_3 = c('calcite_flush_1_nondim', 'gypsum_flush_12_nondim', 
             'gypsum_flush_11_nondim', 'gypsum_flush_9_nondim')
tuning_4 = c('gypsum_flush_10_nondim')
training_data_fctn = training_data[,input_vars_for_emulation_fctn]


for(var_num in output_vars){
  complete_cases    = which(!is.na(training_data[,var_num]))
  hetGP_inputs      = find_reps(X = as.matrix(training_data_fctn[complete_cases,]),
                                Z = training_data[complete_cases,var_num])
  if(new_names[var_num] %in% tuning_2){
    model = mleHetGP( X = list(X0 = hetGP_inputs$X0, Z0 = hetGP_inputs$Z0, mult = hetGP_inputs$mult),
                      noiseControl = list(g_min = 3),
                      Z = hetGP_inputs$Z,
                      settings = list(checkHom = FALSE ) )
  }else if(new_names[var_num] %in% tuning_3){
    model = mleHetGP( X = list(X0 = hetGP_inputs$X0, Z0 = hetGP_inputs$Z0, mult = hetGP_inputs$mult),
                      noiseControl = list(g_min = 3),
                      Z = hetGP_inputs$Z, covtype = "Matern3_2",
                      settings = list(checkHom = FALSE ) )
  } else if(new_names[var_num] %in% tuning_4){
    model = mleHetGP( X = list(X0 = hetGP_inputs$X0, Z0 = hetGP_inputs$Z0, mult = hetGP_inputs$mult),
                      noiseControl = list(g_min = 3),
                      Z = hetGP_inputs$Z, covtype = "Matern5_2",
                      settings = list(checkHom = FALSE ) )
  } else {
    
    model = mleHetGP( X = list(X0 = hetGP_inputs$X0, Z0 = hetGP_inputs$Z0, mult = hetGP_inputs$mult),
                      noiseControl = list(g_min = 3),
                      Z = hetGP_inputs$Z,
                      lower = 1.5,
                      upper = 110,
                      settings = list(checkHom = FALSE ) )
  }
  name_of_model = paste("ftcn_models/het_gp_output_", 
                        new_names[var_num], ".Rdata", sep = "")
  save(model, file=name_of_model)
}



