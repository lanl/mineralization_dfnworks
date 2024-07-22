###############################################################################
### Examine the RMSE for the models created in the joint_emulation folder.
## Author: Alexander C. Murph
## Date: January 2024
library(Matrix)
library(hetGP)
library(ggplot2)
library(latex2exp)
library(gridExtra)
library(ggExtra) 
library(rlang)
setwd("~/GitLab/sa_for_chemistry/variable_dfn/max_values/calibration_verification")

# Load in the full training data:
training_data   = read.csv(file="../joint_emulation/training_data.csv") 
training_data$X = NULL
# Load in the testing data:
testing_data    = read.csv(file="../joint_emulation/testing_data.csv") 
testing_data$X  = NULL

# Variable info:
new_names = c("dfn_seed",                   "curr_inflow_rate",           "curr_diffusion_coef",        "curr_porosity",              "curr_gypsum_rate_constant", 
              "curr_calcite_rate_constant", "curr_gypsum_surface_area",   "curr_calcite_surface_area",  "dfn_p32",                    "dfn_volume",                
              "backbone_volume",            "backbone_p32",               "final_gypsum",               "final_calcite",              "Pe",                        
              "Da_1_gypsum",                "Da_2_gypsum",                "Da_1_calcite",               "Da_2_calcite",               "tau",                       
              "max_calcite",                "min_gypsum")

cols_of_log_vars                 = c(2,3,5,6,8)
cols_of_net_vars                 = c(9,10,11,12)
cols_of_non_fctn_vars            = c(2:14)
cols_of_fctn_vars                = c(15,16,17,18,19,20)
output_vars                      = c(21:length(new_names))
input_vars_for_emulation_nonfctn = c(2:8,10)
input_vars_for_emulation_fctn    = c(cols_of_net_vars[2], cols_of_fctn_vars)
rmses                            = NULL


#################################
##### Non functional variables calibration assessment
output_qois = c()
p           = length(input_vars_for_emulation_nonfctn)
for(output_var_idx in output_vars){
  model_name = paste("../joint_emulation/non_ftcn_models/het_gp_output_",new_names[output_var_idx], ".Rdata", sep = "")
  model      = NULL
  load(model_name)
  
  ####################################
  ## Get the predictions using the model:
  het_GP_inputs_test  = find_reps(X = as.matrix(testing_data[,input_vars_for_emulation_nonfctn]),
                                  Z = testing_data[,output_var_idx])
  Z_test              = het_GP_inputs_test$Z
  
  n_test      = length(Z_test)
  test_preds  = predict(object=model, 
                        x = het_GP_inputs_test$X0, 
                        xprime = het_GP_inputs_test$X0)
  test_true   = Z_test
  
  # Finally, we do the same thing for the testing data:
  error_terms_test = (test_true - as.matrix(test_preds$mean))
  rmse             = sqrt( mean( (error_terms_test)**2, na.rm = T ) )
  
  output_var_parts = unlist(strsplit(new_names[output_var_idx], "[_]"))
  output_var_type  = output_var_parts[1]
  output_var_new   = paste(output_var_parts[3], output_var_parts[1], sep = '_')
  flush_level      = paste(output_var_parts[3], output_var_parts[2], output_var_parts[4], sep = '_')
  
  temp_row         = data.frame(rmse = rmse, output_var = flush_level, input_type = "Non-Functional", output_var_type = output_var_type, group_type = paste("Non-Functional", output_var_type))
  rmses            = rbind(rmses, temp_row)
  output_qois      = c(output_qois, flush_level)
}


#################################
##### Functional variables calibration assessment
for(output_var_idx in output_vars){
  # output_var_idx = 1
  model_name = paste("../joint_emulation/ftcn_models/het_gp_output_",new_names[output_var_idx], ".Rdata", sep = "")
  model      = NULL
  load(model_name)
  
  ####################################
  ## Get the predictions using the model:
  het_GP_inputs_test  = find_reps(X = as.matrix(testing_data[,input_vars_for_emulation_fctn]),
                                  Z = testing_data[,output_var_idx])
  Z_test              = het_GP_inputs_test$Z
  
  n_test      = length(Z_test)
  test_preds  = predict(object=model, 
                        x = het_GP_inputs_test$X0, 
                        xprime = het_GP_inputs_test$X0)
  test_true   = Z_test
  
  # Finally, we do the same thing for the testing data:
  error_terms_test = (test_true - as.matrix(test_preds$mean))
  rmse             = sqrt( mean( (error_terms_test)**2, na.rm = T ) )
  
  output_var_parts = unlist(strsplit(new_names[output_var_idx], "[_]"))
  output_var_type  = output_var_parts[1]
  flush_level      = paste(output_var_parts[3], output_var_parts[2], output_var_parts[4], sep = '_')
  
  temp_row         = data.frame(rmse = rmse, output_var = flush_level, input_type = "Functional", output_var_type = output_var_type, group_type = paste("Functional", output_var_type))
  rmses            = rbind(rmses, temp_row)
}

rmses$output_var = factor(rmses$output_var, levels = unique(output_qois))

rmses

# pdf("rmses_final_model.pdf")
# ggplot(rmses, aes(x = as.factor(output_var), y = rmse, group = interaction(input_type, output_var_type))) + geom_line(aes(linetype=input_type, color = output_var_type)) +
#   xlab("Output QoI") + theme(axis.text.x = element_text(angle = 270, vjust = 0.5, hjust=1)) + ggtitle("RMSEs for Final Models")
# dev.off()


