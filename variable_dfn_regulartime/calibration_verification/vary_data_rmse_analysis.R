###############################################################################
### Trying to assess if I need more data.
## Author: Alexander C. Murph
## Date: January 2024
library(Matrix)
library(hetGP)
library(ggplot2)
library(latex2exp)
library(gridExtra)
library(ggExtra) 
library(rlang)
setwd("~/GitLab/sa_for_chemistry/variable_dfn_regulartime/calibration_verification")

# Load in the full training data:
training_data   = read.csv(file="../joint_emulation/training_data.csv") 
training_data$X = NULL
# Load in the testing data:
testing_data    = read.csv(file="../joint_emulation/testing_data.csv") 
testing_data$X  = NULL

# Variable info:
new_names = c("dfn_seed","log_curr_inflow_rate","log_curr_diffusion_coef","curr_porosity",
              "log_curr_gypsum_rate_constant","log_curr_calcite_rate_constant","curr_gypsum_surface_area",
              "log_curr_calcite_surface_area","dfn_p32","dfn_volume","backbone_volume",
              "backbone_p32","final_gypsum","final_calcite","Pe","Da_1_gypsum","Da_2_gypsum",
              "Da_1_calcite","Da_2_calcite","tau","calcite_flush_1_nondim",
              "gypsum_flush_1_nondim","calcite_flush_10_nondim","gypsum_flush_10_nondim",
              "calcite_flush_50_nondim","gypsum_flush_50_nondim","calcite_flush_100_nondim",
              "gypsum_flush_100_nondim","calcite_flush_200_nondim","gypsum_flush_200_nondim",
              "calcite_flush_500_nondim","gypsum_flush_500_nondim","calcite_flush_1000_nondim",
              "gypsum_flush_1000_nondim","calcite_flush_5000_nondim","gypsum_flush_5000_nondim",
              "calcite_flush_10000_nondim","gypsum_flush_10000_nondim")
cols_of_log_vars                 = c(2,3,5,6,8)
cols_of_net_vars                 = c(9,10,11,12)
cols_of_non_fctn_vars            = c(2:14)
cols_of_fctn_vars                = c(15,16,17,18,19,20)
output_vars                      = c(21:38)
input_vars_for_emulation_nonfctn = c(2:9)
input_vars_for_emulation_fctn    = c(cols_of_net_vars[1], cols_of_fctn_vars)
rmses                            = NULL

# Perform this experiment for increasing data sizes and record rmse.
data_percentages = seq(from = 0.5, to = 1, by = 0.1)
orig_norm_data   = read.csv(file = "../joint_emulation/full_data.csv")
n_orig           = nrow(orig_norm_data)
# Make sure there is no implicit order to the observations:
orig_norm_data   = orig_norm_data[sample(1:n_orig, size = n_orig),]

for(data_percent in data_percentages){
  # Grab the data, subset it properly, and remake the models.
  norm_data        = orig_norm_data[1:(floor(data_percent*n_orig)),]
  testing_data     = norm_data[ceil(nrow(norm_data)*.75):nrow(norm_data),]
  training_data    = norm_data[1:floor(nrow(norm_data)*.75),]
  
  #################################
  ##### Non functional variables calibration assessment
  tuning_2 = c('calcite_flush_10_nondim', 'gypsum_flush_50_nondim', 'gypsum_flush_100_nondim', 
               'gypsum_flush_200_nondim', 'gypsum_flush_500_nondim', 'gypsum_flush_1000_nondim',
               'gypsum_flush_5000_nondim', 'gypsum_flush_10000_nondim')
  
  output_qois           = c()
  p                     = length(input_vars_for_emulation_nonfctn)
  training_data_nonfctn = training_data[,input_vars_for_emulation_nonfctn]
  testing_data_nonfctn  = testing_data[,input_vars_for_emulation_nonfctn]
  for(output_var_idx in output_vars){
    # Organize the data in the proper way to input into mleHetGP.
    hetGP_inputs      = find_reps(X = as.matrix(training_data_nonfctn),
                                  Z = training_data[,output_var_idx])
    # Create the model.
    model = NULL
    
    if(new_names[output_var_idx] %in% tuning_2){
      model = mleHetGP( X = list(X0 = hetGP_inputs$X0, Z0 = hetGP_inputs$Z0, mult = hetGP_inputs$mult),
                        Z = hetGP_inputs$Z,
                        settings = list(checkHom = FALSE ) )
    } else {
      
      model = mleHetGP( X = list(X0 = hetGP_inputs$X0, Z0 = hetGP_inputs$Z0, mult = hetGP_inputs$mult),
                        noiseControl = list(g_min = 3),
                        Z = hetGP_inputs$Z, covtype = "Matern3_2",
                        settings = list(checkHom = FALSE ) )
    }
    
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
    rmse             = sqrt( mean( (error_terms_test)**2 ) )
    
    output_var_parts = unlist(strsplit(new_names[output_var_idx], "[_]"))
    output_var_type  = output_var_parts[1]
    output_var_new   = paste(output_var_parts[3], output_var_parts[1], sep = '_')
    flush_level      = paste(output_var_parts[3], output_var_parts[2], output_var_parts[4], sep = '_')
    
    temp_row         = data.frame(rmse = rmse, output_var = flush_level, 
                                  input_type = "Non-Functional", output_var_type = output_var_type, 
                                  group_type = paste("Non-Functional", output_var_type),
                                  data_percentage = data_percent)
    rmses            = rbind(rmses, temp_row)
    output_qois      = c(output_qois, flush_level)
  }
  
  
  #################################
  ##### Functional variables calibration assessment
  tuning_2 = c('gypsum_flush_1_nondim', 
               'calcite_flush_10_nondim', 'gypsum_flush_50_nondim',
               'gypsum_flush_100_nondim', 
               'gypsum_flush_200_nondim', 
               'gypsum_flush_500_nondim',
               'gypsum_flush_5000_nondim',
               'gypsum_flush_10000_nondim')
  tuning_3 = c('calcite_flush_10000_nondim', 'calcite_flush_5000_nondim', 'calcite_flush_1000_nondim',
               'gypsum_flush_1000_nondim', 
               'calcite_flush_500_nondim', 'calcite_flush_200_nondim', 'calcite_flush_100_nondim','calcite_flush_1_nondim')
  
  p                     = length(input_vars_for_emulation_fctn)
  training_data_fctn    = training_data[,input_vars_for_emulation_fctn]
  testing_data_fctn     = testing_data[,input_vars_for_emulation_fctn]
  for(output_var_idx in output_vars){
    
    hetGP_inputs      = find_reps(X = as.matrix(training_data_fctn),
                                  Z = training_data[,output_var_idx])
    if(new_names[output_var_idx] %in% tuning_2){
      model = mleHetGP( X = list(X0 = hetGP_inputs$X0, Z0 = hetGP_inputs$Z0, mult = hetGP_inputs$mult),
                        Z = hetGP_inputs$Z,
                        settings = list(checkHom = FALSE ) )
    }else if(new_names[output_var_idx] %in% tuning_3){
      model = mleHetGP( X = list(X0 = hetGP_inputs$X0, Z0 = hetGP_inputs$Z0, mult = hetGP_inputs$mult),
                        Z = hetGP_inputs$Z, covtype = "Matern3_2",
                        settings = list(checkHom = FALSE ) )
    } else {
      
      model = mleHetGP( X = list(X0 = hetGP_inputs$X0, Z0 = hetGP_inputs$Z0, mult = hetGP_inputs$mult),
                        noiseControl = list(g_min = 3),
                        Z = hetGP_inputs$Z, covtype = "Matern3_2",
                        settings = list(checkHom = FALSE ) )
    }
    
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
    rmse             = sqrt( mean( (error_terms_test)**2 ) )
    
    output_var_parts = unlist(strsplit(new_names[output_var_idx], "[_]"))
    output_var_type  = output_var_parts[1]
    flush_level      = paste(output_var_parts[3], output_var_parts[2], output_var_parts[4], sep = '_')
    
    temp_row         = data.frame(rmse = rmse, output_var = flush_level, 
                                  input_type = "Functional", output_var_type = output_var_type, 
                                  group_type = paste("Functional", output_var_type),
                                  data_percentage = data_percent)
    rmses            = rbind(rmses, temp_row)
  }
}



rmses$output_var = factor(rmses$output_var, levels = unique(output_qois))



count = 1
lst_p = NULL
for(qoi in unique(output_qois)){
  subset_data = rmses[which(rmses$output_var == qoi),]
  lst_p[[count]] = ggplot(subset_data, aes(x = data_percentage, y = rmse, group = interaction(input_type, output_var_type))) + geom_line(aes(linetype=input_type, color = output_var_type)) +
                          xlab("Output QoI") + theme(axis.text.x = element_text(angle = 270, vjust = 0.5, hjust=1)) + ggtitle(paste("RMSE for Predicting", qoi, "at Varying Data Sizes")) + ylim(c(0, 0.4))
  count = count + 1
}

glist = lapply(lst_p, ggplotGrob)
pdf("rmses_varying_data.pdf")
# do.call("grid.arrange", c(lst_p, nrow=2))
marrangeGrob(grobs=glist, nrow=2, ncol=2)
dev.off()





