###############################################################################
### After fitting the final GP, this script graphs the pivoted cholesky 
### transformation on the prediction errors as a means to assess the model fit.
### I'll likely do this on both the training and testing data.
## Author: Alexander C. Murph
## Date: September 2023
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

#################################
##### Non functional variables calibration assessment
pdf("choleskey_analysis_with_nonfctn_vars.pdf")
for(output_var_idx in output_vars){
  model_name = paste("../joint_emulation/non_ftcn_models/het_gp_output_",new_names[output_var_idx], ".Rdata", sep = "")
  model = NULL
  load(model_name)
  
  ####################################
  ## Get the predictions using the model:
  complete_train      = which(!is.na(training_data[,output_var_idx]))
  het_GP_inputs_train = find_reps(X = as.matrix(training_data[complete_train,input_vars_for_emulation_nonfctn]),
                                  Z = training_data[complete_train,output_var_idx])
  complete_test       = which(!is.na(testing_data[,output_var_idx]))
  het_GP_inputs_test  = find_reps(X = as.matrix(testing_data[complete_test,input_vars_for_emulation_nonfctn]),
                                  Z = testing_data[complete_test,output_var_idx])
  
  if(anyNA(testing_data[complete_test,output_var_idx])) browser()
  
  Z_train = het_GP_inputs_train$Z
  Z_test  = het_GP_inputs_test$Z
  
  n_train     = length(Z_train)
  train_preds = predict(object=model, 
                        x = het_GP_inputs_train$X0, 
                        xprime = het_GP_inputs_train$X0)
  train_true  = as.matrix(Z_train)
  
  n_test      = length(Z_test)
  test_preds  = predict(object=model, 
                        x = het_GP_inputs_test$X0, 
                        xprime = het_GP_inputs_test$X0)
  test_true   = Z_test
  
  ####################################
  ## Get the pivoted cholesky decomposition matrices and rescale prediction errors
  cov_train = train_preds$cov + diag(train_preds$nugs)
  
  # In the following, the decomposition on the cov matrix is done in the order
  # or greatest to least greatest variance.
  chol_train  = chol(cov_train, pivot = TRUE)
  pivot_train = attr(chol_train, "pivot")
  
  # We get the prediction errors:
  pivoted_error_terms_train = (train_true - as.matrix(train_preds$mean))
  
  # Then order them according to the pivott in our cholesky decomposition:
  pivoted_error_terms_train = pivoted_error_terms_train[pivot_train]
  pivoted_cholesky_transformed_errors_train = solve(chol_train, pivoted_error_terms_train) #pivoted_error_terms_train #solve(chol_train)%*%
  
  # Finally, we do the same thing for the testing data:
  cov_test                 = test_preds$cov + diag(test_preds$nugs)
  chol_test                = chol(cov_test, pivot = TRUE)
  pivot_test               = attr(chol_test, "pivot")
  pivoted_error_terms_test = (test_true - as.matrix(test_preds$mean))
  
  pivoted_error_terms_test                 = pivoted_error_terms_test[pivot_test]
  pivoted_cholesky_transformed_errors_test = solve(chol_test, pivoted_error_terms_test) #pivoted_error_terms_test #solve(chol_test)%*%
  
  ####################################
  ## Get the student-t bounds:
  p = 5
  df_train = n_train - p
  upper_t_train = qt(0.975, df_train)
  lower_t_train = qt(0.025, df_train)
  
  df_test = n_test - p
  upper_t_test = qt(0.975, df_test)
  lower_t_test = qt(0.025, df_test)
  
  ####################################
  ## Create data and produce graphs.
  lst_p = list()
  lst_p[[1]] = NULL
  
  # I want a graph for every possible input variable:
  count          = 1
  full_X_testing = NULL
  full_X_testing = testing_data[complete_test,input_vars_for_emulation_nonfctn]
  for(var_num in input_vars_for_emulation_nonfctn){
    
    data_name = paste("pivot_data_for_train", var_num, sep = "")
    eval(call("<-", as.name(data_name), data.frame(X = full_X_testing[,count], Y = as.numeric(as.vector(pivoted_cholesky_transformed_errors_test) )) ))
    
    x = paste("g",var_num, sep = "")
    eval(parse(text = paste(as.name(x), "<-", "ggplot(", as.name(data_name), ', aes(x = X, y = Y)) + geom_point() + 
                geom_hline(yintercept= c(lower_t_test, upper_t_test, 0), linetype=c("dashed", "dashed", "solid"))+  
                ggtitle(TeX("")) + 
                theme_bw() + 
                xlab(TeX(paste("Scaled", new_names[var_num]))) + 
                ylab(TeX("Pivoted Chol Error"))', sep = "") #ylim(c(-3,3))
    ))
    # print(paste("ggMarginal(", as.name(x), ', margins = "y", type = "density", size = 7, fill = "gray", alpha = 0.2)', sep = ""))
    lst_p[[count]] = eval(parse(text = paste("ggMarginal(", as.name(x), ', margins = "y", type = "density", size = 7, fill = "gray", alpha = 0.2)', sep = "") ))
    count          = count + 1
  }
  
  do.call("grid.arrange", c(lst_p, ncol=2,top = paste("Cholesky Errors for Output Variable", new_names[output_var_idx], "on Testing Data")))
}
dev.off()

#################################
##### Functional variables calibration assessment
pdf("choleskey_analysis_with_fctn_vars.pdf")
for(output_var_idx in output_vars){
  model_name = paste("../joint_emulation/ftcn_models/het_gp_output_",new_names[output_var_idx], ".Rdata", sep = "")
  model      = NULL
  load(model_name)
  
  ####################################
  ## Get the predictions using the model:
  complete_train      = which(!is.na(training_data[complete_train,output_var_idx]))
  het_GP_inputs_train = find_reps(X = as.matrix(training_data[complete_train,input_vars_for_emulation_fctn]),
                                  Z = training_data[complete_train,output_var_idx])
  complete_test       = which(!is.na(testing_data[,output_var_idx]))
  het_GP_inputs_test  = find_reps(X = as.matrix(testing_data[complete_test,input_vars_for_emulation_fctn]),
                                  Z = testing_data[complete_test,output_var_idx])
  Z_train             = het_GP_inputs_train$Z
  Z_test              = het_GP_inputs_test$Z
  
  n_train     = length(Z_train)
  train_preds = predict(object=model, 
                        x = het_GP_inputs_train$X0, 
                        xprime = het_GP_inputs_train$X0)
  train_true  = as.matrix(Z_train)
  
  n_test      = length(Z_test)
  test_preds  = predict(object=model, 
                        x = het_GP_inputs_test$X0, 
                        xprime = het_GP_inputs_test$X0)
  test_true   = Z_test
  
  ####################################
  ## Get the pivoted cholesky decomposition matrices and rescale prediction errors
  cov_train = train_preds$cov + diag(train_preds$nugs)
  
  # In the following, the decomposition on the cov matrix is done in the order
  # or greatest to least greatest variance.
  chol_train  = chol(cov_train, pivot = TRUE)
  pivot_train = attr(chol_train, "pivot")
  
  # We get the prediction errors:
  pivoted_error_terms_train = (train_true - as.matrix(train_preds$mean))
  
  # Then order them according to the pivott in our cholesky decomposition:
  pivoted_error_terms_train                 = pivoted_error_terms_train[pivot_train]
  pivoted_cholesky_transformed_errors_train = solve(chol_train)%*%pivoted_error_terms_train
  
  # Finally, we do the same thing for the testing data:
  cov_test                 = test_preds$cov + diag(test_preds$nugs)
  chol_test                = chol(cov_test, pivot = TRUE)
  pivot_test               = attr(chol_test, "pivot")
  pivoted_error_terms_test = (test_true - as.matrix(test_preds$mean))
  
  pivoted_error_terms_test = pivoted_error_terms_test[pivot_test]
  
  # cov_test = cov_test[pivot_test,pivot_test]
  # pivoted_cholesky_transformed_errors_test = solve(chol_test)%*%pivoted_error_terms_test
  # pivoted_cholesky_transformed_errors_test = solve(chol_test, pivoted_error_terms_test)
  pivoted_cholesky_transformed_errors_test = backsolve(chol_test, diag(dim(chol_test)[1]))%*%pivoted_error_terms_test
  
  # if(new_names[output_var_idx] == "calcite_flush_2_nondim") browser()
  
  ####################################
  ## Get the student-t bounds:
  p             = 5
  df_train      = n_train - p
  upper_t_train = qt(0.975, df_train)
  lower_t_train = qt(0.025, df_train)
  
  df_test      = n_test - p
  upper_t_test = qt(0.975, df_test)
  lower_t_test = qt(0.025, df_test)
  
  lst_p      = list()
  lst_p[[1]] = NULL
  
  # I want a graph for every possible input variable:
  count          = 1
  full_X_testing = NULL
  full_X_testing = testing_data[complete_test,input_vars_for_emulation_fctn]
  for(var_num in input_vars_for_emulation_fctn){
    
    data_name = paste("pivot_data_for_train", var_num, sep = "")
    eval(call("<-", as.name(data_name), data.frame(X = full_X_testing[,count], Y = as.numeric(as.vector(pivoted_cholesky_transformed_errors_test) )) ))
    
    x = paste("g",var_num, sep = "")
    eval(parse(text = paste(as.name(x), "<-", "ggplot(", as.name(data_name), ', aes(x = X, y = Y)) + geom_point() + 
                geom_hline(yintercept= c(lower_t_test, upper_t_test, 0), linetype=c("dashed", "dashed", "solid"))+  
                ggtitle(TeX("")) + 
                theme_bw() + 
                xlab(TeX(paste("Scaled", new_names[var_num]))) + 
                ylab(TeX("Pivoted Chol Error"))', sep = "") #ylim(c(-3,3))
    ))
    # print(paste("ggMarginal(", as.name(x), ', margins = "y", type = "density", size = 7, fill = "gray", alpha = 0.2)', sep = ""))
    lst_p[[count]] = eval(parse(text = paste("ggMarginal(", as.name(x), ', margins = "y", type = "density", size = 7, fill = "gray", alpha = 0.2)', sep = "") ))
    count          = count + 1
  }
  
  do.call("grid.arrange", c(lst_p, ncol=2,top = paste("Cholesky Errors for Output Variable", new_names[output_var_idx], "on Testing Data")))
}
dev.off()

