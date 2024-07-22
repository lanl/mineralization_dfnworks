###############################################################################
### Script to analysis the GPs fit during the sequential design, but remaking
### the GPs WITHOUT the upper outliers.
## Author: Alexander C. Murph
## Date: August 2023
library(hetGP)
library(ggplot2)
library(latex2exp)
setwd("/Users/murph/dfnworks_variability/gp_analysis")
source("GP_fit_analysis_helpers.R")

# Make sure that they following is the same as in the sequential_draws file:
num_of_replications = 10
num_of_experiments  = 50
num_of_additional_input_locations = 50
total_num_from_initial_experiment = num_of_replications*num_of_experiments
num_of_additional_replicates      = 5
percentile = 0.1

## Data locations
test_data_loc = "/Users/murph/dfnworks_variability/test_data"
models_loc = "/Users/murph/dfnworks_variability/sequential_design_10th_percentile"
orig_data_loc = "/Users/murph/dfnworks_variability/dfnworks_drivers"

## Baseline data
# orig_norm_data           = get_normalized_full_data_woOutliers(num_of_replications,
#                                                                data_loc = orig_data_loc,
#                                                                num_of_experiments, percentile)
# save(orig_norm_data, file = "sparse_lhs_normalized_data_woOutliers.Rdata")
# load("sparse_lhs_normalized_data_woOutliers.Rdata")

# Make sure you load this data in from the ES servers so that you're working with precisely the same data.
load("../sequential_design_10th_percentile/sparse_lhs_normalized_data.Rdata")
orig_norm_data = norm_data

# norm_data = orig_norm_data
# for(irep in 1:num_of_additional_input_locations){
#   # We begin by fitting a hetGP on this space and using the closed-form of the ISMPE
#   # to determine new input points.
#   new_job_num       = (total_num_from_initial_experiment +
#                          num_of_additional_replicates * irep)
# 
#   norm_data = get_updated_normalized_full_data_woOutliers(norm_data, (total_num_from_initial_experiment +
#                                                              num_of_additional_replicates * (irep-1)),
#                                                              num_of_additional_replicates,
#                                                              data_loc = models_loc,
#                                                              percentile = percentile)
# 
# }
# save(norm_data, file = "seq_design_data_woOutliers.Rdata")
load(file = "seq_design_data_woOutliers.Rdata")

# testing_data = get_updated_normalized_test_data_woOutliers(norm_data, test_data_loc,
#                                                             percentile = percentile)
# save(testing_data, file = "compiled_testing_data_woOutliers.Rdata")
load("compiled_testing_data_woOutliers.Rdata")

full_X_testing = NULL
for(val_i in 1:nrow(testing_data$X0)){
  temp_rows = matrix(rep(unlist(as.vector(testing_data$X0[val_i,])), times = testing_data$mult[val_i]),
                     ncol = ncol(testing_data$X0), nrow = testing_data$mult[val_i], byrow = TRUE)
  full_X_testing    = rbind(full_X_testing, temp_rows)
}
Z_test = testing_data$Z



######################################################
#######################################
##########################
#############
## From here, I should be able to return to the original GP analysis...
# I am no longer able to use the models that I fit at the original run of this experiment.
# They were run on a very poor approximation of the empirical PDF (rather than just
# taking a quantile).

# To perform the analysis that I want to do here (prior to re-running the Sequential
# Design), I'll need to essentially refit every GP.
# Kinda sucks.  Hopefully I can code this on the bus and run it overnight.

# Should be fine.

# We start by thus fitting the GP on the base lhs design data.
GP_analysis_data    = NULL
norm_data           = orig_norm_data
curr_min            = norm_data$new_min
curr_max            = norm_data$new_max
X0                  = as.matrix(norm_data$X0)
Z0                  = norm_data$Z0
mult                = norm_data$mult
Z                   = norm_data$Z
X                   = list(X0 = X0, Z0 = Z0, mult = mult)

load("../sequential_design_10th_percentile/first_model_in_design.Rdata")

# model = mleHetGP( X = list(X0 = X0, Z0 = Z0, mult = mult),
#                   noiseControl = list(g_min = 3),
#                   Z = Z, covtype = "Matern3_2",
#                   settings = list(checkHom = FALSE ) )
# Now I should have the model and the data from this point in the sequential design.
# Use these, with the test data, to calculate the metrics I am interested in.
full_X = NULL
for(val_i in 1:nrow(norm_data$X0)){
  temp_rows = matrix(rep(unlist(as.vector(norm_data$X0[val_i,])), times = norm_data$mult[val_i]),
                     ncol = ncol(norm_data$X0), nrow = norm_data$mult[val_i], byrow = TRUE)
  full_X    = rbind(full_X, temp_rows)
}
Z_train = norm_data$Z

sc.train.score = scores(model=model, Xtest = full_X, Ztest = Z_train)
sc.train.rmse  = scores(model=model, Xtest = full_X, Ztest = Z_train, return.rmse = TRUE)$rmse
sc.test.score  = scores(model=model, Xtest = full_X_testing, Ztest = Z_test)
sc.test.rmse   = scores(model=model, Xtest = full_X_testing, Ztest = Z_test, return.rmse = TRUE)$rmse

# According to the Diagnostics paper, I should be normalizing these errors according to the predicted
# variance at that point in the input space.
train_preds          = predict(object=model, x = full_X)
train_true           = as.matrix(Z_train)
train_pred_var       = train_preds$sd2 + train_preds$nugs
scaled_rmse_train    = sqrt( mean( (train_true - as.matrix(train_preds$mean))^2 / train_pred_var ) )

test_preds          = predict(object=model, x = full_X_testing)
test_true           = as.matrix(Z_test)
test_pred_var       = test_preds$sd2 + test_preds$nugs
scaled_rmse_test    = sqrt( mean( (test_true - as.matrix(test_preds$mean))^2 / test_pred_var ) )

# the optimal IMSPE at this location:
imspe_for_new_value = IMSPE_optim(model, h = 0)$value

# Collect data and record
temp_row = data.frame(score_train   = sc.train.score, rmse_train = sc.train.rmse, scaled_rmse_train = scaled_rmse_train,
                      score_test    = sc.test.score,  rmse_test  = sc.test.rmse, scaled_rmse_test = scaled_rmse_test,
                      imspe_new_val = imspe_for_new_value)

GP_analysis_data = rbind(GP_analysis_data, temp_row)

for(irep in 2:num_of_additional_input_locations) {
  
  model = NULL
  load(file=paste(models_loc, "/models/model_at_rep_", irep, ".Rdata", sep = ""))
  
  X0                  = as.matrix(norm_data$X0)
  Z0                  = norm_data$Z0
  mult                = norm_data$mult
  Z                   = norm_data$Z
  X                   = list(X0 = X0, Z0 = Z0, mult = mult)
  
  model = mleHetGP( X = list(X0 = X0, Z0 = Z0, mult = mult),
                    noiseControl = list(g_min = 3),
                    Z = Z, covtype = "Matern3_2",
                    settings = list(checkHom = FALSE ) )
  
  
  
  # This gets the new data and inputs (we already ran this in our sequential design).
  norm_data = get_updated_normalized_full_data_woOutliers(norm_data, (total_num_from_initial_experiment +
                                                             num_of_additional_replicates * (irep-1)),
                                               num_of_additional_replicates,
                                               models_loc, 
                                               percentile = percentile)
  print("the new inputs looks like.")
  print(tail(norm_data$full_inputs))
  
  # Now I should have the model and the data from this point in the (fixed) sequential design.
  # Use these, with the test data, to calculate the metrics I am interested in.
  full_X = NULL
  for(val_i in 1:nrow(norm_data$X0)){
    temp_rows = matrix(rep(unlist(as.vector(norm_data$X0[val_i,])), times = norm_data$mult[val_i]),
                       ncol = ncol(norm_data$X0), nrow = norm_data$mult[val_i], byrow = TRUE)
    full_X    = rbind(full_X, temp_rows)
  }
  Z_train = norm_data$Z
  
  sc.train.score = scores(model=model, Xtest = full_X, Ztest = Z_train)
  sc.train.rmse  = scores(model=model, Xtest = full_X, Ztest = Z_train, return.rmse = TRUE)$rmse
  sc.test.score  = scores(model=model, Xtest = full_X_testing, Ztest = Z_test)
  sc.test.rmse   = scores(model=model, Xtest = full_X_testing, Ztest = Z_test, return.rmse = TRUE)$rmse
  
  print(paste("rmse for this iteration is:", sc.test.rmse))
  
  # According to the Diagnostics paper, I should be normalizing these errors according to the predicted
  # variance at that point in the input space.
  train_preds          = predict(object=model, x = full_X)
  train_true           = as.matrix(Z_train)
  train_pred_var       = train_preds$sd2 + train_preds$nugs
  scaled_rmse_train    = sqrt( mean( (train_true - as.matrix(train_preds$mean))^2 / train_pred_var ) )
  
  test_preds          = predict(object=model, x = full_X_testing)
  test_true           = as.matrix(Z_test)
  test_pred_var       = test_preds$sd2 + test_preds$nugs
  scaled_rmse_test    = sqrt( mean( (test_true - as.matrix(test_preds$mean))^2 / test_pred_var ) )
  
  # the optimal IMSPE at this location:
  imspe_for_new_value = IMSPE_optim(model, h = 0)$value
  
  print(paste("new IMSPE for this iteration is:", imspe_for_new_value))
  
  # Collect data and record
  temp_row = data.frame(score_train   = sc.train.score, rmse_train = sc.train.rmse, scaled_rmse_train = scaled_rmse_train,
                        score_test    = sc.test.score,  rmse_test  = sc.test.rmse, scaled_rmse_test = scaled_rmse_test,
                        imspe_new_val = imspe_for_new_value)
  
  GP_analysis_data = rbind(GP_analysis_data, temp_row)
}
GP_analysis_data$seq_draw = 1:nrow(GP_analysis_data)
write.csv(GP_analysis_data, "GP_analysis_data_woOutliers.csv")


###############################################################################
###############################################################################
###############################################################################
## Visualization of the data.
analysis_data    = read.csv(file="GP_analysis_data_woOutliers.csv")
GP_analysis_data = analysis_data
analysis_data$X = NULL

########################################
## Graphing RMSE over sequential design.
n_rows = nrow(analysis_data)
rmse_data = data.frame(rmse = c(GP_analysis_data$rmse_train, GP_analysis_data$rmse_test),
                       data = c(rep("train",times=n_rows), rep("test",times=n_rows)),
                       seq_draw = c(GP_analysis_data$seq_draw, GP_analysis_data$seq_draw))
ggplot(rmse_data, aes(x = seq_draw, y = rmse, color = data)) + geom_line() + 
  ggtitle(TeX("")) + 
  theme_bw() + 
  # xlab(TeX("Scaled $\\alpha$")) + 
  # ylab(TeX("Pivoted Chol Error")) + 
  theme(axis.title = element_text(size = 15), plot.title = element_text(size = 15))

# That is kinda weird.  Let's look just at the testing data:
ggplot(analysis_data, aes(x = seq_draw, y = rmse_train)) + geom_line() + 
  ggtitle("RMSE on train set with additional Seq. Design points") + 
  theme_bw() + 
  xlab(TeX("Sequential Design Iteration")) + 
  ylab(TeX("RMSE")) + 
  theme(axis.title = element_text(size = 15), plot.title = element_text(size = 15))

ggplot(analysis_data, aes(x = seq_draw, y = rmse_test)) + geom_line() + 
  ggtitle("RMSE on test set with additional Seq. Design points") + 
  theme_bw() + 
  xlab(TeX("Sequential Design Iteration")) + 
  ylab(TeX("RMSE")) + 
  theme(axis.title = element_text(size = 15), plot.title = element_text(size = 15)) + xlim(10,50) + ylim(0.7,0.75)

# That is kinda weird.  Let's look just at the testing data:
ggplot(analysis_data, aes(x = seq_draw, y = scaled_rmse_train)) + geom_line() + 
  ggtitle("Scaled RMSE on train set with additional Seq. Design points") + 
  theme_bw() + 
  xlab(TeX("Sequential Design Iteration")) + 
  ylab(TeX("RMSE")) + 
  theme(axis.title = element_text(size = 15), plot.title = element_text(size = 15))

ggplot(analysis_data, aes(x = seq_draw, y = scaled_rmse_test)) + geom_line() + 
  ggtitle("Scaled RMSE on test set with additional Seq. Design points") + 
  theme_bw() + 
  xlab(TeX("Sequential Design Iteration")) + 
  ylab(TeX("RMSE")) + 
  theme(axis.title = element_text(size = 15), plot.title = element_text(size = 15))


########################################
## Graphing Score over sequential design.
ggplot(analysis_data, aes(x = seq_draw, y = score_train)) + geom_line() + 
  ggtitle("Score on train set with additional Seq. Design points")

ggplot(analysis_data, aes(x = seq_draw, y = score_test)) + geom_line() + 
  ggtitle("Score on test set with additional Seq. Design points")
# This last one looks somewhat good....but why is it so negative?

########################################
## Graphing Each IMPSE
ggplot(analysis_data, aes(x = seq_draw, y = imspe_new_val)) + geom_line() + 
  ggtitle("Model optimum IMSPE additional Seq. Design points")
# It would appear that we did, in truth, reduce this value as we calculated
# further and further design points.


# Note to discuss w Justin and Kelly -- love the scores and stuff.  I'm not liking
# the poor RMSE on the test set.  We may want to discuss this.







