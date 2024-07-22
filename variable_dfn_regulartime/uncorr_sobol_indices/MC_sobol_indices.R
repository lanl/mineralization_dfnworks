###############################################################################
### Script to perform Monte Carlo estimation for the Sobol Indices for each
### of our Gaussian processes (mean response and dispersion).
### Using Global SA: A Primer by A. Saltelli
## Author: Alexander C. Murph
## Date: August 2023
library(lhs)

estimate_sobol_indices = function(gp_model, num_samples, input_dimension, 
                                  var_of_response, var_names = NULL){
  # Note that this method assumes that inputs are scaled to be on the unit
  # hypercube.
  
  # Saltelli 2010
  # suggests quasi-random, so I'll use a LHS design. (the precise quasi-
  # random design used is not explicitly given NOR CITED UGH).  
  A_matrix = randomLHS(num_samples, input_dimension)
  B_matrix = randomLHS(num_samples, input_dimension)
  
  ############################################################
  ########## Sobol indices on the mean GP: ###################
  # Get the baseline predictions for these two random samples. 
  orig_response_A = predict(x = A_matrix, object = gp_model)
  orig_response_B = predict(x = B_matrix, object = gp_model)
  
  f_A_vector = orig_response_A$mean
  f_B_vector = orig_response_B$mean
  
  # Now, we must replace each column of A with a column of B, and recalculate
  # these predictions.  We'll hold these all in a list, since I believe that
  # we can use these to calculate both first-order and total-effect indices.
  f_A_sub_B_of_i     = list()
  # If I want all the same design points for these estimates, I must also do the sd2
  # matrix shuffling here (but I won't use this until later).
  f_A_sub_B_of_i_sd2 = list()
  for(i in 1:input_dimension){
    A_sub_B_of_i            = A_matrix
    # Saltelli also suggest a 'radial design' here, but I believe that this may
    # only be for total-effect.  When doing both, it seems more efficient to
    # do it the classical way.
    A_sub_B_of_i[,i]        = B_matrix[,i]
    # Note that every element of this list is a N vector of values.
    temp_orig_AB_model      = predict(x = A_sub_B_of_i, object = gp_model)
    f_A_sub_B_of_i[[i]]     = temp_orig_AB_model$mean
    f_A_sub_B_of_i_sd2[[i]] = temp_orig_AB_model$nugs
  }
  
  # With all these model evaluations, we calculate the estimates for the Sobol
  # indices.
  first_order  = list()
  total_effect = list()
  for(i in 1:input_dimension){
    for(j in 1:num_samples){
      if(j == 1){
        first_order[[i]]  = 0
        total_effect[[i]] = 0
      }
      jth_row_of_f_A    = f_A_vector
      # I've been having some issues with negative estimates for first-order.
      # I'm going to try a different estimation method:
      first_order[[i]]  = first_order[[i]] + (1/(num_samples*var_of_response)) *
                           f_B_vector[j] * (f_A_sub_B_of_i[[i]][j] - f_A_vector[j])
      # first_order[[i]]  = first_order[[i]] + (1/(2*num_samples*var_of_response)) *
      #                       (f_B_vector[j] - f_A_sub_B_of_i[[i]][j])**2
      total_effect[[i]] = total_effect[[i]] + 0.5 * (1/(num_samples*var_of_response)) * 
                           (f_A_vector[j] - f_A_sub_B_of_i[[i]][j]) ** 2
    }
    # first_order[[i]]  = 1 - first_order[[i]]
  }
  first_order_mean  = first_order
  total_effect_mean = total_effect
  
  ############################################################
  ########## Sobol indices on the sd2 GP: ###################
  first_order    = NULL
  total_effect   = NULL
  f_A_vector     = NULL
  f_B_vector     = NULL
  f_A_sub_B_of_i = NULL
  
  # Get the baseline predictions for these two random samples. 
  f_A_vector = orig_response_A$nugs
  f_B_vector = orig_response_B$nugs
  
  # Now, we must replace each column of A with a column of B, and recalculate
  # these predictions.  We'll hold these all in a list, since I believe that
  # we can use these to calculate both first-order and total-effect indices.
  # Recall that I did this for sd2 at the same time as for the mean.
  f_A_sub_B_of_i = f_A_sub_B_of_i_sd2
  
  # With all these model evaluations, we calculate the estimates for the Sobol
  # indices.
  first_order  = list()
  total_effect = list()
  for(i in 1:input_dimension){
    for(j in 1:num_samples){
      if(j == 1){
        first_order[[i]]  = 0
        total_effect[[i]] = 0
      }
      jth_row_of_f_A    = f_A_vector
      # I've been having some issues with negative estimates for first-order.
      # I'm going to try a different estimation method:
      # first_order[[i]]  = first_order[[i]] + (1/(num_samples*var_of_response)) *
      #                      f_B_vector[j] * (f_A_sub_B_of_i[[i]][j] - f_A_vector[j])
      first_order[[i]]  = first_order[[i]] + (1/(2*num_samples*var_of_response)) *
        (f_B_vector[j] - f_A_sub_B_of_i[[i]][j])**2
      total_effect[[i]] = total_effect[[i]] + 0.5 * (1/(num_samples*var_of_response)) * 
        (f_A_vector[j] - f_A_sub_B_of_i[[i]][j]) ** 2
    }
    first_order[[i]]  = 1 - first_order[[i]]
  }
  
  first_order_sd2  = first_order
  total_effect_sd2 = total_effect
  
  if(!is.null(var_names)){
    names(first_order_mean)  = var_names
    names(total_effect_mean) = var_names
    names(first_order_sd2)   = var_names
    names(total_effect_sd2)  = var_names
  }
  
  # Return all of the effects:
  return( list(first_order_mean = first_order_mean, total_effect_mean = total_effect_mean,
               first_order_sd2  = first_order_sd2,   total_effect_sd2 = total_effect_sd2) )
}


estimate_sobol_index = function(gp_model, num_samples, input_dimension, 
                                input_var_idx, var_of_response, var_names = NULL){
  # Note that this method assumes that inputs are scaled to be on the unit
  # hypercube.
  
  # Saltelli 2010
  # suggests quasi-random, so I'll use a LHS design. (the precise quasi-
  # random design used is not explicitly given NOR CITED UGH).  
  A_matrix = randomLHS(num_samples, input_dimension)
  B_matrix = randomLHS(num_samples, input_dimension)
  
  ############################################################
  ########## Sobol indices on the mean GP: ###################
  # Get the baseline predictions for these two random samples. 
  orig_response_A = predict(x = A_matrix, object = gp_model)
  orig_response_B = predict(x = B_matrix, object = gp_model)
  
  f_A_vector = orig_response_A$mean
  f_B_vector = orig_response_B$mean
  
  # Now, we must replace each column of A with a column of B, and recalculate
  # these predictions.  We'll hold these all in a list, since I believe that
  # we can use these to calculate both first-order and total-effect indices.
  f_A_sub_B_of_i     = list()
  # If I want all the same design points for these estimates, I must also do the sd2
  # matrix shuffling here (but I won't use this until later).
  f_A_sub_B_of_i_sd2 = list()
  
  # for(i in 1:input_dimension){
  i                       = input_var_idx
  A_sub_B_of_i            = A_matrix
  # Saltelli also suggest a 'radial design' here, but I believe that this may
  # only be for total-effect.  When doing both, it seems more efficient to
  # do it the classical way.
  A_sub_B_of_i[,i]        = B_matrix[,i]
  # Note that every element of this list is a N vector of values.
  temp_orig_AB_model      = predict(x = A_sub_B_of_i, object = gp_model)
  f_A_sub_B_of_i[[1]]     = temp_orig_AB_model$mean
  f_A_sub_B_of_i_sd2[[1]] = temp_orig_AB_model$nugs
  # }
  
  # With all these model evaluations, we calculate the estimates for the Sobol
  # indices.
  first_order  = list()
  total_effect = list()
  # for(i in 1:input_dimension){
  for(j in 1:num_samples){
    if(j == 1){
      first_order[[1]]  = 0
      total_effect[[1]] = 0
    }
    jth_row_of_f_A    = f_A_vector
    # first_order[[1]]  = first_order[[1]] + (1/num_samples) *
    #   f_B_vector[j] * (f_A_sub_B_of_i[[1]][j] - f_A_vector[j])
    first_order[[1]]  = first_order[[1]] + (1/(num_samples*var_of_response)) *
      f_B_vector[j] * (f_A_sub_B_of_i[[1]][j] - f_A_vector[j])
    total_effect[[1]] = total_effect[[1]] + 0.5 * (1/num_samples*var_of_response) * 
      (f_A_vector[j] - f_A_sub_B_of_i[[1]][j]) ** 2
  }
  # }
  first_order_mean  = first_order
  total_effect_mean = total_effect
  
  ############################################################
  ########## Sobol indices on the sd2 GP: ###################
  first_order    = NULL
  total_effect   = NULL
  f_A_vector     = NULL
  f_B_vector     = NULL
  f_A_sub_B_of_i = NULL
  
  # Get the baseline predictions for these two random samples. 
  f_A_vector = orig_response_A$nugs
  f_B_vector = orig_response_B$nugs
  
  # Now, we must replace each column of A with a column of B, and recalculate
  # these predictions.  We'll hold these all in a list, since I believe that
  # we can use these to calculate both first-order and total-effect indices.
  # Recall that I did this for sd2 at the same time as for the mean.
  f_A_sub_B_of_i = f_A_sub_B_of_i_sd2
  
  # With all these model evaluations, we calculate the estimates for the Sobol
  # indices.
  first_order  = list()
  total_effect = list()
  # for(i in 1:input_dimension){
  for(j in 1:num_samples){
    if(j == 1){
      first_order[[1]]  = 0
      total_effect[[1]] = 0
    }
    jth_row_of_f_A    = f_A_vector
    # first_order[[1]]  = first_order[[1]] + (1/num_samples) *
    #   f_B_vector[j] * (f_A_sub_B_of_i[[1]][j] - f_A_vector[j])
    # total_effect[[1]] = total_effect[[1]] + 0.5 * (1/num_samples) * 
    #   (f_A_vector[j] - f_A_sub_B_of_i[[1]][j]) ** 2
    
    first_order[[1]]  = first_order[[1]] + (1/(2*num_samples*var_of_response)) *
      (f_B_vector[j] - f_A_sub_B_of_i[[1]][j])**2
    total_effect[[1]] = total_effect[[1]] + 0.5 * (1/(num_samples*var_of_response)) * 
      (f_A_vector[j] - f_A_sub_B_of_i[[1]][j]) ** 2
  
    }
  # }
  
  first_order_sd2  = first_order
  total_effect_sd2 = total_effect
  
  if(!is.null(var_names)){
    names(first_order_mean)  = var_names[input_var_idx]
    names(total_effect_mean) = var_names[input_var_idx]
    names(first_order_sd2)   = var_names[input_var_idx]
    names(total_effect_sd2)  = var_names[input_var_idx]
  }
  
  # Return all of the effects:
  return( list(first_order_mean = first_order_mean, total_effect_mean = total_effect_mean,
               first_order_sd2  = first_order_sd2,   total_effect_sd2 = total_effect_sd2) )
}





