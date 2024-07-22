# Helper functions for the sequential design.
# This is an updated code suite where I update some issues I had created
# by using bkde unnecessarily.  For this file, I just use the quantile
# function with 'percentile.'
# Since I'm working with the quantile function, I'm also going to stop
# struggling over this 0-1 scale.
## Author: Alexander Murph
## Fate: August 2023
library(data.table)
library(KernSmooth)
library(pracma)

min_alpha_semicorr = 1*10**-10
max_alpha_semicorr = 1*10**-8
min_beta_semicorr  = 0.1
max_beta_semicorr  = 1.2
min_sigma_semicorr = 0.5
max_sigma_semicorr = 1
min_alpha_TPL_exp  = 2.25
max_alpha_TPL_exp  = 3.5
min_p32            = 5e-02
max_p32            = 0.2

# Convert pdf to cdf
pdf_to_cdf <- function(pdf, x_grid, norm=TRUE){
  cdf = cumtrapz(x_grid, pdf)
  if(norm){ 
    cdf = cdf/cdf[length(cdf)] 
  }
  return(cdf)
}

# Convert cdf to pdf
cdf_to_pdf <- function(cdf, x_grid){
  pdf = gradient(as.numeric(cdf), x_grid)
  return(pdf)
}

# Convert cdf to quantile function
# Uses linear interpolation
cdf_to_quant <- function(cdf, x_grid){
  interp_fn = approxfun(cdf, x_grid)
  # quant = interp_fn(x_grid)
  quant = interp_fn
  return(quant)
}

get_normalized_full_data = function(num_of_replications, num_of_experiments,
                                    percentile = NULL){
  # The data are normalized to a 0-1 scale over the ENTIRE simulation.  Because
  # of this, whenever we add a new point that has a max/min BT time above/below
  # the current spread, we must reload the raw data and re-normalize everything.
  ## Note: I could do a slight speedup by not reloading from the files every time,
  ## but it seems that fread is so fast that this might not be necessary.
  
  # `percentile` is the percentile of the BT pdf to report.  If this remains
  # NULL, we return the entire pdf.
  
  # We begin by collecting all of the data values into a single DF.  
  dfn_HF_data = NULL
  print(paste("iterating through all", (num_of_replications*num_of_experiments), "of the dfnWorks log files..."))
  for(ijob in 0:(num_of_replications*num_of_experiments-1)){
    exp_num       = ijob%%num_of_experiments
    cwd           = getwd()
    base_dir      = cwd
    data_location = paste(base_dir, '/../dfnworks_drivers/compiled_data/output_', exp_num, '/', sep = "")
    file_name     = paste(data_location, 'partime', ijob, '.dat', sep = "")
    if(file.exists(file_name)){
      print(paste("file number", ijob, "exists.  Reading in the data..."))
      HF_flowdata                       = fread(file_name, fill=TRUE)
      HF_flowdata[,9:ncol(HF_flowdata)] = NULL
      names(HF_flowdata)                = c("num_of_time_steps", "flux_weights", 
                                            "total_travel_time", "x_", "y_", 
                                            "z_final_pos", "beta", "total_length_m" )
      
      temp_row      = data.frame(sim_num            = rep(ijob, times = nrow(HF_flowdata)), 
                                 total_travel_time  = as.numeric(HF_flowdata$total_travel_time))
      temp_row      = temp_row[complete.cases(temp_row),]
      dfn_HF_data   = rbind(dfn_HF_data, temp_row)
    }
  }
  
  # I need to normalize over the entire simulation, which means I'll need to 
  # save these values and potentially update this normalization as I get new
  # values (which is a huge pain and hopefully not too computationally cumbersome).
  min_BT = min(dfn_HF_data$total_travel_time)
  max_BT = max(dfn_HF_data$total_travel_time)
  
  # Now, we gather the input location for every raw BT, and transform this raw BT
  # to a PDF.
  sim_inputs         = NULL
  if(is.null(percentile)){
    sim_response       = NULL
  } else {
    sim_response       = c()
  }
 
  num_of_grid_points = 50
  bw                 = 0.03
  print(paste("gathering all", length(unique(dfn_HF_data$sim_num)), "of the unique inputs"))
  for(sim_num in unique(dfn_HF_data$sim_num)){
    exp_num       = sim_num%%num_of_experiments
    cwd           = getwd()
    base_dir      = cwd
    data_location = paste(base_dir, '/../dfnworks_drivers/compiled_data/output_', exp_num, '/', sep = "")
    
    # First, we'll grab the response:
    temp_data     = NULL
    temp_data     = dfn_HF_data[which(dfn_HF_data$sim_num == sim_num),]
    if(nrow(temp_data)==0){
      next
    }

    if(is.null(percentile)){
      pdf_feats     = bkde(temp_data$total_travel_time, kernel="normal", bandwidth=bw,
                           gridsize=num_of_grid_points, range.x=c(0,1))
      pdf_feats     = pdf_feats$y
      sim_response  = rbind(sim_response, pdf_feats)
    } else {
      pdf_feats     = quantile(log(temp_data$total_travel_time), probs = percentile)
      # cdf           = pdf_to_cdf(pdf_feats$y, pdf_feats$x)
      # quant         = cdf_to_quant(cdf, pdf_feats$x)
      # pdf_feats     = quant(percentile)
      sim_response  = c(sim_response, pdf_feats)
    }
    
    # Now grab the inputs:
    file_name     = paste(data_location, 'SA_params_', sim_num, '.dat', sep = "")
    param_data    = read.csv(file_name)
    print("Checking to see if the parameter data added to the full data is correct.  Raw SA_params file inputs are:")    
    print(param_data)

    # Big note: My original coding of how these variables are placed in the input file is different than later
    #           versions.  This is definitely going to lead to issues...
    temp_row      = data.frame(alpha_semi = round((param_data$alpha_semicorr[1] - min_alpha_semicorr) /
                                            (max_alpha_semicorr - min_alpha_semicorr), digits = 15), 
                               beta_semi  = round((param_data$beta_semicorr[1]- min_beta_semicorr) /
                                            (max_beta_semicorr - min_beta_semicorr), digits = 15),
                               sigma_semi = round((param_data$sigma_semicorr[1]- min_sigma_semicorr) /
                                            (max_sigma_semicorr - min_sigma_semicorr), digits = 15), 
                               alpha_TPL  = round((param_data$alpha_radius[1] - min_alpha_TPL_exp) /
                                            (max_alpha_TPL_exp - min_alpha_TPL_exp), digits = 15),
                               p32        = round((param_data$p32[1] - min_p32) / 
                                            (max_p32 - min_p32), digits = 15),
   			                       sim_num    = sim_num )
    print("temp row added to data is:")
    print(temp_row)
    sim_inputs    = rbind(sim_inputs, temp_row)
  }
  
  # Make sure all our data objects are the correct size.
  if(is.null(percentile)){
    if(nrow(sim_inputs)!=nrow(sim_response)) stop("All inputs must map to a unique output.")
  } else {
    if(nrow(sim_inputs)!=length(sim_response)) stop("All inputs must map to a unique output.")
  }
  
  print("creating indicator variable for the inputs...")  
  # Lastly, in the case where we have replicates, hetGP requires a specific form
  # for the input.  We'll do that transformation here.
  concat_values  = function(x){paste(x, collapse = ", ")}
  input_ind      = apply(sim_inputs[,1:5], 1, concat_values)
  sim_inputs$ind = input_ind
  # The following variable names are based on the hetGP package.
  # Note that hetGP assumes a 1-dim output (insists on it being a vector).
  # TODO: update this to handle the full PDF someday.
  X0            = NULL
  mult          = c()
  Z             = c()
  Z0            = c()
  for(unique_input in unique(input_ind)){
    temp_inputs    = sim_inputs[which(input_ind == unique_input), 1:5]
    mult           = c(mult, nrow(temp_inputs))
    X0             = rbind(X0, temp_inputs[1,])
    if(is.null(percentile)){
      temp_responses = sim_response[which(input_ind == unique_input),]
      Z              = rbind(Z, temp_responses)
      Z0             = rbind(Z, colMeans(temp_responses))
    } else {
      temp_responses = sim_response[which(input_ind == unique_input)]
      Z              = c(Z, temp_responses)
      Z0             = c(Z0, mean(temp_responses))
    }
  }
  
  return(list(full_inputs  = sim_inputs, full_outputs = sim_response,
              new_min      = min_BT,     new_max      = max_BT,
              X0           = X0,         Z0           = Z0,
              mult         = mult,       Z            = Z,
              raw_sim_data = dfn_HF_data))
}

get_updated_normalized_full_data = function(normalized_curr_data, 
                                            num_of_experiments_finished,
                                            num_of_new_replicates,
                                            percentile = NULL){
  # See get_normalized_full_data for further details.  The 0-1 scale may need
  # to be updated depending on the new values.
  
  # Start by collecting SPECIFICALLY the new data values:
  dfn_HF_newdata = NULL
  cat("attempting to add new data points...\n", file = "full_output.txt", append = TRUE)
  for(ijob in (num_of_experiments_finished+5):(num_of_experiments_finished+num_of_new_replicates-1+5)){
    exp_num       = ijob
    cwd           = getwd()
    base_dir      = cwd
    data_location = paste(base_dir, '/compiled_data/output_', exp_num, '/', sep = "")
    
    file_name       = paste(data_location, 'partime', ijob, '.dat', sep = "")
    if(file.exists(file_name)){
      HF_flowdata                       = fread(file_name, fill=TRUE)
      HF_flowdata[,9:ncol(HF_flowdata)] = NULL
      names(HF_flowdata)                = c("num_of_time_steps", "flux_weights", 
                                            "total_travel_time", "x_", "y_", 
                                            "z_final_pos", "beta", "total_length_m" )
      
      temp_row         = data.frame(sim_num            = rep(ijob, times = nrow(HF_flowdata)), 
                                    total_travel_time  = as.numeric(HF_flowdata$total_travel_time))
      temp_row         = temp_row[complete.cases(temp_row),]
      dfn_HF_newdata   = rbind(dfn_HF_newdata, temp_row)
      cat("succeeded in adding a new data point! \n", file = "full_output.txt", append = TRUE)
    }
  }

  if(is.null(dfn_HF_newdata)){
	cat("no new data were created. \n", file = "full_output.txt", append = TRUE)
	return(list(full_inputs           = normalized_curr_data$full_inputs, full_outputs = normalized_curr_data$full_outputs,
        	    new_min               = normalized_curr_data$new_min,     new_max      = normalized_curr_data$new_max,
             	    X0                    = normalized_curr_data$X0,         Z0           = normalized_curr_data$Z0,
             	    mult                  = normalized_curr_data$mult,       Z            = normalized_curr_data$Z,
             	    raw_sim_data          = normalized_curr_data$raw_sim_data,
             	    updated_normalization = NULL,
             	    X_new                 = NULL,      Z_new        = NULL))
  }
  
  # I need to see if the new data values require a re-normalization.
  min_BT      = min(dfn_HF_newdata$total_travel_time)
  max_BT      = max(dfn_HF_newdata$total_travel_time)
  dfn_HF_data = NULL
  
  dfn_HF_data                      = rbind(normalized_curr_data$raw_sim_data,
                                           dfn_HF_newdata)
  updated_normalization            = FALSE

  cat("updated the data normalization \n", file = "full_output.txt", append = TRUE)
  
  # TODO: I don't really have to recalculate the pdfs unless I renormalize, 
  #       I'm just being lazy here. if I have a lot of runtime issues this 
  #       would be a place to fix that.
  # Now, we gather the input location for every raw BT, and transform this raw BT
  # to a PDF.
  sim_inputs         = NULL
  if(is.null(percentile)){
    sim_response       = NULL
  } else {
    sim_response       = c()
  }
  
  cat("recalculating pdfs...\n", file = "full_output.txt", append = TRUE)
  Z_new              = c()
  X_new              = NULL
  num_of_grid_points = 50
  bw                 = 0.03
  for(sim_num in unique(dfn_HF_data$sim_num)){
    exp_num       = sim_num%%num_of_experiments
    cwd           = getwd()
    base_dir      = cwd

    # First, we'll grab the response:
    temp_data     = NULL
    temp_data     = dfn_HF_data[which(dfn_HF_data$sim_num == sim_num),]
    if(nrow(temp_data)==0){
      next
    }
    
    if(is.null(percentile)){
      pdf_feats     = bkde(temp_data$total_travel_time, kernel="normal", bandwidth=bw,
                           gridsize=num_of_grid_points, range.x=c(0,1))
      pdf_feats     = pdf_feats$y
      sim_response  = rbind(sim_response, pdf_feats)
    } else {
      pdf_feats     = quantile(log(temp_data$total_travel_time), probs = percentile)
      sim_response  = c(sim_response, pdf_feats)
    }

    # The fast update to the hetGP requires JUST the new inputs and new reponses.  I'll record
    # this here and return it in the full data list.
    full_input_data = norm_data$full_inputs
    if( !(sim_num %in% full_input_data$sim_num) ){
    	# Note that this assumes that we have a percentile that we're using
    	# TODO: update this for functional output.
    	Z_new = c(Z_new, pdf_feats)
    }
    
    # Now grab the inputs:
    # TODO: This is another thing that doesn't need to happen again.  I want to 
    #       debug other things, but maybe return to this.
    if(sim_num %in% full_input_data$sim_num){
	      temp_row     = full_input_data[which(full_input_data$sim_num == sim_num),]
        temp_row$ind = NULL
    } else {
	      data_location = paste(base_dir, '/compiled_data/output_', sim_num, '/', sep = "")
	      file_name     = paste(data_location, 'SA_params_', sim_num, '.dat', sep = "")
        param_data    = read.csv(file_name)

	#cat("Checking to see if the parameter data added to the full data is correct.  Raw SA_params file inputs are:\n", file = "full_output.txt", append = TRUE)
    	#cat(param_data, file = "full_output.txt", append = TRUE)
	temp_row      = data.frame(alpha_semi = round((param_data$alpha_semicorr[1] - min_alpha_semicorr) /
                                        (max_alpha_semicorr - min_alpha_semicorr), digits = 15),
                           beta_semi  = round((param_data$beta_semicorr[1]- min_beta_semicorr) /
                                        (max_beta_semicorr - min_beta_semicorr), digits = 15),
                           sigma_semi = round((param_data$sigma_semicorr[1]- min_sigma_semicorr) /
                                        (max_sigma_semicorr - min_sigma_semicorr), digits = 15),
                           alpha_TPL  = round((param_data$alpha_radius[1] - min_alpha_TPL_exp) /
                                        (max_alpha_TPL_exp - min_alpha_TPL_exp), digits = 15),
                           p32        = round((param_data$p32[1] - min_p32) /
                                        (max_p32 - min_p32), digits = 15),
                           sim_num    = sim_num )

       # temp_row      = data.frame(alpha_semi = round((param_data[1] - min_alpha_semicorr) /
        #                             (max_alpha_semicorr - min_alpha_semicorr), digits = 15),
         #                          beta_semi  = round((param_data[1,2]- min_beta_semicorr) /
          #                           (max_beta_semicorr - min_beta_semicorr), digits = 15),
           #                        sigma_semi = round((param_data[1,3]- min_sigma_semicorr) /
            #                         (max_sigma_semicorr - min_sigma_semicorr), digits = 15),
             #                      alpha_TPL  = round((param_data[1,4] - min_alpha_TPL_exp) /
              #                       (max_alpha_TPL_exp - min_alpha_TPL_exp), digits = 15),
               #                    p32        = round((param_data[1,5] - min_p32) /
                #                     (max_p32 - min_p32), digits = 15),
       	#			                     sim_num = sim_num)
	      X_new         = rbind(X_new, temp_row[,1:5])
    	#cat("temp row added to data is:\n", file = "full_output.txt", append = TRUE)
    	#cat(temp_row, file = "full_output.txt", append = TRUE)
    }
    sim_inputs    = rbind(sim_inputs, temp_row)
  }

  cat("finished recalculating pdfs + adding new input location logs \n", file = "full_output.txt", append = TRUE)
  
  # Make sure all our data objects are the correct size.
  if(is.null(percentile)){
    if(nrow(sim_inputs)!=nrow(sim_response)) stop("All inputs must map to a unique output.")
  } else {
    if(nrow(sim_inputs)!=length(sim_response)) stop("All inputs must map to a unique output.")
  }
  
  # Lastly, in the case where we have replicates, hetGP requires a specific form
  # for the input.  We'll do that transformation here.
  concat_values  = function(x){paste(x, collapse = ", ")}
  input_ind      = apply(sim_inputs[,1:5], 1, concat_values)
  sim_inputs$ind = input_ind
  # The following variable names are based on the hetGP package.
  # Note that hetGP assumes a 1-dim output (insists on it being a vector).
  # TODO: update this to handle the full PDF someday.
  X0            = NULL
  mult          = c()
  Z             = c()
  Z0            = c()
  for(unique_input in unique(input_ind)){
    temp_inputs    = sim_inputs[which(input_ind == unique_input), 1:5]
    mult           = c(mult, nrow(temp_inputs))
    X0             = rbind(X0, temp_inputs[1,])
    if(is.null(percentile)){
      temp_responses = sim_response[which(input_ind == unique_input),]
      Z              = rbind(Z, temp_responses)
      Z0             = rbind(Z, colMeans(temp_responses))
    } else {
      temp_responses = sim_response[which(input_ind == unique_input)]
      Z              = c(Z, temp_responses)
      Z0             = c(Z0, mean(temp_responses))
    }
  }
  
  return(list(full_inputs           = sim_inputs, full_outputs = sim_response,
              new_min               = min_BT,     new_max      = max_BT,
              X0                    = X0,         Z0           = Z0,
              mult                  = mult,       Z            = Z,
              raw_sim_data          = dfn_HF_data, 
              updated_normalization = updated_normalization,
	            X_new                 = X_new,      Z_new        = Z_new))
}



