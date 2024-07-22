# Exploratory Data Analysis on the Chemistry data.  I did some of this 
# already in the sequential_draws.R file; I'm going through it more formally 
# here.
## Author: Alexander C. Murph
## Date: January 2024
library(hetGP)
library(gridExtra)
library(ggplot2)
library(reshape)
library(tidyverse)
setwd("~/GitLab/sa_for_chemistry/variable_dfn/max_values/eda")

# Data location:
data_path = "~/GitLab/sa_for_chemistry/variable_dfn/max_values/data"

# Load data:
norm_data    = read.csv(paste(data_path, "/combined_data.csv", sep = ""))
norm_data$X  = NULL

# Variable info:
new_names = c("dfn_seed",                   "curr_inflow_rate",           "curr_diffusion_coef",        "curr_porosity",              "curr_gypsum_rate_constant", 
              "curr_calcite_rate_constant", "curr_gypsum_surface_area",   "curr_calcite_surface_area",  "dfn_p32",                    "dfn_volume",                
              "backbone_volume",            "backbone_p32",               "final_gypsum",               "final_calcite",              "Pe",                        
              "Da_1_gypsum",                "Da_2_gypsum",                "Da_1_calcite",               "Da_2_calcite",               "tau",                       
              "max_calcite",                "min_gypsum")

# Variable updates after the EDA.
# taopoints_to_remove = c("calcite_flush_0_nondim", "calcite_flush_1_nondim", 
#                         "calcite_flush_2_nondim", "calcite_flush_3_nondim",
#                         "gypsum_flush_0_nondim", "gypsum_flush_1_nondim",
#                         "gypsum_flush_2_nondim", "gypsum_flush_3_nondim")
# # early_gypsum_points = c("gypsum_flush_4_nondim")
# new_names           = new_names[which(!(new_names%in%taopoints_to_remove))]

cols_of_log_vars                 = c(2,3,5,6,8)
cols_of_net_vars                 = c(9,10,11,12)
cols_of_non_fctn_vars            = c(2:14)
cols_of_fctn_vars                = c(15,16,17,18,19,20)
output_vars                      = c(21:length(new_names))
input_vars_for_emulation_nonfctn = c(2:8,10)
input_vars_for_emulation_fctn    = c(cols_of_net_vars[2], cols_of_fctn_vars)


### Justification for putting the functional variables on the log scale.
## lol just look at them:
pairs(norm_data[,input_vars_for_emulation_fctn])
# All the values for the variables other than dfn_p32 here are heavily concentrated towards zero.
# EDA on the functional variables shows that they should really be considered on the log-scale.
# This will lead to a change needed in the sobol calculation.  I need to get the explicit
# functions used from Jeffrey.
norm_data[,input_vars_for_emulation_fctn] = log(norm_data[,input_vars_for_emulation_fctn])
# and now they look a lot better:
pairs(norm_data[,input_vars_for_emulation_fctn])

# From here, we rescale.  A 0-1 scale will not affect the 
# sobol calculations later.  We then split and save the data.
# zero_one         = function(x){(x-min(x, na.rm = T))/(max(x, na.rm = T)-min(x, na.rm = T))}
# norm_data        = as.data.frame(apply(norm_data, 2, zero_one))

##########
## Correlated Inputs (non-ftcn vars)
# The inputted variables are meant to be uncorrelated.  This means that I'll need to
# handle the functional inputs separately.  The input variables that we sample
# via a LHS are made to be uncorrelated.  The remaining question is are there
# (original) input variables that are correlated with the graph-measured vars.
pairs(norm_data[,cols_of_non_fctn_vars])
# Inflow rate seems to drive final gypsum and final calcite.  
# The graph-measured variables are all strongly positively correlated with
# one another.
# I am going to take only 1 graph-measured variable (dfn_volume).  I am going
# to remove final_calcite and final_gypsum.
# Note that I'm taking p32 out of this for the fixed-DFN experiment.
input_vars_for_emulation_nonfctn = c(2:8)

##########
## Correlated Inputs (w/  nonftcn vars)
pairs(norm_data[,input_vars_for_emulation_nonfctn])
# Hmmm...I think I need to get these relationships explicitly from Jeffrey.
# Also everything but dfn_p32 seems like it needs to be on the log-scale before 
# we emulate.


######
# Check for heteroskedasaticity with the QoIs.
# Let's look over the non-functional variables first.
lst_p = NULL
count = 1
for(var_idx in input_vars_for_emulation_nonfctn){
  for(output_var_idx in output_vars){
    lst_p[[count]] = ggplot(norm_data, aes_string(x = new_names[var_idx], y = new_names[output_var_idx])) + geom_point() + ggtitle(paste("Values of QoI", new_names[output_var_idx]))
    count          = count + 1
  }
}

glist = lapply(lst_p, ggplotGrob)
pdf("heterskedastic_analysis_nonfctn.pdf")
marrangeGrob(grobs=glist, nrow=2, ncol=2)
dev.off()

# Now for the fctn vars
lst_p = NULL
count = 1
for(var_idx in input_vars_for_emulation_fctn){
  for(output_var_idx in output_vars){
    lst_p[[count]] = ggplot(norm_data, aes_string(x = new_names[var_idx], y = new_names[output_var_idx])) + 
                      geom_point() + ggtitle(paste("Values of QoI", new_names[output_var_idx]))
    count          = count + 1
  }
}

glist = lapply(lst_p, ggplotGrob)
pdf("heterskedastic_analysis_fctn.pdf")
marrangeGrob(grobs=glist, nrow=2, ncol=2)
dev.off()


###### Correlations with outputs:
pairs(norm_data[,c(input_vars_for_emulation_fctn, output_vars[11:20])])
pairs(norm_data[,c(input_vars_for_emulation_nonfctn, output_vars[1:8])])


#### Let's visualize the outputs.  I'm seeing some weird stuff on the earlier end
#### for these values on the tao time scale.
# temp_seed        = norm_data$dfn_seed
norm_data_temp   = norm_data
# zero_one         = function(x){(x-min(x, na.rm = T))/(max(x, na.rm = T)-min(x, na.rm = T))}
# norm_data        = apply(norm_data, 2, zero_one)
gypsum_data      = as.data.frame(norm_data[,gypsum_outputs])
xx               = melt(gypsum_data)
# names(xx) = c('idx','variable','value')
ggplot(xx, aes(x = factor(variable), y = value)) + geom_boxplot()+
  theme(axis.text.x = element_text(angle = 270, vjust = 0.5, hjust=1))

gypsum_data$sim_num     = factor(1:nrow(gypsum_data))
yy                      = melt(gypsum_data, id = "sim_num")
yy$variable             = factor(yy$variable)
ggplot(yy, aes(x = factor(variable), y = value, color = sim_num, group = sim_num)) + geom_line()+
  theme(legend.position = 'none',axis.text.x = element_text(angle = 270, vjust = 0.5, hjust=1))




calcite_outputs_modified = c("calcite_flush_2_nondim", "calcite_flush_3_nondim", "calcite_flush_4_nondim",
                             "calcite_flush_5_nondim",     "calcite_flush_6_nondim",
                             "calcite_flush_7_nondim",    
                             "calcite_flush_8_nondim",     "calcite_flush_9_nondim",      
                             "calcite_flush_10_nondim",    "calcite_flush_11_nondim",  "calcite_flush_12_nondim")

calcite_data = as.data.frame(norm_data[,calcite_outputs_modified])
xx = melt(calcite_data)
ggplot(xx, aes(x = factor(variable), y = value)) + geom_boxplot()+
  theme(axis.text.x = element_text(angle = 270, vjust = 0.5, hjust=1))

calcite_data$sim_num = factor(1:nrow(calcite_data))
yy                      = melt(calcite_data, id = "sim_num")
yy$variable             = factor(yy$variable)
ggplot(yy, aes(x = factor(variable), y = value, color = sim_num, group = sim_num)) + geom_line()+
  theme(legend.position = 'none',axis.text.x = element_text(angle = 270, vjust = 0.5, hjust=1))

# How much does each experiment vary?
zz = c()
for(i in 1:nrow(calcite_data)){
  zz = c(zz, sd(yy[which( yy$sim_num == i),3], na.rm = T ))
}


# Okay so big thought about this experiment.  Maybe what we actually care about is
# 1) the CHANGE in total gypsum/calcite that happens at each time step.
# 2) normalizing this change by the total....?
gypsum_data      = as.data.frame(norm_data[,gypsum_outputs])
final_gypsum     = norm_data[,"final_gypsum"]
# This is often zero....why might gypsum go down?

f = function(x){max(x,na.rm=T)}
max_values = as.numeric(apply(gypsum_data, 1, f))
ff = function(x){
  x / max_values
  }
for(col_idx in 1:ncol(gypsum_data)){
  gypsum_data[,col_idx] = gypsum_data[,col_idx] / max_values
}
xx               = melt(gypsum_data)

gypsum_data$sim_num     = factor(1:nrow(gypsum_data))
yy                      = melt(gypsum_data, id = "sim_num")
yy$variable             = factor(yy$variable)
ggplot(yy, aes(x = factor(variable), y = value, color = sim_num, group = sim_num)) + geom_line()+
  theme(legend.position = 'none',axis.text.x = element_text(angle = 270, vjust = 0.5, hjust=1))


## Let's try this for calcite:
calcite_data      = as.data.frame(norm_data[,calcite_outputs])
final_calcite     = norm_data[,"final_calcite"]
# This is often zero....why might calcite go down?

f = function(x){max(x,na.rm=T)}
max_values = as.numeric(apply(calcite_data, 1, f))
# max_values = final_calcite
ff = function(x){
  x / max_values
}
for(col_idx in 1:ncol(calcite_data)){
  calcite_data[,col_idx] = calcite_data[,col_idx] / max_values
}
xx               = melt(calcite_data)
ggplot(xx, aes(x = factor(variable), y = value)) + geom_boxplot()+
  theme(axis.text.x = element_text(angle = 270, vjust = 0.5, hjust=1))

calcite_data$sim_num     = factor(1:nrow(calcite_data))
yy                      = melt(calcite_data, id = "sim_num")
yy$variable             = factor(yy$variable)
ggplot(yy, aes(x = factor(variable), y = value, color = sim_num, group = sim_num)) + geom_line()+
  theme(legend.position = 'none',axis.text.x = element_text(angle = 270, vjust = 0.5, hjust=1))


# OKAY game plan: let's focus on calcite and do this normalization by maximum term idea.






