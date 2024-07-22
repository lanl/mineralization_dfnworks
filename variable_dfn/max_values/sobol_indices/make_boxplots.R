###############################################################################
### Script to make boxplots from bootstrap Sobol' indices data.
## Author: Alexander C. Murph
## Date: January 2024

## Let's make some box-plots :)
library(ggplot2)
library(reshape)
library(latex2exp)
library(gridExtra)
library(grid)
library(tidyverse)
library(patchwork)
setwd("~/GitLab/sa_for_chemistry/variable_dfn/max_values/sobol_indices")

###############################################################################
###############################################################################
###############################################################################
####### Graphs for the non-functional variables.
sobol_data_for_nonfctn = NULL
seed_data_for_nonfctn  = NULL
for(file_name in list.files('bootstrap_data_nonfctn')){
  temp_data              = read.csv(paste('bootstrap_data_nonfctn/', file_name, sep = ''))
  temp_data$X            = NULL
  
  # if((file_name=="calcite_flush_2_nondim.csv")|(file_name=="seed_variable__calcite_flush_2_nondim__.csv")|
  #    (file_name=="seed_variable__calcite_flush_4_nondim__1.csv") | (file_name=="calcite_flush_4_nondim_1.csv")) next
  
  split_string = unlist(strsplit(file_name, '__'))
  if(split_string[1] == 'seed_variable'){
    temp_data$output_var   = split_string[2]
    seed_data_for_nonfctn  = rbind(seed_data_for_nonfctn, temp_data)
  } else {
    sobol_data_for_nonfctn = rbind(sobol_data_for_nonfctn, temp_data)
  }
}

new_names = c("dfn_seed",                   "curr_inflow_rate",           "curr_diffusion_coef",        "curr_porosity",              "curr_gypsum_rate_constant", 
              "curr_calcite_rate_constant", "curr_gypsum_surface_area",   "curr_calcite_surface_area",  "dfn_p32",                    "dfn_volume",                
              "backbone_volume",            "backbone_p32",               "final_gypsum",               "final_calcite",              "Pe",                        
              "Da_1_gypsum",                "Da_2_gypsum",                "Da_1_calcite",               "Da_2_calcite",               "tau",                       
              "max_calcite",                "min_gypsum")


gypsum_outputs  = c("min_gypsum")
outputs_abrv    = c("flush_0", "flush_1",
                    "flush_2", "flush_3", "flush_4",     
                    "flush_5", "flush_6", 
                    "flush_7", "flush_8", "flush_9",     
                    "flush_10", "flush_11",   
                    "flush_12")
calcite_outputs = c("max_calcite")

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# First Order Indices on the Mean Response GP, for calcite variables:
bootstrapped_FOmean            = sobol_data_for_nonfctn[which((sobol_data_for_nonfctn$sobol_type == 'FOmean')&(sobol_data_for_nonfctn$output_var %in% calcite_outputs)),]
bootstrapped_FOmean$boot_num   = NULL
bootstrapped_FOmean$sobol_type = NULL
graph_data                     = melt(bootstrapped_FOmean, id. = "output_var")

g2 = ggplot(graph_data, aes(x = variable, y = value, fill = output_var)) + geom_boxplot(lwd=0.25, outlier.size = 0.2, width = 0.9) + 
  ggtitle(TeX("")) +theme_bw()+ 
  xlab("") + 
  ylab("") +
  theme(axis.text=element_text(size=15), axis.title = element_text(size = 15)) + 
  theme(plot.title = element_blank(),
        axis.text.x = element_text(angle = 270, vjust = 0.5, hjust=1),
        panel.spacing = unit(0, "cm"),
        plot.margin = margin(0, 0, 0, 0, "cm"), 
        plot.caption = element_blank()) + labs(fill="Flush Time") + guides(fill="none")  + scale_fill_manual(name = "max_calcite",values = c("grey"))


# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Total-Effect Indices on the Mean Response GP, for calcite variables:
bootstrapped_TEmean            = sobol_data_for_nonfctn[which((sobol_data_for_nonfctn$sobol_type == 'TEmean')&(sobol_data_for_nonfctn$output_var %in% calcite_outputs)),]
bootstrapped_TEmean$boot_num   = NULL
bootstrapped_TEmean$sobol_type = NULL
bootstrapped_seed              = seed_data_for_nonfctn[(seed_data_for_nonfctn$output_var %in% calcite_outputs),]
bootstrapped_TEmean$epsilon    = bootstrapped_seed$TEseed
graph_data                     = melt(bootstrapped_TEmean, id = "output_var")

g4 = ggplot(graph_data, aes(x = variable, y = value, fill = output_var)) + geom_boxplot(lwd=0.25, outlier.size = 0.2, width = .9) + 
  ggtitle(TeX("")) + 
  xlab("Model Parameters") + theme_bw() + 
  ylab("") + 
  theme(axis.text=element_text(size=15), axis.title = element_text(size = 15)) +
  theme(plot.title = element_blank(),
        axis.text.x = element_text(angle = 270, vjust = 0.5, hjust=1),
        panel.spacing = unit(0, "cm"),
        plot.margin = margin(0, 0, 0, 0, "cm"), 
        plot.caption = element_blank()) + labs(fill="Flush Time")+ guides(fill="none")  + scale_fill_manual(name = "max_calcite",values = c("grey"))

# Editing figures together for paper:
g2 = g2 + scale_x_discrete(labels = TeX(c('curr_inflow_rate' = r"(Q)", 
                                          'curr_diffusion_coef' = r"(D)", 
                                          'curr_porosity' = r"($\phi$)",
                                          'curr_gypsum_rate_constant' = r"($k_{gyp.}$)",
                                          'curr_calcite_rate_constant' = r"($k_{cal.}$)",
                                          'curr_gypsum_surface_area' = r"(Gyp. SA)",
                                          'curr_calcite_surface_area' = r"(Cal. SA)",
                                          'dfn_volume' = r"(dfn vol.)")) )
g4 = g4 + scale_x_discrete(labels = TeX(c('curr_inflow_rate' = r"(Q)", 
                                          'curr_diffusion_coef' = r"(D)", 
                                          'curr_porosity' = r"($\phi$)",
                                          'curr_gypsum_rate_constant' = r"($k_{gyp.}$)",
                                          'curr_calcite_rate_constant' = r"($k_{cal.}$)",
                                          'curr_gypsum_surface_area' = r"(Gyp. SA)",
                                          'curr_calcite_surface_area' = r"(Cal. SA)",
                                          'dfn_volume' = r"(dfn vol.)",
                                          'epsilon' = r"($\epsilon$)")) )

FO_calcite_nonfctn = g2
TE_calcite_nonfctn = g4


# # To get the functional plots, source the make_boxplots stuff from the uncorrelated simulation files.
# source("../uncorr_sobol_indices/make_boxplots.R")
# 
# TE_calcite_fctn    = TE_calcite_fctn + ylim(0,0.5) #+
# FO_calcite_fctn    = FO_calcite_fctn + ylim(-0.1,0.2) #+
TE_calcite_nonfctn = TE_calcite_nonfctn + ylim(0,0.5) #+
FO_calcite_nonfctn = FO_calcite_nonfctn + ylim(-0.1,0.2) #+

row1 <- ggplot() + annotate(geom = 'text', x=1, y=1, label="MC Samples of First-\n Order Sobol' Indices", angle = 90, size = 5) + theme_void() 
row2 <- ggplot() + annotate(geom = 'text', x=1, y=1, label="MC Samples of Total-\n Effect Sobol' Indices", angle = 90, size = 5) + theme_void() 

col1 <- ggplot() + annotate(geom = 'text', x=1, y=1, label="") + theme_void() 
col2 <- ggplot() + annotate(geom = 'text', x=1, y=1, label="") + theme_void() 

# layoutplot <- "
# aeeeeeeeeeeffffffffff
# aeeeeeeeeeeffffffffff
# bgggggggggghhhhhhhhhh
# bgggggggggghhhhhhhhhh
# "
# 
# plotlist <- list(a = row1, b = row2, e= FO_calcite_nonfctn, f=FO_calcite_fctn, g=TE_calcite_nonfctn, h=TE_calcite_fctn)
# 
# wrap_plots(plotlist, guides = 'collect', design = layoutplot)


row1 <- ggplot() + annotate(geom = 'text', x=1, y=1, label="MC Samples of First-\n Order Sobol' Indices", angle = 90, size = 5) + theme_void() 
row2 <- ggplot() + annotate(geom = 'text', x=1, y=1, label="MC Samples of Total-\n Effect Sobol' Indices", angle = 90, size = 5) + theme_void() 

col1 <- ggplot() + annotate(geom = 'text', x=1, y=1, label="") + theme_void() 
col2 <- ggplot() + annotate(geom = 'text', x=1, y=1, label="") + theme_void() 

# layoutplot <- "
# aeeeeeeeeeeffffffffff
# aeeeeeeeeeeffffffffff
# bgggggggggghhhhhhhhhh
# bgggggggggghhhhhhhhhh
# "

layoutplot <- "
aeeeeeeeeeeeeeeeeeeeeee
aeeeeeeeeeeeeeeeeeeeeee
bgggggggggggggggggggggg
bgggggggggggggggggggggg
"

plotlist <- list(a = row1, b = row2, e= FO_calcite_nonfctn, g=TE_calcite_nonfctn)

wrap_plots(plotlist, guides = 'collect', design = layoutplot)






