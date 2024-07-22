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
library(ggpubr)
library(tidyverse)
setwd("~/GitLab/sa_for_chemistry/variable_dfn/max_values/uncorr_sobol_indices")

###############################################################################
###############################################################################
###############################################################################
new_names = c("dfn_seed",                   "curr_inflow_rate",           "curr_diffusion_coef",        "curr_porosity",              "curr_gypsum_rate_constant", 
              "curr_calcite_rate_constant", "curr_gypsum_surface_area",   "curr_calcite_surface_area",  "dfn_p32",                    "dfn_volume",                
              "backbone_volume",            "backbone_p32",               "final_gypsum",               "final_calcite",              "Pe",                        
              "Da_1_gypsum",                "Da_2_gypsum",                "Da_1_calcite",               "Da_2_calcite",               "tau",                       
              "max_calcite",                "min_gypsum")
outputs_abrv    = c("flush_4",     
                    "flush_5", "flush_6", 
                    "flush_7", "flush_8", "flush_9")
calcite_outputs = c("max_calcite")

# ###############################################################################
# ###############################################################################
# ###############################################################################
# Functional Variables

sobol_data_for_fctn = NULL
seed_data_for_fctn  = NULL
for(file_name in list.files('bootstrap_data_fctn')){
  temp_data              = read.csv(paste('bootstrap_data_fctn/', file_name, sep = ''))
  temp_data$X            = NULL
  
  # if((file_name=="calcite_flush_2_nondim.csv")|(file_name=="seed_variable__calcite_flush_2_nondim__.csv")|
  #    (file_name=="seed_variable__calcite_flush_4_nondim__1.csv") | (file_name=="calcite_flush_4_nondim_1.csv")) next
  
  split_string = unlist(strsplit(file_name, '__'))
  if(split_string[1] == 'seed_variable'){
    temp_data$output_var   = split_string[2]
    names(temp_data)       = c('value', 'output_var', 'input_var', 'sobol_type', 'boot_num')
    seed_data_for_fctn     = rbind(seed_data_for_fctn, temp_data)
  } else {
    names(temp_data)       = c('value', 'output_var', 'input_var', 'sobol_type', 'boot_num')
    sobol_data_for_fctn    = rbind(sobol_data_for_fctn, temp_data)
  } 
}

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# First Order Indices on the Mean Response GP, for calcite variables:
bootstrapped_FOmean            = sobol_data_for_fctn[which((sobol_data_for_fctn$sobol_type == 'FOmean')&
                                                             (sobol_data_for_fctn$output_var %in% calcite_outputs)),]
bootstrapped_FOmean$boot_num   = NULL
bootstrapped_FOmean$sobol_type = NULL
graph_data                     = bootstrapped_FOmean
graph_data$output_var          = factor(graph_data$output_var, levels = calcite_outputs)
levels(graph_data$output_var)  = outputs_abrv

g2 = ggplot(graph_data, aes(x = input_var, y = value, fill = output_var)) + geom_boxplot(lwd=0.25, outlier.size = 0.2, width = .9) + 
  ggtitle(TeX("")) + 
  theme_bw() + 
  xlab("") + 
  ylab("") + 
  theme(axis.text=element_text(size=15), axis.title = element_text(size = 15)) + 
  theme(plot.title = element_blank(),
        axis.text.x = element_text(angle = 270, vjust = 0.5, hjust=1),
        panel.spacing = unit(0, "cm"),
        plot.margin = margin(0, 0, 0, 0, "cm"), 
        plot.caption = element_blank()) + labs(fill="Flush Time")+ guides(fill="none")  + scale_fill_manual(name = "max_calcite",values = c("grey"))

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Total-Effect Indices on the Mean Response GP, for calcite variables:
bootstrapped_TEmean            = sobol_data_for_fctn[which((sobol_data_for_fctn$sobol_type == 'TEmean')&
                                                             (sobol_data_for_fctn$output_var %in% calcite_outputs)),]
bootstrapped_TEmean$boot_num   = NULL
bootstrapped_TEmean$sobol_type = NULL
bootstrapped_seed              = seed_data_for_fctn[(seed_data_for_fctn$output_var %in% calcite_outputs),]
bootstrapped_seed$input_var    = "epsilon"
bootstrapped_seed$sobol_type   = NULL
bootstrapped_seed$boot_num     = NULL
bootstrapped_TEmean            = rbind(bootstrapped_TEmean, bootstrapped_seed)

# bootstrapped_TEmean$epsilon    = bootstrapped_seed$TEseed
# graph_data                     = melt(bootstrapped_TEmean, id = "output_var")
graph_data                     = bootstrapped_TEmean
graph_data$output_var          = factor(graph_data$output_var, levels = calcite_outputs)
levels(graph_data$output_var)  = outputs_abrv
graph_data$input_var           = factor(graph_data$input_var, levels = c('Da_1_gypsum', 'Da_1_calcite',
                                                                         'Da_2_calcite', 'Da_2_gypsum', 
                                                                         'dfn_volume', 'Pe',
                                                                         'tau',
                                                                         'epsilon'))


g4 = ggplot(bootstrapped_TEmean, aes(x = factor(input_var, levels = c('Da_1_gypsum', 'Da_1_calcite',
                                                                                 'Da_2_calcite', 'Da_2_gypsum', 
                                                                                 'dfn_volume', 'Pe',
                                                                                 'tau',
                                                                                 'epsilon')), 
                                     y = value, fill = output_var)) + geom_boxplot(lwd=0.25, outlier.size = 0.2, width = .9) + 
  ggtitle(TeX("")) + 
  theme_bw() + 
  xlab("Dimensionless Parameters") + 
  ylab("") + 
  theme(axis.text=element_text(size=15), axis.title = element_text(size = 15)) +
  theme(plot.title = element_blank(),
        axis.text.x = element_text(angle = 270, vjust = 0.5, hjust=1),
        panel.spacing = unit(0, "cm"),
        plot.margin = margin(0, 0, 0, 0, "cm"), 
        plot.caption = element_blank()) + labs(fill="Flush Time")+ guides(fill="none")  + scale_fill_manual(name = "max_calcite",values = c("grey"))

g2 = g2 + scale_x_discrete(labels = TeX(c('Pe' = r"(Pe)", 
                                          'Da_1_gypsum' = r"($Da_{I}$ Gyp.)", 
                                          'Da_2_gypsum' = r"($Da_{II}$ Gyp.)",
                                          'Da_1_calcite' = r"($Da_{I}$ Cal.)",
                                          'Da_2_calcite' = r"($Da_{II}$ Cal.)",
                                          'tau' = r"($\tau$)",
                                          'dfn_volume' = r"(dfn vol.)")) )
g4 = g4 + scale_x_discrete(labels = TeX(c('Pe' = r"(Pe)", 
                                          'Da_1_gypsum' = r"($Da_{I}$ Gyp.)", 
                                          'Da_2_gypsum' = r"($Da_{II}$ Gyp.)",
                                          'Da_1_calcite' = r"($Da_{I}$ Cal.)",
                                          'Da_2_calcite' = r"($Da_{II}$ Cal.)",
                                          'tau' = r"($\tau$)",
                                          'dfn_volume' = r"(dfn vol.)",
                                          'epsilon' = r"($\epsilon$)")) )
# grid.arrange(g2, g4, nrow = 1, top = textGrob("Sobol' Indices for Calcite Precipitation at Different Flush Times for Functional Variables", gp=gpar(fontsize=15,font=8)))

TE_calcite_fctn = g4
FO_calcite_fctn = g2 

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

plotlist <- list(a = row1, b = row2, e= FO_calcite_fctn, g=TE_calcite_fctn)

wrap_plots(plotlist, guides = 'collect', design = layoutplot)






