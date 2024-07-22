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
library(patchwork)
setwd("~/GitLab/sa_for_chemistry/variable_dfn_regulartime/uncorr_sobol_indices")

cvi_colours = list(
  jeffrey_colors =  c('#1f77b4',
                      '#ff7f0e',
                      '#2ca02c',
                      '#d62728',
                      '#9467bd',
                      '#8c564b',
                      '#e377c2',
                      '#7f7f7f',
                      '#bcbd22',
                      '#17becf')
)


cvi_palettes = function(name, n, all_palettes = cvi_colours, type = c("discrete", "continuous")) {
  palette = all_palettes[[name]]
  if (missing(n)) {
    n = length(palette)
  }
  type = match.arg(type)
  out = switch(type,
               continuous = grDevices::colorRampPalette(palette)(n),
               discrete = palette[1:n]
  )
  structure(out, name = name, class = "palette")
}

###############################################################################
###############################################################################
###############################################################################

new_names = c( "dfn_seed",                   "curr_inflow_rate",           "curr_diffusion_coef",        
               "curr_porosity",              "curr_gypsum_rate_constant",
               "curr_calcite_rate_constant", "curr_gypsum_surface_area",   "curr_calcite_surface_area",  
               "dfn_p32",                    "dfn_volume",                
               "backbone_volume",            "backbone_p32",               "final_gypsum",               
               "final_calcite",              "Pe",                        
               "Da_1_gypsum",                "Da_2_gypsum",                "Da_1_calcite",               
               "Da_2_calcite",               "tau",                       
               "calcite_flush_0_nondim",     "gypsum_flush_0_nondim",      "calcite_flush_1_nondim",     
               "gypsum_flush_1_nondim",      "calcite_flush_2_nondim",    
               "gypsum_flush_2_nondim",      "calcite_flush_3_nondim",     "gypsum_flush_3_nondim",      
               "calcite_flush_4_nondim",     "gypsum_flush_4_nondim",     
               "calcite_flush_5_nondim",     "gypsum_flush_5_nondim",      "calcite_flush_6_nondim",     
               "gypsum_flush_6_nondim",      "calcite_flush_7_nondim",    
               "gypsum_flush_7_nondim",      "calcite_flush_8_nondim",     "gypsum_flush_8_nondim",      
               "calcite_flush_9_nondim",     "gypsum_flush_9_nondim",     
               "calcite_flush_10_nondim",    "gypsum_flush_10_nondim",     
               "calcite_flush_0",            "gypsum_flush_0",             "calcite_flush_1",           
               "gypsum_flush_1",             "calcite_flush_2",            "gypsum_flush_2",             
               "calcite_flush_3",            "gypsum_flush_3",            
               "calcite_flush_4",            "gypsum_flush_4",             "calcite_flush_5",            
               "gypsum_flush_5",             "calcite_flush_6",           
               "gypsum_flush_6",             "max_calcite",                "min_gypsum" )

# Variable updates after the EDA.
taopoints_to_remove = c("calcite_flush_0_nondim",     "gypsum_flush_0_nondim",      "calcite_flush_1_nondim",     
                        "gypsum_flush_1_nondim",      "calcite_flush_2_nondim",    
                        "gypsum_flush_2_nondim",      "calcite_flush_3_nondim",     "gypsum_flush_3_nondim",      
                        "calcite_flush_4_nondim",     "gypsum_flush_4_nondim",     
                        "calcite_flush_5_nondim",     "gypsum_flush_5_nondim",      "calcite_flush_6_nondim",     
                        "gypsum_flush_6_nondim",      "calcite_flush_7_nondim",    
                        "gypsum_flush_7_nondim",      "calcite_flush_8_nondim",     "gypsum_flush_8_nondim",      
                        "calcite_flush_9_nondim",     "gypsum_flush_9_nondim",     
                        "calcite_flush_10_nondim",    "gypsum_flush_10_nondim",
                        "gypsum_flush_0", 
                        "gypsum_flush_1",   "gypsum_flush_2",             
                        "gypsum_flush_3",            
                        "gypsum_flush_4",  
                        "gypsum_flush_5",     
                        "gypsum_flush_6",             "max_calcite",                "min_gypsum")
new_names           = new_names[which(!(new_names%in%taopoints_to_remove))]

gypsum_outputs  = c("gypsum_flush_0_nondim", "gypsum_flush_1_nondim",
                    "gypsum_flush_2_nondim", "gypsum_flush_3_nondim", "gypsum_flush_4_nondim",     
                    "gypsum_flush_5_nondim", "gypsum_flush_6_nondim", 
                    "gypsum_flush_7_nondim", "gypsum_flush_8_nondim", "gypsum_flush_9_nondim",     
                    "gypsum_flush_10_nondim", "gypsum_flush_11_nondim",   
                    "gypsum_flush_12_nondim")
outputs_abrv    = c("flush_2", "flush_3", 
                    "flush_4", "flush_5", "flush_6")
calcite_outputs = c("calcite_flush_2", 
                    "calcite_flush_3", "calcite_flush_4", "calcite_flush_5", 
                    "calcite_flush_6")

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
bootstrapped_FOmean            = NULL
bootstrapped_FOmean            = sobol_data_for_fctn[which((sobol_data_for_fctn$sobol_type == 'FOmean')&
                                                             (sobol_data_for_fctn$output_var %in% calcite_outputs)),]
bootstrapped_FOmean$boot_num   = NULL
bootstrapped_FOmean$sobol_type = NULL
graph_data                     = bootstrapped_FOmean
graph_data$output_var          = factor(graph_data$output_var, levels = calcite_outputs)
levels(graph_data$output_var)  = outputs_abrv

g2 = ggplot(graph_data, aes(x = input_var, y = value, fill = output_var)) + geom_boxplot(lwd=0.25, outlier.size = 0.2, 
                                                                                         alpha = 0.5, width = .9) + 
  ggtitle(TeX("")) + 
  theme_bw() + 
  xlab("") + 
  ylab("") + 
  theme(axis.text=element_text(size=15), axis.title = element_text(size = 15)) + 
  theme(plot.title = element_blank(), axis.text.x = element_text(angle = 270, vjust = 0.5, hjust=1), 
        plot.caption = element_blank()) + 
  labs(fill="Time [y]")#+ guides(fill="none")  

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Total-Effect Indices on the Mean Response GP, for calcite variables:
bootstrapped_TEmean            = NULL
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
graph_data                     = bootstrapped_TEmean
graph_data$output_var          = factor(graph_data$output_var, levels = calcite_outputs)
levels(graph_data$output_var)  = outputs_abrv
graph_data$input_var           = factor(graph_data$input_var, levels = c('Da_1_gypsum', 'Da_1_calcite',
                                                                         'Da_2_calcite', 'Da_2_gypsum', 
                                                                         'Pe', 'dfn_volume',
                                                                         'tau',
                                                                         'epsilon'))

g4 = ggplot(graph_data, aes(x = factor(input_var, levels = c('Da_1_gypsum', 'Da_1_calcite',
                                                                      'Da_2_calcite', 'Da_2_gypsum', 
                                                                      'dfn_volume', 'Pe',
                                                                      'tau',
                                                                      'epsilon')), 
                                     y = value, fill = output_var)) + geom_boxplot(lwd=0.25, outlier.size = 0.2, 
                                                                                   alpha = 0.5, width = .9) + 
  ggtitle(TeX("")) + 
  theme_bw() + 
  xlab("Dimensionless Parameters") + 
  ylab("") + 
  theme(axis.text=element_text(size=15), axis.title = element_text(size = 15)) +
  theme(plot.title = element_blank(), axis.text.x = element_text(angle = 270, vjust = 0.5, hjust=1), 
        plot.caption = element_blank()) + 
  labs(fill="Timepoint")+ guides(fill="none")


# Editing figures together for paper:
g2 = g2 + scale_x_discrete(labels = TeX(c('Pe' = r"(Pe)", 
                                          'Da_1_gypsum' = r"($Da_{I}$ Gyp.)", 
                                          'Da_2_gypsum' = r"($Da_{II}$ Gyp.)",
                                          'Da_1_calcite' = r"($Da_{I}$ Cal.)",
                                          'Da_2_calcite' = r"($Da_{II}$ Cal.)",
                                          'dfn_volume' = 'dfn vol.',
                                          'tau' = r"($\tau$)")) )
g4 = g4 + scale_x_discrete(labels = TeX(c('Pe' = r"(Pe)", 
                                          'Da_1_gypsum' = r"($Da_{I}$ Gyp.)", 
                                          'Da_2_gypsum' = r"($Da_{II}$ Gyp.)",
                                          'Da_1_calcite' = r"($Da_{I}$ Cal.)",
                                          'Da_2_calcite' = r"($Da_{II}$ Cal.)",
                                          'dfn_volume' = 'dfn vol.',
                                          'tau' = r"($\tau$)",
                                          'epsilon' = r"($\epsilon$)")) )

g2 = g2 + scale_fill_manual(labels = TeX(c('flush_0' = r"($10^{-7}$)",
                                           'flush_1' = r"($10^{-6}$)",
                                           'flush_2' = r"($10^{-5}$)",
                                           'flush_3' = r"($10^{-4}$)",
                                           'flush_4' = r"($10^{-3}$)",
                                           'flush_5' = r"($10^{-2}$)",
                                           'flush_6' = r"($10^{-1}$)")),
                            values = cvi_palettes("jeffrey_colors", type = "discrete"))
g4 = g4 + scale_fill_manual(labels = TeX(c('flush_0' = r"($10^{-7}$)",
                                           'flush_1' = r"($10^{-6}$)",
                                           'flush_2' = r"($10^{-5}$)",
                                           'flush_3' = r"($10^{-4}$)",
                                           'flush_4' = r"($10^{-3}$)",
                                           'flush_5' = r"($10^{-2}$)",
                                           'flush_6' = r"($10^{-1}$)")),
                            values = cvi_palettes("jeffrey_colors", type = "discrete"))

FO_calcite_fctn = g2
TE_calcite_fctn = g4

layoutplot <- "
aeeeeeeeeeeeeeeeeeeeeee
aeeeeeeeeeeeeeeeeeeeeee
bgggggggggggggggggggggg
bgggggggggggggggggggggg
"

row1 <- ggplot() + annotate(geom = 'text', x=1, y=1, label="MC Samples of First-\n Order Sobol' Indices", angle = 90, size = 5) + theme_void() 
row2 <- ggplot() + annotate(geom = 'text', x=1, y=1, label="MC Samples of Total-\n Effect Sobol' Indices", angle = 90, size = 5) + theme_void() 

col1 <- ggplot() + annotate(geom = 'text', x=1, y=1, label="") + theme_void() 
col2 <- ggplot() + annotate(geom = 'text', x=1, y=1, label="") + theme_void() 

plotlist <- list(a = row1, b = row2, e= FO_calcite_fctn, g=TE_calcite_fctn)

wrap_plots(plotlist, guides = 'collect', design = layoutplot)





