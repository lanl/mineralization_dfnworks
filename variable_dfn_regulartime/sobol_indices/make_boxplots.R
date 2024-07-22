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
setwd("~/GitLab/sa_for_chemistry/variable_dfn_regulartime/sobol_indices")

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
####### Graphs for the non-functional variables.
sobol_data_for_nonfctn = NULL
seed_data_for_nonfctn  = NULL
for(file_name in list.files('bootstrap_data_nonfctn')){
  # if( (file_name == "seed_variable__calcite_flush_2_nondim__1.csv")|(file_name == "calcite_flush_2_nondim_1.csv") ) next
  temp_data              = read.csv(paste('bootstrap_data_nonfctn/', file_name, sep = ''))
  temp_data$X            = NULL
  
  split_string = unlist(strsplit(file_name, '__'))
  if(split_string[1] == 'seed_variable'){
    temp_data$output_var   = split_string[2]
    seed_data_for_nonfctn  = rbind(seed_data_for_nonfctn, temp_data)
  } else {
    sobol_data_for_nonfctn = rbind(sobol_data_for_nonfctn, temp_data)
  }
}

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

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# First Order Indices on the Mean Response GP, for calcite variables:
bootstrapped_FOmean            = sobol_data_for_nonfctn[which((sobol_data_for_nonfctn$sobol_type == 'FOmean')&(sobol_data_for_nonfctn$output_var %in% calcite_outputs)),]
bootstrapped_FOmean$boot_num   = NULL
bootstrapped_FOmean$sobol_type = NULL
graph_data                     = melt(bootstrapped_FOmean, id. = "output_var")
graph_data$output_var          = factor(graph_data$output_var, levels = calcite_outputs)
levels(graph_data$output_var)  = outputs_abrv

g2 = ggplot(graph_data, aes(x = variable, y = value, fill = output_var)) + geom_boxplot(lwd=0.25, outlier.size = 0.2, 
                                                                                        alpha = 0.5, width = .9) + 
  ggtitle(TeX("")) + theme_bw() + 
  xlab("") + 
  ylab("") + 
  theme(axis.text=element_text(size=15), axis.title = element_text(size = 15)) + 
  theme(plot.title = element_blank(),
        axis.text.x = element_text(angle = 270, vjust = 0.5, hjust=1),
        panel.spacing = unit(0, "cm"),
        plot.margin = margin(0, 0, 0, 0, "cm"), 
        plot.caption = element_blank()) + labs(fill="Time [y]")  +
  scale_fill_manual(values = cvi_palettes("jeffrey_colors", type = "discrete")) #+
  # guides(fill="none")  


# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Total-Effect Indices on the Mean Response GP, for calcite variables:
bootstrapped_TEmean            = sobol_data_for_nonfctn[which((sobol_data_for_nonfctn$sobol_type == 'TEmean')&(sobol_data_for_nonfctn$output_var %in% calcite_outputs)),]
bootstrapped_TEmean$boot_num   = NULL
bootstrapped_TEmean$sobol_type = NULL
bootstrapped_seed              = seed_data_for_nonfctn[(seed_data_for_nonfctn$output_var %in% calcite_outputs),]
bootstrapped_TEmean$epsilon    = bootstrapped_seed$TEseed
graph_data                     = melt(bootstrapped_TEmean, id = "output_var")
graph_data$output_var          = factor(graph_data$output_var, levels = calcite_outputs)
levels(graph_data$output_var)  = outputs_abrv

g4 = ggplot(graph_data, aes(x = variable, y = value, fill = output_var)) + geom_boxplot(lwd=0.25, 
                                                                                        outlier.size = 0.2, 
                                                                                        alpha = 0.5, width = .9) + 
  ggtitle(TeX("")) + 
  xlab("Model Parameters") + 
  ylab("") + theme_bw()+ 
  theme(axis.text=element_text(size=15), axis.title = element_text(size = 15)) +
  theme(plot.title = element_blank(),
        axis.text.x = element_text(angle = 270, vjust = 0.5, hjust=1),
        panel.spacing = unit(0, "cm"),
        plot.margin = margin(0, 0, 0, 0, "cm"), 
        plot.caption = element_blank()) + labs(fill="Time [y]") #+
  # guides(fill="none")  

# Editing figures together for paper:
g2 = g2 + scale_x_discrete(labels = TeX(c('curr_inflow_rate' = r"(Q)", 
                                          'curr_diffusion_coef' = r"(D)", 
                                          'curr_porosity' = r"($\phi$)",
                                          'curr_gypsum_rate_constant' = r"($k_{gyp.}$)",
                                          'curr_calcite_rate_constant' = r"($k_{cal.}$)",
                                          'curr_gypsum_surface_area' = r"(Gyp. SA)",
                                          'dfn_volume' = r'(dfn vol.)',
                                          'curr_calcite_surface_area' = r"(Cal. SA)")) )
g4 = g4 + scale_x_discrete(labels = TeX(c('curr_inflow_rate' = r"(Q)", 
                                          'curr_diffusion_coef' = r"(D)", 
                                          'curr_porosity' = r"($\phi$)",
                                          'curr_gypsum_rate_constant' = r"($k_{gyp.}$)",
                                          'curr_calcite_rate_constant' = r"($k_{cal.}$)",
                                          'curr_gypsum_surface_area' = r"(Gyp. SA)",
                                          'curr_calcite_surface_area' = r"(Cal. SA)",
                                          'dfn_volume' = r'(dfn vol.)',
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

FO_calcite_nonfctn = g2
TE_calcite_nonfctn = g4


# source("../uncorr_sobol_indices/make_boxplots.R")
# 
# TE_calcite_fctn    = TE_calcite_fctn + ylim(-0.1,0.75) #+
# FO_calcite_fctn    = FO_calcite_fctn + ylim(-0.05,0.4) #+
TE_calcite_nonfctn = TE_calcite_nonfctn + ylim(-0.1,0.75) #+
FO_calcite_nonfctn = FO_calcite_nonfctn + ylim(-0.1,0.4) #+

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

layoutplot <- "
aeeeeeeeeeeeeeeeeeeeeee
aeeeeeeeeeeeeeeeeeeeeee
bgggggggggggggggggggggg
bgggggggggggggggggggggg
"

plotlist <- list(a = row1, b = row2, e= FO_calcite_nonfctn, g=TE_calcite_nonfctn)

wrap_plots(plotlist, guides = 'collect', design = layoutplot)


