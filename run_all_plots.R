# Creates all plots.
## Author: Alexander C. Murph
## Date: May 2024

###########
## Max Values plots
# Figure 5a
rm(list = ls())
dpi_set = 700
width_set = 14
height_set = 7
source("~/GitLab/sa_for_chemistry/fixed_dfn/max_values/sobol_indices_fixedDFN/make_boxplots.R")
wrap_plots(plotlist, guides = 'collect', design = layoutplot)
ggsave(filename = "fixed_dfn_max_values_nonfctn.png", 
       path = "~/GitLab/sa_for_chemistry/paper_stuff/max_values", 
       width = width_set, height = height_set, device='png', dpi=dpi_set)

# Max Calcite Figure in the Appendix
rm(list = ls())
dpi_set = 700
width_set = 14
height_set = 7
source("~/GitLab/sa_for_chemistry/fixed_dfn/max_values/uncorr_sobol_indices_fixedDFN/make_boxplots.R")
wrap_plots(plotlist, guides = 'collect', design = layoutplot)
ggsave(filename = "fixed_dfn_max_values_fctn.png", 
       path = "~/GitLab/sa_for_chemistry/paper_stuff/max_values", 
       width = width_set, height = height_set, device='png', dpi=dpi_set)

# Figure 5b
rm(list = ls())
dpi_set = 700
width_set = 14
height_set = 7
source("~/GitLab/sa_for_chemistry/variable_dfn/max_values/sobol_indices/make_boxplots.R")
wrap_plots(plotlist, guides = 'collect', design = layoutplot)
ggsave(filename = "variable_dfn_max_values_nonfctn.png", 
       path = "~/GitLab/sa_for_chemistry/paper_stuff/max_values", 
       width = width_set, height = height_set, device='png', dpi=dpi_set)

# Max Calcite Figure in the Appendix
rm(list = ls())
dpi_set = 700
width_set = 14
height_set = 7
source("~/GitLab/sa_for_chemistry/variable_dfn/max_values/uncorr_sobol_indices/make_boxplots.R")
wrap_plots(plotlist, guides = 'collect', design = layoutplot)
ggsave(filename = "variable_dfn_max_values_fctn.png", 
       path = "~/GitLab/sa_for_chemistry/paper_stuff/max_values", 
       width = width_set, height = height_set, device='png', dpi=dpi_set)


########
## Fixed DFN experiment reparam time.
# Figure 9a
rm(list = ls())
dpi_set = 700
width_set = 14
height_set = 7
source("~/GitLab/sa_for_chemistry/fixed_dfn/sobol_indices_fixedDFN/make_boxplots.R")
wrap_plots(plotlist, guides = 'collect', design = layoutplot)
ggsave(filename = "fixedDFN_reparameterizedTime_nonfctn.png", 
       path = "~/GitLab/sa_for_chemistry/paper_stuff", 
       width = width_set, height = height_set, device='png', dpi=dpi_set)

# Fixed DFN figure in the Appendix
rm(list = ls())
dpi_set = 700
width_set = 14
height_set = 7
source("~/GitLab/sa_for_chemistry/fixed_dfn/uncorr_sobol_indices_fixedDFN/make_boxplots.R")
wrap_plots(plotlist, guides = 'collect', design = layoutplot)
ggsave(filename = "fixedDFN_reparameterizedTime_fctn.png", 
       path = "~/GitLab/sa_for_chemistry/paper_stuff", 
       width = width_set, height = height_set, device='png', dpi=dpi_set)

########
## All variable DFN studies not max value
# Figure 7
rm(list = ls())
dpi_set = 700
width_set = 14
height_set = 7
source("~/GitLab/sa_for_chemistry/variable_dfn/sobol_indices/make_boxplots.R")
wrap_plots(plotlist, guides = 'collect', design = layoutplot)
ggsave(filename = "variableDFN_reparameterizedTime_nonfctn.png", 
       path = "~/GitLab/sa_for_chemistry/paper_stuff", 
       width = width_set, height = height_set, device='png', dpi=dpi_set)

# Variable DFN figure in the Appendix
rm(list = ls())
dpi_set = 700
width_set = 14
height_set = 7
source("~/GitLab/sa_for_chemistry/variable_dfn/uncorr_sobol_indices/make_boxplots.R")
wrap_plots(plotlist, guides = 'collect', design = layoutplot)
ggsave(filename = "variableDFN_reparameterizedTime_fctn.png", 
       path = "~/GitLab/sa_for_chemistry/paper_stuff", 
       width = width_set, height = height_set, device='png', dpi=dpi_set)

# Figure 9b
rm(list = ls())
dpi_set = 700
width_set = 14
height_set = 7
source("~/GitLab/sa_for_chemistry/variable_dfn_regulartime/sobol_indices/make_boxplots.R")
wrap_plots(plotlist, guides = 'collect', design = layoutplot)
ggsave(filename = "variableDFN_regularTime_nonfctn.png", 
       path = "~/GitLab/sa_for_chemistry/paper_stuff", 
       width = width_set, height = height_set, device='png', dpi=dpi_set)

# Variable DFN figure for uncorrelated variables in the Appendix
rm(list = ls())
dpi_set = 700
width_set = 14
height_set = 7
source("~/GitLab/sa_for_chemistry/variable_dfn_regulartime/uncorr_sobol_indices/make_boxplots.R")
wrap_plots(plotlist, guides = 'collect', design = layoutplot)
ggsave(filename = "variableDFN_regularTime_fctn.png", 
       path = "~/GitLab/sa_for_chemistry/paper_stuff", 
       width = width_set, height = height_set, device='png', dpi=dpi_set)


#####
## Run all the EDA Stuff:
# Eda graphic in the Appendix
rm(list = ls())
dpi_set = 700
width_set = 14
height_set = 7
source("~/GitLab/sa_for_chemistry/variable_dfn/eda/data_analysis.R")
source("~/GitLab/sa_for_chemistry/fixed_dfn/eda/data_analysis.R")
grid.arrange(g1, g2, ncol = 1)
ggsave(filename = "eda_graphic.png", 
       path = "~/GitLab/sa_for_chemistry/paper_stuff", 
       width = width_set, height = height_set, device='png', dpi=dpi_set)

# Eda graphic in the Appendix
rm(list = ls())
dpi_set = 700
width_set = 14
height_set = 7
source("~/GitLab/sa_for_chemistry/variable_dfn_regulartime/eda/data_analysis.R")
g2
ggsave(filename = "eda_graphic.png", 
       path = "~/GitLab/sa_for_chemistry/paper_stuff", 
       width = width_set, height = height_set, device='png', dpi=dpi_set)

