# File to run 100 instances of the simulation study in parallel.
library(parallel)

# This is the number of output variables.
lower_sim_num = 5
upper_sim_num = 6
num_of_sims   = (7*7*300)/5
num_cores     = 62

sim_study = function(x){
  system(paste("/n/local_linux/R/4.3.0/bin/Rscript uncorr_sobol_on_seq_design_runs_fctn.R", x))
}

par_vector = mclapply(( (lower_sim_num-1)*num_of_sims + 1 ):((lower_sim_num)*num_of_sims), sim_study, mc.cores = num_cores)
