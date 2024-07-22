# File to run 100 instances of the simulation study in parallel.
library(parallel)

# This is the number of output variables.
num_of_sims = 8*300
num_cores   = 62

sim_study = function(x){
  system(paste("/n/local_linux/R/4.3.0/bin/Rscript sobol_on_seq_design_runs_nonfctn.R", x))
}

par_vector = mclapply(1:num_of_sims, sim_study, mc.cores = num_cores)
