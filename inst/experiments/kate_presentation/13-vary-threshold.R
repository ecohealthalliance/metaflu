library(metaflu)
library(ggplot2)
library(dplyr)
library(doMC)
library(tidyr)
library(purrr)
library(gridExtra)
library(abind)
registerDoMC(cores = 20)
library(doRNG)
library(grid)
set.seed(123)
source(rprojroot::find_rstudio_root_file("R","summarizing_functions.R"))

# Scripts to produce results of varying 1 or 2 parameters for final presentation.
# Run either 1D or 2D case (2D currently commented out); change title of saved results

# Setup (1D)
# For the presentation, cull_time varied from 1 to 20, pi_report from 0.05 to 1, and I_crit from 1 to 10.
threshold = c(1:20)
varied_param = "cull_time"

results_list <- vary_params(param_name = varied_param, param_values = threshold, num_of_chickens = 1000)

saveRDS(results_list, "size1000.rds")

print("I AM DONE!")

#####################################################
# Setup (2D)
# For the presentation, farm size varied from (50, 200, 100) and omega from (0.01, 0.03, 0.05)
# num_of_farms = 50 #200 is default
# threshold_1 = c(1:20)
# varied_param_1 = "cull_time"
#
# threshold_2 = 0.01
# varied_param_2 = "omega"
#
# results_list <- vary_params_2d(param1_name = varied_param_1, param1_values = threshold_1, param2_name = varied_param_2,
#                                param2_values = threshold_2)
#
# saveRDS(results_list, "small_omega_small_size.rds")
#
# print("I AM DONE!")
