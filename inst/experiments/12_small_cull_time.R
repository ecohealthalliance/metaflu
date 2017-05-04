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

#Setup
sims = 100
num = 200
size = 50
cull_times = seq(0.1, 1, 0.1)
parms = list(
  beta = 1.44456,   #contact rate for direct transmission
  gamma = 0.167,  #recovery rate
  mu = 0,         #base mortality rate
  alpha = 0.4,      #disease mortality rate
  phi = 0,  #infectiousness of environmental virions
  eta = 0,     #degradation rate of environmental virions
  nu =  0.00,    #uptake rate of environmental virion
  sigma = 0,      #virion shedding rate
  omega = 0.03,   #movement rate (look at varying this too--inversely dependent with farm size)
  rho = 0.85256,        #contact  nonlinearity 0=dens-dependent, 1=freq-dependent
  lambda = 0,     #force of infection from external sources
  tau_crit = 5,   #critical suveillance time
  I_crit = 1,     #threshold for reporting
  pi_report = 1, #reporting probability (vary this too!)
  pi_detect = 1, #detection probability (1 for cleanliness)
  cull_time = 1,   #time to detect, which will be changed in this simulation
  network_type = "smallworld",
  network_parms = list(dim = 1, size = num, nei = 2.33, p = 0.0596, multiple = FALSE, loops = FALSE),
  stochastic_network = TRUE
)


vary_params(param_name = "cull_time", param_values = cull_times, sims = 1000, num_of_farms = 200, num_of_chickens = 40, parms = parms)


