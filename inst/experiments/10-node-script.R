library(metaflu)
library(ggplot2)
library(dplyr)
library(doMC)
library(tidyr)
library(purrr)
library(gridExtra)
library(abind)
#registerDoMC(cores = 20)
#devtools::install_github("renozao/doRNG", force = TRUE)
library(doRNG)
library(feather)
P <- rprojroot::find_package_root_file

#Set number of farms and ~ number of chickens per farm

farm_number <- 100
farm_size <- 40

parms = list(
  beta = 1.44456,   #contact rate for direct transmission
  gamma = 0.167,  #recovery rate
  mu = 0,         #base mortality rate
  alpha = 0.4,      #disease mortality rate
  phi = 0,  #infectiousness of environmental virions
  eta = 0,     #degradation rate of environmental virions
  nu =  0.00,    #uptake rate of environmental virion
  sigma = 0,      #virion shedding rate
  omega = 0.03,   #movement rate
  rho = 0.85256,        #contact  nonlinearity 0=dens-dependent, 1=freq-dependent
  lambda = 0,     #force of infection from external sources
  tau_crit = 5,   #critical suveillance time
  I_crit = 1,     #threshold for reporting
  pi_report = 0.9, #reporting probability
  pi_detect = 0.9, #detection probability
  cull_time = 5,   #time to detect
  network_type = "smallworld",
  network_parms = list(dim = 1, nei = 2.33, p = 0.0596, multiple = FALSE, loops = FALSE),
  stochastic_network = TRUE
)


farm_size_vec <- c(50, 100, 200, 500, 1000, 1500)
farm_number_vec <- c(20, 50, 100, 200, 500, 750)

sims <- 1000

set.seed(123)

size_list <- lapply(farm_size_vec, function(x){
  farm_size <- x
  number_list <- lapply(farm_number_vec, function(y){
    farm_number <- y
    print(paste("Farm Size:",x,"    Farm Number:",y))
    parms$network_parms$size <- farm_number
    r_list <- mclapply(seq_len(sims), function(z){
      patches <- basic_patches(farm_size, farm_number)
      i_patches <- seed_initial_infection(patches)
      return(mf_sim(init = i_patches, parameters = parms, times = 1:365, n_sims = 1))
    }, mc.cores = 20)
    saveRDS(r_list, paste0(x,"chickens",y,"farms.rds"))
    print("Saved.")
    return(y)
  })
  return(x)
})
