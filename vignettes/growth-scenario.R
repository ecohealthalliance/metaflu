#Poisson Networks

#2nx = (1-a)nx + 10anx
#n = 1/9
#So about 11% need to increase by 10

library(metaflu)
library(ggplot2)
library(dplyr)
library(doMC)
library(tidyr)

set.seed(17)



basic <- function(farm_size, farm_number){
  initial_cond <- basic_nodes(farm_size, farm_number)
  infected_patches <- sample(seq_len(nrow(initial_cond)), 1)
  initial_cond[infected_patches, 2] <- 1
  initial_cond[infected_patches, 1] <- initial_cond[infected_patches, 1] - 1

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
    tau_crit = 0,   #critical suveillance time
    I_crit = 0,     #threshold for reporting
    pi_report = 0, #reporting probability
    pi_detect = 0, #detection probability
    chi = make_net(network_type = "smallworld",
                   network_parms = list(dim = 1, size = farm_number, nei = 2.33, p = 0.0596, multiple = FALSE, loops = FALSE)),
    stochastic_network = TRUE
  )
  x <- mf_sim(init = initial_cond, parameters = parms, times=0:100, n_sims = 20, return_array=TRUE)
  return(x)
}

#registerDoMC(cores=20)

basic_results <- basic(15, 100)

saveRDS(basic_results, "basic_results.rds")

grown <- function(farm_size, farm_number){
  set.seed(123)
  initial_cond <- growth_nodes(farm_size, farm_number)
  infected_patches <- sample(seq_len(nrow(initial_cond)), 1)
  initial_cond[infected_patches, 2] <- 1
  initial_cond[infected_patches, 1] <- initial_cond[infected_patches, 1] - 1

  parms = list(
    beta = 0.0081,   #contact rate for direct transmission
    gamma = 0.167,  #recovery rate
    mu = 0,         #base mortality rate
    alpha = 0.4,      #disease mortality rate
    phi = 0,  #infectiousness of environmental virions
    eta = 0,     #degradation rate of environmental virions
    nu =  0.00,    #uptake rate of environmental virion
    sigma = 0,      #virion shedding rate
    omega = 0.03,   #movement rate
    rho = 0,        #contact  nonlinearity 0=dens-dependent, 1=freq-dependent
    lambda = 0,     #force of infection from external sources
    tau_crit = 0,   #critical suveillance time
    I_crit = 0,     #threshold for reporting
    pi_report = 0, #reporting probability
    pi_detect = 0, #detection probability
    chi = make_net(network_type = "smallworld",
                   network_parms = list(dim = 1, size = farm_number, nei = 2.33, p = 0.0596, multiple = FALSE, loops = FALSE)),
    stochastic_network = TRUE
  )
  x <- mf_sim(init = initial_cond, parameters = parms, times=0:100, n_sims = 100, return_array=FALSE)
  return(x)
}

grown_results2 <- grown(15,100)
saveRDS(grown_results, "grown_results.rds")
