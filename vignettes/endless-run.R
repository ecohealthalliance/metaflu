library(metaflu)
library(doMC)


basic_nodes <- function(farm_size, farm_number){
  initial_cond <- cbind(rpois(farm_number, farm_size), matrix(0, ncol = 3, nrow = farm_number))
  return(initial_cond)
}

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
  x <- mf_sim(init = initial_cond, parameters = parms, times=0:100, n_sims = 3, return_array=TRUE)
  return(x)
}

registerDoMC(cores=3)

basic_results <- basic(15, 100)