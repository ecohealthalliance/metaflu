library(metaflu)
library(ggplot2)
library(dplyr)
library(doMC)
library(tidyr)
library(purrr)
library(gridExtra)
library(abind)
#devtools::install_github("renozao/doRNG", force = TRUE)
library(doRNG)
set.seed(123)

#Function to apply multiple values of a parameter
#inputs: parameter name, vector of parameter values
vary_params <- function(param_value, param_vector){

  #create list of results for varying the given parameter
  results_list <- lapply(param_vector, function(x){
    parms[[param_value]] <- x
    sims <- 1000
    g_list <- lapply(seq_len(sims), function(y){
      patches <- grow_patches_clustered(basic_patches(40,100))
      i_patches <- seed_initial_infection(patches)
      return(mf_sim(init = i_patches, parameters = parms, times=1:365, n_sims = 1))
    })

    return(do.call("abind", g_list))
  })

  #graph the duration of the outbreak for each param value
  dlist <- lapply(results_list, function(x) get_duration_array(x))
  dmeans <- sapply(dlist, function(x) mean(x$duration))
  duration_df <- data.frame(values = param_vector, mean_durations = dmeans)
  ggplot(data = duration_df) + geom_point(aes(x = values, y = mean_durations)) +
    labs(x = "parameter values", y = "mean duration") + theme_classic()
  #+ scale_y_log10() + scale_x_log10()

  #graph severity of outbreak for each param value
  farm_num <- lapply(results_list, function(x) get_number_farms(x))

  outbreak_list <- unlist(lapply(farm_num, function(x){
    return(sum(x > 1)/length(x))
  }))

  outbreak_df <- data.frame(values = param_vector, props = outbreak_list)
  ggplot(data = outbreak_df) + geom_point(aes(x = values, y = props)) +
    labs(x = "parameter", y = "proportion of outbreaks >1 farm")

  #unclutter environment
  rm(results_list)
}


#Set parameters for network: always 100 farms, with 50 chickens each for the small run
farm_number <- 100

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
  network_parms = list(dim = 1, size = farm_number, nei = 2.33, p = 0.0596, multiple = FALSE, loops = FALSE),
  stochastic_network = TRUE
)

#For small farms
farm_size <- 50

#vary culling time
cull_time_vector <- c(1:5)
cull_time_results <- vary_params("cull_time", cull_time_vector)
print(cull_time) #to check if it resets to 1

break

#vary reporting/detection probability
prob_vector <- seq(0.1,1,0.1)
prob_results <- vary_params("pi_report", prob_vector) #pi-detect is set to 1, so there's only 1 real probability

#vary movement rate
rate_vector <- seq(0.01,0.1, 0.01)
rate_results <- vary_params("omega", rate_vector)


