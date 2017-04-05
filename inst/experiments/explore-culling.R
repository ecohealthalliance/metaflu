# SETUP
knitr::opts_chunk$set(echo = FALSE, fig.width=10, fig.height=6, warning = FALSE)
devtools::load_all()
library(metaflu)
library(ggplot2)
library(dplyr)
library(doMC)
library(tidyr)
library(purrr)
library(gridExtra)
library(abind)
set.seed(123)
registerDoMC(cores = 20)

#FUNCTIONS
basic_patches <- function(farm_size, farm_number){
  initial_cond <- cbind(rpois(farm_number, farm_size), matrix(0, ncol = 3, nrow = farm_number))
  return(initial_cond)
}

create_initial_condition <- function(patches){
  infected_patches <- sample(seq_len(nrow(patches)), 1)
  patches[infected_patches, 2] <- 1
  patches[infected_patches, 1] <- patches[infected_patches, 1] - 1
  return(patches)
}



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
  tau_crit = 0,   #critical suveillance time
  I_crit = 0,     #threshold for reporting
  pi_report = 0, #reporting probability
  pi_detect = 0, #detection probability
  network_type = "smallworld",
  network_parms = list(dim = 1, size = farm_number, nei = 2.33, p = 0.0596, multiple = FALSE, loops = FALSE),
  stochastic_network = TRUE
)

sims <- 1000

get_medians <- function(val){
parms["pi_report"] <- val
sim_list <- mclapply(seq_len(sims), function(num){
  conditions <- create_initial_condition(basic_patches(farm_size, farm_number))
  return(mf_sim(init = conditions, parameters = parms, times=1:365, n_sims = 1))
}, mc.cores = 20 )
basic_results <- do.call("abind", sim_list)

# RETURNS
durations <- get_duration_array(basic_results)
med_durations <- median(durations$duration)

farm_no <- get_number_farms(basic_results)
med_farm <- median(farm_no)

tot_i <- get_tot_infections_array(basic_results)
med_i <- median(tot_i$total_i)

df <- data.frame(value = val, duration = med_durations, farms = med_farm, infections = med_i)

return(df)
}



param_values <- seq(from = 0.1, to = 1.0, by = 0.05)

list_of_medians <- lapply(param_values, function(x) get_medians(x))

median_df <- bind_rows(list_of_medians)

