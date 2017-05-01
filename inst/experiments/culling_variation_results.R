#Necessary functions
library(metaflu)
library(ggplot2)
library(dplyr)
library(doMC)
library(tidyr)
library(purrr)
library(gridExtra)
library(abind)
registerDoMC(cores = 20)
#devtools::install_githhub("renozao/doRNG", force = TRUE)
library(doRNG)
set.seed(123)

#Set up parameters
farm_size <- 50
farm_number <- 200

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
  cull_time = 2,   #time to detect, which will be changed in this simulation
  network_type = "smallworld",
  network_parms = list(dim = 1, size = farm_number, nei = 2.33, p = 0.0596, multiple = FALSE, loops = FALSE),
  stochastic_network = TRUE
)

#vary culling time: 0 to 3 weeks, by day
cull_time_vector <- c(1:21)
sims <- 1000

lapply(cull_time_vector, function(a){
  parms[["cull_time"]] <- a
  results_list <- mclapply(seq_len(sims), function (x){
    patches <- grow_patches_clustered(basic_patches(50,200))
    i_patches <- seed_initial_infection(patches)
    array <- mf_sim(init = i_patches, parameters = parms, times=1:365, n_sims = 1)
    return(array)
  }, mc.cores = 20)
  #results <- do.call("abind",results_list)
  saveRDS(results_list, paste0("cull_time",a,".rds"))
  print(paste("Cull time", a, "saved!"))
  })

#below: the code used to create the culling graphs in presentation_file
library(abind)
cull_time_vector <- c(1:21)
results_list <- lapply (cull_time_vector, function(a){
  results <- readRDS(paste0("inst/experiments/cull_time",a,".rds"))
  results <- do.call("abind", results)
  results_df <- create_pres_df(results, a)
  print(paste0("done with run", a))
  return(results_df)
})

final_df <- bind_rows(results_list)

ggplot(data = final_df) +
  geom_point(aes(x = scenario, y = mean_prop_loss)) +
  labs(title = "Proportion of Loss by Time to Culling", x = "Cull Time", y = "Days") +
  theme_minimal()

ggplot(data = final_df) +
  geom_point(aes(x = scenario, y = mean_proportion_farms)) +
  labs(title = "Proportion of Infected Farms by Time to Culling", x = "Cull Time", y = "Days") +
  theme_minimal()

ggplot(data = final_df) +
  geom_point(aes(x = scenario, y = mean_duration)) +
  labs(title = "Duration of Epidemic by Time to Culling", x = "Cull Time", y = "Days") +
  theme_minimal()

ggplot(data = final_df) +
  geom_point(aes(x = scenario, y = mean_fraction_exposure)) +
  labs(title = "Fraction Exposed by Time to Culling", x = "Cull Time", y = "Days") +
  theme_minimal()

