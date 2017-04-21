
library(metaflu)
library(ggplot2)
library(dplyr)
library(doMC)
library(tidyr)
library(purrr)
library(gridExtra)
library(abind)
registerDoMC(cores = 20)
#devtools::install_github("renozao/doRNG", force = TRUE)
library(doRNG)
set.seed(123)





#Set number of farms and ~ number of chickens per farm
farm_number <- 100
farm_size <- 40
print(farm_number)
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
    cull_time = 1,   #time to detect, which will be changed in this simulation
    network_type = "smallworld",
    network_parms = list(dim = 1, size = farm_number, nei = 2.33, p = 0.0596, multiple = FALSE, loops = FALSE),
    stochastic_network = TRUE
    )


cull_time_vector <- c(7, 365, 1095, 1825, 2555, 3285, 3650) #1 week; 1, 3, 5, 7, 9, 10 years

results_list <- lapply(cull_time_vector, function(x){
  print(x)
  parms["cull_time"] <- x
  sims <- 1000
  g_list <- mclapply(seq_len(sims), function(y){
    patches <- grow_patches_clustered(basic_patches(40,100))
    i_patches <- seed_initial_infection(patches)
    return(mf_sim(init = i_patches, parameters = parms, times=1:365, n_sims = 1))
  }, mc.cores = 20 )

  return(do.call("abind", g_list))
})

saveRDS(results_list, "cull_time.rds")

print("graphing infections")
tot_i_list <- lapply(results_list, function(x) get_tot_infections_array(x))
inf_means <- sapply(tot_i_list, function(x) mean(x$total_i))
infection_df <- data.frame(cull_time = cull_time_vector, mean_infections = inf_means)
ggplot(data = infection_df) + geom_point(aes(x = cull_time, y = mean_infections)) +
  labs(x = "days to culling", y = "mean number of infections")

print("graphing duration")
duration_list <- lapply(results_list, function(x) get_duration_array(x))
dmeans <- sapply(duration_list, function(x) mean(x$duration))
duration_df <- data.frame(cull_time = cull_time_vector, mean_durations = dmeans)
ggplot(data = duration_df) + geom_point(aes(x = cull_time, y = mean_durations)) +
  labs(x = "days to culling", y = "duration of epidemic")




