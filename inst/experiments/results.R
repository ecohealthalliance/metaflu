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
#devtools::install_github("renozao/doRNG", force = TRUE)
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
  cull_time = 1,   #time to detect, which will be changed in this simulation
  network_type = "smallworld",
  network_parms = list(dim = 1, size = farm_number, nei = 2.33, p = 0.0596, multiple = FALSE, loops = FALSE),
  stochastic_network = TRUE
)

#vary culling time: 0 to 3 weeks, by day
cull_time_vector <- c(1)

results_list <- lapply(cull_time_vector, function(x){
  parms[["cull time"]] <- x
  sims <- 1000
  g_list <- mclapply(seq_len(sims), function(y){
    patches <- grow_patches_clustered(basic_patches(40,100))
    i_patches <- seed_initial_infection(patches)
    return(mf_sim(init = i_patches, parameters = parms, times=1:365, n_sims = 1))
  }, mc.cores = 20)

  return(do.call("abind", g_list))
})

#temporary failsafe
#saveRDS(results_list, file == results, compress = FALSE)

temp <- lapply(seq_along(results_list), function(a){
  df.d <- get_duration_array(results_list[[a]])
  df.e <- get_exposure(results_list[[a]])
  df.ef <- get_exposure_fraction(results_list[[a]])
  df.f <- get_number_farms(results_list[[a]])
  df.p <- get_proportion_loss(results_list[[a]])
  df <- full_join(df.d, df.e, by = "sim") %>%
    full_join(df.ef, by = "sim") %>%
    full_join(df.f, by = "sim") %>%
    full_join(df.p, by = "sim")
  df$cull_time <- cull_time_vector[a]
  return(df)
})

full_df <- bind_rows(temp)

#graph duration of outbreak for each value
dlist <- lapply(seq_along(results_list), function(x) get_duration_array(results_list[[x]]))
duration_df <- data.frame(values = cull_time_vector, mean_durations = dmeans)
duration <- ggplot(data = duration_df) +
  geom_point(aes(x = values, y = mean_durations)) +
  labs(x = "cull time", y = "mean duration") +
  theme_classic()
plot(duration)

#graph severity of outbreak for each param value
farm_num <- lapply(results_list, function(x) get_number_farms(x))

outbreak_list <- unlist(lapply(farm_num, function(x){
  return(sum(x > 5)/length(x))
}))

outbreak_df <- data.frame(values = cull_time_vector, props = outbreak_list)
outbreak <- ggplot(data = outbreak_df) + geom_point(aes(x = values, y = props)) +
  labs(x = "cull time", y = "proportion of outbreaks >5 farms")
plot(outbreak)

#unclutter environment
rm(results_list)


