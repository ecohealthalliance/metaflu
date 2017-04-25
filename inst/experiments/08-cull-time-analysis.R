
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
    pi_report = 1, #reporting probability
    pi_detect = 1, #detection probability (1 for cleanliness)
    cull_time = 1,   #time to detect, which will be changed in this simulation
    network_type = "smallworld",
    network_parms = list(dim = 1, size = farm_number, nei = 2.33, p = 0.0596, multiple = FALSE, loops = FALSE),
    stochastic_network = TRUE
    )


cull_time_vector <- c(1:20)

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


#saveRDS(results_list, "cull_time.rds", compress = FALSE)
#delete before committing, as this is a lot of space

#print("graphing infections")
tot_i_list <- lapply(results_list, function(x) get_tot_infections_array(x))
inf_means <- sapply(tot_i_list, function(x) mean(x$total_i))
infection_df <- data.frame(cull_time = cull_time_vector, mean_infections = inf_means)
ggplot(data = infection_df) + geom_point(aes(x = cull_time, y = mean_infections)) +
  labs(x = "days to culling", y = "mean number of infections") + scale_y_log10() + scale_x_log10()

#get oubreak proportions
farm_num <- lapply(results_list, function(x) get_number_farms(x))

outbreak_list <- unlist(lapply(farm_num, function(x){
  return(sum(x > 1)/length(x))
}))

outbreak_df <- data.frame(cull_time = cull_time_vector, props = outbreak_list)
ggplot(data = outbreak_df) + geom_point(aes(x = cull_time, y = props)) +
  labs(x = "days to culling", y = "proportion of outbreaks >1 farm") +
  scale_x_log10() + scale_y_log10()

makegraphs <- 0 #currently set so that the next block, which creates summarizing graphs for each paragraph, won't run. Set to 1 to see these graphs.

create_graph_panel <- function(result_array, value){

  gdurations <- get_duration_array(result_array)

  gs.df <- get_all_sims("S", result_array)
  gi.df <- get_all_sims("I", result_array)
  gr.df <- get_all_sims("R", result_array)

  gfarm_no <- get_number_farms(result_array)

  gtot_i <- get_tot_infections_array(result_array)

  gsir.df <- ggplot() +
    geom_line(data = gs.df, aes(x = time, y = pop, group = sim), alpha = 0.05, color = "blue") +
    geom_line(data = gi.df, aes(x = time, y = pop, group = sim), alpha = 0.05, color = "red") +
    geom_line(data = gr.df, aes(x = time, y = pop, group = sim), alpha = 0.05, color = "green") +
    labs(title = paste("Time to Culling:",value, "days"), x = "Time", y = "Population") +
    scale_colour_manual(name = "Compartment", values=c(S = "blue", I = "red", R = "green")) +
    theme_minimal()
#plot(gsir.df)

  gtotal_i <- ggplot() +
    geom_histogram(data = gtot_i, aes(x = total_i), bins = 100) +
    labs(title = "Infections", x = "Infections", y = "Number of Simulations") +
    theme_minimal()

  gfarm_numbers <- ggplot() +
    geom_histogram(aes(x = gfarm_no), bins = 100) +
    labs(title = "Farm Spread", x = "Farms Infected", y = "Number of Simulations") +
    coord_flip() +
    theme_minimal()

  gduration <- ggplot() +
    geom_histogram(data = gdurations, aes(x = duration), bins = 100) +
    theme_minimal() +
    labs(title = "Duration of Epidemic", x = "Duration in Days", y = "Number of Simulations")

  lay <- rbind(c(1,1,1,4),
               c(1,1,1,4),
               c(2,2,3,3))

  grid.arrange(gsir.df,gtotal_i,gduration,gfarm_numbers, layout_matrix = lay)

}

if (makegraphs > 0){
  graph_panel <- purrr::map2(results_list,cull_time_vector, function(x,y) create_graph_panel(x,y))
}
