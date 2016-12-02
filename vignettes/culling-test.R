
library(metaflu)
library(ggplot2)
library(dplyr)
library(doMC)
library(tidyr)


#function to get duration of epidemic
get_duration <- function(results){
  duration <- results %>%
    filter(class == "I") %>%
    group_by(sim, time) %>%
    summarize(infectious = sum(population))%>%
    group_by(sim) %>%
    summarize(days_greater = sum(infectious > 0 ))
  duration
}

#function to get cross-sectional infections over time
get_infections <- function(results){
  infections <- results %>%
    filter(class == "I") %>%
    group_by(sim, time) %>%
    summarize(infectious = sum(population)) %>%
    group_by(time) %>%
    summarize(middle = median(infectious), lower = quantile(infectious, probs = c(0.025)), upper = quantile(infectious, probs = c(0.975)))
  infections
}

#function to get cross-sectional susceptible over time
get_susceptibles<- function(results){
  susceptibles <- results %>%
    filter(class == "S") %>%
    group_by(sim, time) %>%
    summarize(susceptibles = sum(population)) %>%
    group_by(time) %>%
    summarize(middle = median(susceptibles), lower = quantile(susceptibles, probs = c(0.025)), upper = quantile(susceptibles, probs = c(0.975)))
  susceptibles
}

# function to get cross-sectional recovered over time
get_recovered<- function(results){
  recovered <- results %>%
    filter(class == "R") %>%
    group_by(sim, time) %>%
    summarize(recovereds = sum(population)) %>%
    group_by(time) %>%
    summarize(middle = median(recovereds), lower = quantile(recovereds, probs = c(0.025)), upper = quantile(recovereds, probs = c(0.975)))
  recovered
}

#function to get epidemic failure rate
get_failure<-function(results){
  check_fails <- function(time, population, patch, class){
    initial<- which(time == 1 & class == "I" & population > 0)
    all<-which(class == "I" & population > 0)
    identical(patch[initial],unique(patch[all]))
  }
  failures <- results %>%
    group_by(sim) %>%
    summarize(failed = check_fails(time, population, patch, class))
  failures
}

proportion_failed <- function(failure_results){
  sum(failure_results$failed)/length(failure_results$failed)
}

#function to get total number of infections
get_tot_infections <- function(results){
  tot_infections <- results %>%
    filter(class == "S") %>%
    group_by(sim, time) %>%
    summarize(tots = sum(population)) %>%
    group_by(sim) %>%
    summarize(total_i = max(tots) - min(tots))
}


set.seed(17)

figS9 <- function(farm_size, farm_number){
  initial_cond <- matrix(c(farm_size, 0, 0, 0), nrow=farm_number, ncol=4, byrow=TRUE)
  infected_patches <- sample(seq_len(nrow(initial_cond)), 2)
  initial_cond[infected_patches, 2] <- 1
  initial_cond[infected_patches, 1] <- initial_cond[infected_patches, 1] - 1

  figS9parms = list(
    beta = 0.004,   #contact rate for direct transmission
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
    tau_crit = 1,   #critical suveillance time
    I_crit = 5,     #threshold for reporting
    pi_report = .9, #reporting probability
    pi_detect = .9, #detection probability
    chi = make_net(network_type = "smallworld",
                   network_parms = list(dim = 1, size = farm_number, nei = 2.33, p = 0.0596, multiple = FALSE, loops = FALSE)),
    stochastic_network = TRUE
  )

  x <- mf_sim(init = initial_cond, parameters = figS9parms, times=0:1000, n_sims = 100)
  x
}

registerDoMC(cores=4)

recreate_data <- function(funcname){
  a <- c(40,50,64,80,100,125,160,200,256,320,400,500,640,800)
  parMat <- cbind(a, rev(a))
  results <- apply(parMat, 1, function(x) get_outputs(funcname,x))
  results
}

get_outputs <- function(funcname, matrow){
  print(paste("Starting",matrow[1],"by",matrow[2],"simulation."))
  results <- do.call(funcname, list(matrow[1], matrow[2]))
  duration <- get_duration(results)
  lower_d <- quantile(unlist(duration$days_greater), probs = c(0.025))
  upper_d <- quantile(unlist(duration$days_greater), probs = c(0.975))
  median_d <- quantile(unlist(duration$days_greater), probs = c(0.5))
  failure <- proportion_failed(get_failure(results))
  abundance <- get_tot_infections(results)
  lower_a <- quantile(unlist(abundance$total_i), probs = c(0.025))
  upper_a <- quantile(unlist(abundance$total_i), probs = c(0.975))
  median_a <- quantile(unlist(abundance$total_i), probs = c(0.5))
  df <- data.frame(failure, lower_d, median_d, upper_d, lower_a, median_a, upper_a)
  df
}


figS9results <- bind_rows(recreate_data("figS9"))

full_recreate_data <- function(funcname){
  a <- c(5,8,10,16,20,25,32,40,50,64,80,100,125,160,200,256,320,400,500,640,800,1000,1280,1600,2000,3200,4000,6400)
  parMat <- cbind(a, rev(a))
  results <- apply(parMat, 1, function(x) get_outputs(funcname, x))
  results
}

fullfigS9results <- bind_rows(full_recreate_data("figS9"))

profvis::profvis({
  x <- figS9(200,160)
})

#longest are separate and arrange

plot(x = seq_along(fullfigS9results$median_d), y = fullfigS9results$median_d, lty = 1)

