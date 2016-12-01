
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
  
  x <- mf_sim(init = initial_cond, parameters = figS9parms, times=0:1000, n_sims = 10)
  x
}

registerDoMC(cores=4)

recreate_data <- function(funcname){
  a <- c(40,50,64,80,100,125,160,200,256,320,400,500,640,800)
  parMat <- cbind(a, rev(a))
  results <- apply(parMat, 1, function(x) do.call(funcname, list(x[1], x[2])))
  results
}

figS9results <- recreate_data("figS9")


proportions <- unlist(lapply(figS9results, function(x) proportion_failed(get_failure(x))))

df1 <- data.frame(id = c(1:7), fails = proportions)

ggplot(data = df1) +
  geom_line(aes(x = id, y = fails)) +
  #scale_x_continuous(name = "Chickens:Farm", labels = c("40:800", "80:400", "160:200","200:160", "320:100","500:64","640:50")) + 
  scale_y_continuous(name = "Proportion of Epidemic Failures") +
  theme_bw()

duration_list <- lapply(figS9results, get_duration)
lower <- sapply(duration_list,function(x) quantile(unlist(x$days_greater), probs = c(0.025)))
upper <- sapply(duration_list, function(x) quantile(unlist(x$days_greater), probs = c(0.975)))
median <- sapply(duration_list, function(x) quantile(unlist(x$days_greater), probs = c(0.5)))

df2 <- data.frame(lower, median, upper)

ggplot(data = df2) +
  geom_line(aes(x = 1:7, y = median)) +
  geom_line(aes(x = 1:7, y = lower), linetype = "dashed") +
  geom_line(aes(x = 1:7, y = upper), linetype = "dashed") + 
  #scale_x_continuous(name = "Chickens:Farm", labels = c("80:400", "160:200","200:160", "320:100","500:64")) + 
  scale_y_continuous(name = "Median Duration (days)") +
  theme_bw()




res <- figS9(100, 3200)

#recreate_data <- function(funcname){
#  a <- c(5,8,10,16,20,25,32,40,50,64,80,100,125,160,200,256,320,400,500,640,800,1000,1280,1600,2000,3200,4000,6400)
#  parMat <- cbind(a, rev(a))
#  results <- apply(parMat, 1, function(x) do.call(funcname, list(x[1], x[2])))
#  results
#}


