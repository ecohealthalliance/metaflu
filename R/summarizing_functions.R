
#function to get duration of epidemic
get_duration_array <- function(results){
  sims <- seq_len(dim(results)[4])
  durations <- lapply(sims, function(x) {
  infections <- results["I", , , x]    #gets matrix of patches by times
  i_by_time <- colSums(infections)       #gets total infections by time
  zero_i_col <- min(which(i_by_time == 0)) #gets first column there are zero infections
  duration <- as.numeric(attr(i_by_time, "names")[zero_i_col])
  return(duration - 1)
  })
  return(unlist(durations))
}

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

