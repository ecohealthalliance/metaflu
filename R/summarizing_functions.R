
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
  return(data.frame(sim = sims, duration = unlist(durations)))
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
get_infectious_array <- function(results){
  times <- seq_len(dim(results)[3])
  summaries <- sapply(times, function(x){
  t0 <- results["I", , x, ]
  i_by_sim <- colSums(t0)
  time_summary <- quantile(i_by_sim, probs = c(0.025, 0.5, 0.975), names = FALSE)
  return(time_summary)
  })
  df <- data.frame(time = times-1, lower = summaries[1,], median = summaries[2,], upper = summaries[3,])
  return(df)
}

get_infectious <- function(results){
  infections <- results %>%
    filter(class == "I") %>%
    group_by(sim, time) %>%
    summarize(infectious = sum(population)) %>%
    group_by(time) %>%
    summarize(middle = median(infectious), lower = quantile(infectious, probs = c(0.025)), upper = quantile(infectious, probs = c(0.975)))
  infections
}

#function to get cross-sectional susceptible over time
get_susceptibles_array <- function(results){
  times <- seq_len(dim(results)[3])
  summaries <- sapply(times, function(x){
    t0 <- results["S", , x, ]
    s_by_sim <- colSums(t0)
    time_summary <- quantile(s_by_sim, probs = c(0.025, 0.5, 0.975), names = FALSE)
    return(time_summary)
  })
  df <- data.frame(time = times-1, lower = summaries[1,], median = summaries[2,], upper = summaries[3,])
  return(df)
}

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
get_recovered_array <- function(results){
  times <- seq_len(dim(results)[3])
  summaries <- sapply(times, function(x){
    t0 <- results["R", , x, ]
    r_by_sim <- colSums(t0)
    time_summary <- quantile(r_by_sim, probs = c(0.025, 0.5, 0.975), names = FALSE)
    return(time_summary)
  })
  df <- data.frame(time = times-1, lower = summaries[1,], median = summaries[2,], upper = summaries[3,])
  return(df)
}

get_recovered<- function(results){
  recovered <- results %>%
    filter(class == "R") %>%
    group_by(sim, time) %>%
    summarize(recovereds = sum(population)) %>%
    group_by(time) %>%
    summarize(middle = median(recovereds), lower = quantile(recovereds, probs = c(0.025)), upper = quantile(recovereds, probs = c(0.975)))
  recovered
}

#function to get epidemic failure T/F -- true if infectious do not spread beyond one patch (when seeded)
get_failure_array <- function(results){
  sims <- seq_len(dim(results)[4])
  f_info <- sapply(sims, function(x){
    i <- results["I",,,x]
    initials <- which(i[,1] > 0)
    new_i <- i[-initials,]
    failure <- max(colSums(new_i)) == 0
  })
  df <- data.frame(sim = sims, failed = f_info)
  return(df)
}


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

#function to get total number of infections per simulation
get_tot_infections_array <- function(results){
  sims <- seq_len(dim(results)[4])
  tot_i <- sapply(sims, function(x){
    s <- results["S",,,x]
    s_by_t <- colSums(s)
    total_i <- max(s_by_t) - min(s_by_t)
    return(total_i)
  })
  df <- data.frame(sim = sims, total_i = tot_i)
  return(df)
}

get_tot_infections <- function(results){
  tot_infections <- results %>%
    filter(class == "S") %>%
    group_by(sim, time) %>%
    summarize(tots = sum(population)) %>%
    group_by(sim) %>%
    summarize(total_i = max(tots) - min(tots))
  return(tot_infections)
}


get_all_sims <- function(compartment, results){
  c.results <- lapply(seq_len(dim(results)[4]), function(x) results[compartment,,,x])
  c.reduced <- t(sapply(c.results, function(x) colSums(x)))
  df <- data.frame(c.reduced)
  colnames(df) <- seq_len(dim(df)[2])
  df$sim <- seq_len(dim(df)[1])
  c.df <- gather(df, time, pop, -sim, convert= TRUE)
  return(c.df)
}

get_number_farms <- function(results){
  farms <- sapply(seq_len(dim(results)[4]), function(x){
    tot_i <- rowSums(results["I",,,x])
    return(sum(tot_i > 0))
  })
  return(farms)
}

get_number_culls <- function(results){
  farms <- sapply(seq_len(dim(results)[4]), function(x){
    tot_culls <- rowSums(results["C",,,x])
    return(sum(tot_culls > 0))
  })
}

basic_patches <- function(farm_size, farm_number){
  initial_cond <- cbind(rpois(farm_number, farm_size), matrix(0, ncol = 4, nrow = farm_number))
  return(initial_cond)
}

create_initial_condition <- function(patches){
  infected_patches <- sample(seq_len(nrow(patches)), 1)
  patches[infected_patches, 2] <- 1
  patches[infected_patches, 1] <- patches[infected_patches, 1] - 1
  return(patches)
}

