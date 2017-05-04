
#' Returns a dataframe of simulation number and epidemic duration, the number of days with at least one infection present.
#' @export
#' @param results 4-dimensional array of results (compartment, patch, time, simulation)
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

#' Returns a dataframe of simulation number and lower/median/upper number of infected chickens
#' @export
#' @param results 4-dimensional array of results (compartment, patch, time, simulation)
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

#' Returns a dataframe of simulation number and lower/median/upper number of susceptible chickens
#' @export
#' @param results 4-dimensional array of results (compartment, patch, time, simulation)
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

#' Returns a dataframe of simulation number and lower/median/upper number of recovered chickens
#' @export
#' @param results 4-dimensional array of results (compartment, patch, time, simulation)
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

#' Returns a dataframe of simulation number and whether the outbreak failed (defined as at most 1 infection)
#' @export
#' @param results 4-dimensional array of results (compartment, patch, time, simulation)
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

#' Returns array
#' @export
#' @param results 4-dimensional array of results (compartment, patch, time, simulation)
reduce_epi_array <- function(results){
  non_failure_tf <- !get_failure_array(results)[,2]
  non_failures <- results[,,,non_failure_tf]
  return(non_failures)
}


#' Returns the proportion of epidemic failures, define as the number of failures/total results
#' @export
#' @param failure_results
proportion_failed <- function(failure_results){
  sum(failure_results$failed)/length(failure_results$failed)
}

#' Returns a dataframe of the simulation number and the proportion of chickens lost to disease or culling
#' @export
#' @param results 4-dimensional array of results (compartment, patch, time, simulation)
get_proportion_loss <- function(results){
  sims <- seq_len(dim(results)[4])
  loss <- sapply(sims, function(x){
    s <- results["S",,,x]
    s_by_t <- colSums(s)
    total_i <- max(s_by_t) - min(s_by_t)
    prop_loss <- total_i/max(s_by_t)
    return(prop_loss)
  })
  df <- data.frame(sim = sims, total_i = loss)
  return(df)
}

#' Returns a dataframe of the simulation number and the number of birds infected
#' @export
#' @param results 4-dimensional array of results (compartment, patch, time, simulation)
get_exposure <- function(results){
  infections <- apply(results,4,function(x){
    i <- x["I",,]
    return(sum(colSums(i)))
  })
  df <- data.frame(sim = seq_along(infections), inf_exp = infections)
  return(df)
}

#' Returns a dataframe of the simulation number and the fraction of birds infected
#' @export
#' @param results 4-dimensional array of results (compartment, patch, time, simulation)
get_exposure_fraction <- function(results){
  infections <- apply(results,4,function(x){
    i <- x["I",,]
    initial_s <- sum(x["S",,1])
    times <- dim(x)[3]
    return(sum(colSums(i))/(initial_s*times))
  })
  df <- data.frame(sim = seq_along(infections), inf_exp = infections)
  return(df)
}

#' Returns a dataframe of simulation number, time, and number of birds in the given compartment
#' @export
#' @param compartment character representing which compartment (S, I, R, V, or C) to analyze
#' @param results 4-dimensional array of results (compartment, patch, time, simulation)
get_all_sims <- function(compartment, results){
  c.results <- lapply(seq_len(dim(results)[4]), function(x) results[compartment,,,x])
  c.reduced <- t(sapply(c.results, function(x) colSums(x)))
  df <- data.frame(c.reduced)
  colnames(df) <- seq_len(dim(df)[2])
  df$sim <- seq_len(dim(df)[1])
  c.df <- gather(df, time, pop, -sim, convert= TRUE)
  return(c.df)
}

#' Returns a dataframe of the simulation number and the number of farms affected by the outbreak
#' @export
#' @param results 4-dimensional array of results (compartment, patch, time, simulation)
get_number_farms <- function(results){
  farms <- sapply(seq_len(dim(results)[4]), function(x){
    tot_i <- rowSums(results["I",,,x])
    return(sum(tot_i > 0))
  })
  f.df <- data.frame(sim = seq_len(dim(results)[4]), num_farms = farms)
  return(f.df)
}

#' Returns a vector of the number of culling events in each simulation
#' @export
#' @param results 4-dimensional array of results (compartment, patch, time, simulation)
get_number_culls <- function(results){
  farms <- sapply(seq_len(dim(results)[4]), function(x){
    tot_culls <- rowSums(results["C",,,x])
    return(sum(tot_culls > 0))
  })
}

#' Produces graphs showing how bird loss, farms affected, outbreak duration, and number of birds exposed changes
#' as the given parameter varies along the given range
#' @export
#' @importFrom parallel detectCores
#' @param param_value the name of the parameter to vary
#' @param param_vector the range over which the parameter will vary
#' @param sims the number of simulations to run for each value of the parameter
#' @param farm_num the number of farms in the network
#' @param farm_size the typical farm size
#' @param parms the list of parameters
vary_params <- function(param_name, param_values, sims, num_of_farms, num_of_chickens, parms){
  results_list <- lapply(param_values, function(x){
    parms[[param_name]] <- x
    g_list <- mclapply(seq_len(sims), function(y){
      patches <- grow_patches_clustered(basic_patches(num_of_chickens, num_of_farms))
      i_patches <- seed_initial_infection(patches)
      return(mf_sim(init = i_patches, parameters = parms, times=1:365, n_sims = 1))
    }, mc.cores = detectCores()/2)
    return(do.call("abind", g_list))
  })

  create_pres_df <- function(res_array, scenario){ #results, cull_time (string)
    mean_prop_loss <- mean(get_proportion_loss(res_array)[,2])
    mean_proportion_farms <- mean(get_number_farms(res_array)[,2]/num_of_farms) #denominator: # of farms
    mean_duration <- mean(get_duration_array(res_array)[,2])
    mean_fraction_exposure <- mean(get_exposure_fraction(res_array)[,2])
    prop_failure <- proportion_failed(get_failure_array(res_array))
    df <- data.frame(scenario, mean_prop_loss, mean_proportion_farms, mean_duration,
                     mean_fraction_exposure, prop_failure, stringsAsFactors = FALSE)
    return(df)
  }

  summarized_runs <- lapply(seq_along(results_list), function(x) create_pres_df(results_list[[x]], param_values[x]))

  final_df <- bind_rows(summarized_runs)

  loss <- ggplot(data = final_df) +
    geom_point(aes(x = param_values, y = mean_prop_loss)) +
    labs(title = "Proportion of Loss", x = "Culling Time (Days)", y = "Proportion of Chickens Lost") +
    theme_minimal() +
    scale_x_log10() + scale_y_log10()

  farms <- ggplot(data = final_df) +
    geom_point(aes(x = param_values, y = mean_proportion_farms)) +
    labs(title = "Proportion of Infected Farms", x = "Culling Time (Days)", y = "Proportion of Farms Infected") +
    theme_minimal() +
    scale_x_log10() + scale_y_log10()

  duration <- ggplot(data = final_df) +
    geom_point(aes(x = param_values, y = mean_duration)) +
    labs(title = "Duration of Epidemic ", x = "Culling Time (Days)", y = "Days") +
    theme_minimal() +
    scale_x_log10() + scale_y_log10()

  exposure <- ggplot(data = final_df) +
    geom_point(aes(x = param_values, y = mean_fraction_exposure)) +
    labs(title = "Exposure Index", x = "Culling Time (Days)", y = "Exposure Index") +
    theme_minimal() +
    scale_x_log10() + scale_y_log10()

  lay <- rbind(c(1,2),
               c(3,4))

  grid.arrange(loss, farms, duration, exposure, layout_matrix = lay, top=textGrob("Cull Time Parameter Analysis", gp = gpar(fontsize = 16)))

}


