
#' @export
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

#' @export
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

#' @export
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


#' @export
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

#' @export
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

#' @export
reduce_epi_array <- function(results){
  non_failure_tf <- !get_failure_array(results)[,2]
  non_failures <- results[,,,non_failure_tf]
  return(non_failures)
}



#' @export
proportion_failed <- function(failure_results){
  sum(failure_results$failed)/length(failure_results$failed)
}

#' @export
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

#' @export
get_exposure <- function(results){
  infections <- apply(results,4,function(x){
    i <- x["I",,]
    return(sum(colSums(i)))
  })
  df <- data.frame(sim = seq_along(infections), inf_exp = infections)
  return(df)
}

#' @export
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


#' @export
get_all_sims <- function(compartment, results){
  c.results <- lapply(seq_len(dim(results)[4]), function(x) results[compartment,,,x])
  c.reduced <- t(sapply(c.results, function(x) colSums(x)))
  df <- data.frame(c.reduced)
  colnames(df) <- seq_len(dim(df)[2])
  df$sim <- seq_len(dim(df)[1])
  c.df <- gather(df, time, pop, -sim, convert= TRUE)
  return(c.df)
}

#' @export
get_number_farms <- function(results){
  farms <- sapply(seq_len(dim(results)[4]), function(x){
    tot_i <- rowSums(results["I",,,x])
    return(sum(tot_i > 0))
  })
  f.df <- data.frame(sim = seq_len(dim(results)[4]), num_farms = farms)
  return(f.df)
}

#' @export
get_number_culls <- function(results){
  farms <- sapply(seq_len(dim(results)[4]), function(x){
    tot_culls <- rowSums(results["C",,,x])
    return(sum(tot_culls > 0))
  })
}

#' @export
vary_params <- function(param_value, param_vector, sims, farm_num, farm_size, parms){

  results_list <- lapply(param_vector, function(x){
    parms[[param_value]] <- x

    g_list <- mclapply(seq_len(sims), function(y){
      patches <- grow_patches_clustered(basic_patches(farm_size, farm_num))
      i_patches <- seed_initial_infection(patches)
      return(mf_sim(init = i_patches, parameters = parms, times=1:365, n_sims = 1))
    }, mc.cores = 20)

    return(do.call("abind", g_list))
  })

  final_df <- bind_rows(results_list)

  closs <- ggplot(data = final_df) +
    geom_point(aes(x = scenario, y = mean_prop_loss)) +
    labs(title = "Proportion of Loss", x = "Culling Time (Days)", y = "Proportion of Chickens Lost") +
    theme_minimal()

  cfarms <- ggplot(data = final_df) +
    geom_point(aes(x = scenario, y = mean_proportion_farms)) +
    labs(title = "Proportion of Infected Farms", x = "Culling Time (Days)", y = "Proportion of Farms Infected") +
    theme_minimal()

  cduration <- ggplot(data = final_df) +
    geom_point(aes(x = scenario, y = mean_duration)) +
    labs(title = "Duration of Epidemic ", x = "Culling Time (Days)", y = "Days") +
    theme_minimal()

  cexposure <- ggplot(data = final_df) +
    geom_point(aes(x = scenario, y = mean_fraction_exposure)) +
    labs(title = "Exposure Index", x = "Culling Time (Days)", y = "Exposure Index") +
    theme_minimal()

  lay <- rbind(c(1,2),
               c(3,4))

  grid.arrange(closs, cfarms, cduration, cexposure, layout_matrix = lay, top=textGrob("Cull Time Parameter Analysis", gp = gpar(fontsize = 16)))

  rm(results_list)
}


