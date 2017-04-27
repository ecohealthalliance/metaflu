#' @export
basic_patches <- function(farm_size, farm_number){
  initial_cond <- cbind(rpois(farm_number, farm_size), matrix(0, ncol = 4, nrow = farm_number))
  return(initial_cond)
}

#' @export
create_initial_condition <- function(patches){
  infected_patches <- sample(seq_len(nrow(patches)), 1)
  patches[infected_patches, 2] <- 1
  patches[infected_patches, 1] <- patches[infected_patches, 1] - 1
  return(patches)
}

#' @export
seed_initial_infection <- function(patches, multiplier = 1, risk_group = NULL){
  temp_patches <- patches
  if (!is.null(risk_group)){
    temp_patches[risk_group,1] <- temp_patches[risk_group,1]*multiplier
  }
  s <- cumsum(temp_patches[,1])
  random <- sample.int(max(s), 1)
  index <- which.max(s > random)
  patches[index, 2] <- 1
  patches[index, 1] <- patches[index, 1] - 1
  return(patches)
}

#' @export
grow_patches_clustered <- function(old_condition){
  s_num <- round(nrow(old_condition)*1/9)
  old_condition[1:s_num,1] <- old_condition[1:s_num,1]*10
  return(old_condition)
}

#' @export
grow_patches_random <- function(old_condition){
  s_num <- round(nrow(old_condition)*1/9)
  random_numbers <- sample.int(nrow(old_condition), s_num)
  old_condition[random_numbers,1] <- old_condition[random_numbers,1]*10
  return(old_condition)
}
