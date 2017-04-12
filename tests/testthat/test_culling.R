library(abind)
library(doMC)
library(metaflu)

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
  tau_crit = 2,   #critical suveillance time
  I_crit = 3,     #threshold for reporting
  pi_report = 1, #reporting probability
  pi_detect = 1, #detection probability
  cull_time = 1, #detection cull time
  network_type = "smallworld",
  network_parms = list(dim = 1, size = farm_number, nei = 2.33, p = 0.0596, multiple = FALSE, loops = FALSE),
  stochastic_network = TRUE
)

sims <- 101

sim_list <- mclapply(seq_len(sims), function(num){
  conditions <- create_initial_condition(basic_patches(farm_size, farm_number))
  return(mf_sim(init = conditions, parameters = parms, times=1:365, n_sims = 1))
}, mc.cores = detectCores())
basic_results <- do.call("abind", sim_list)

non_zeros <- apply(basic_results,4, function(x) {
  culled <- rowSums(x["C",,])
  patches <- which(culled > 0)
  return(x[,patches,])
})

cull_mat <- do.call("abind",list(non_zeros, along = 2))

context("Culling")


test_that("Cull status numbers are correct", {
  truths <- apply(cull_mat,2,function(x){
    correct_numbers <- all(x["C",] %in% c(0,1,2))
    mono_inc <- all(x["C",] == cummax(x["C",]))
    return(c(correct_numbers, mono_inc))
  })
  expect_true(all(truths))
})


test_that("Culled farms remain empty", {
  truths <- apply(cull_mat, 2, function(x){
    first2 <- which.max(x["C",] == 2)
    return(all(x[1:4,first2:dim(x)[2]] == 0))
  })
  expect_true(all(truths))
})

move_parms = list(
  beta = 1.44456,   #contact rate for direct transmission
  gamma = 0.167,  #recovery rate
  mu = 0,         #base mortality rate
  alpha = 0.4,      #disease mortality rate
  phi = 0,  #infectiousness of environmental virions
  eta = 0,     #degradation rate of environmental virions
  nu =  0.00,    #uptake rate of environmental virion
  sigma = 0,      #virion shedding rate
  omega = 0.70,   #movement rate
  rho = 0.85256,        #contact  nonlinearity 0=dens-dependent, 1=freq-dependent
  lambda = 0,     #force of infection from external sources
  tau_crit = 2,   #critical suveillance time
  I_crit = 3,     #threshold for reporting
  pi_report = 1, #reporting probability
  pi_detect = 1, #detection probability
  cull_time = 1, #detection cull time
  network_type = "smallworld",
  network_parms = list(dim = 1, size = farm_number, nei = 2.33, p = 0.0596, multiple = FALSE, loops = FALSE),
  stochastic_network = TRUE
)

sims <- 101

move_list <- mclapply(seq_len(sims), function(num){
  conditions <- create_initial_condition(basic_patches(farm_size, farm_number))
  return(mf_sim(init = conditions, parameters = move_parms, times=1:365, n_sims = 1))
}, mc.cores = detectCores())

move_results <- do.call("abind", sim_list)

test_that("Chickens are not disappearing", {
  truths <- apply(move_results, 4, function(x){
    total_r <- colSums(x["R",,])
    inconsistencies <- which(total_r < dplyr::lag(total_r,1))
    if (!length(inconsistencies)){
      return(TRUE)
    }
    culled <- rowSums(x["C",,])
    patch_nums <- as.numeric(names(which(culled > 0)))
    culled_patches <- x["C",patch_nums,]
    if (is.null(dim(culled_patches)[1])){
      return(TRUE)
    }
    times <- apply(culled_patches, 1, function(y) which.max(y == 2))
    return(inconsistencies %in% times)
  })
  expect_true(all(unlist(truths)))
})
