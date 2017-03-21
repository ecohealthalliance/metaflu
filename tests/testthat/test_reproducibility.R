suppressWarnings(suppressMessages({
  library(dplyr)
  library(foreach)
  library(doMC)
}))

try_parms = list(
  beta = 0.004,                     #contact rate for direct transmission
  gamma = 0.167,                    #recovery rate
  mu = 0,                           #base mortality rate
  alpha = 0,                        #disease mortality rate
  phi = 1.96e-4,                    #infectiousness of environmental virions
  eta = 0.14,                       #degradation rate of environmental virions
  nu =  0.001,                      #uptake rate of environmental virion
  sigma = 0,                        #virion shedding rate
  omega = 0,                        #movement rate
  rho = 0,                          #contact rate nonlinearity 0=density-dependent, 1=frequency-dependent
  lambda = 0,                       #force of infection from external sources
  tau_crit = 0,                     #critical suveillance time
  I_crit = 0,                       #threshold for reporting
  pi_report = 0,                    #reporting probability
  pi_detect = 0,                    #detection probability
  chi = NULL,
  network_type = "smallworld",
  network_parms = list(dim = 1, size = 20, nei = 2.33, p = 0.0596, multiple = FALSE, loops = FALSE),
  stochastic_network = TRUE
)

set.seed(123)
initial_cond <- matrix(c(100, 0, 0, 0), nrow=20, ncol=4, byrow=TRUE)
infected_patches <- sample(seq_len(nrow(initial_cond)), 2)
initial_cond[infected_patches, 2] <- 1
initial_cond[infected_patches, 1] <- initial_cond[infected_patches, 1] - 1

# registerDoSEQ()
# set.seed(123)
# rdsOutput <- mf_sim(init = initial_cond, parameters = try_parms, times=0:1000, n_sims = 2)
# saveRDS(rdsOutput,"tests/testthat/stored_sim_results.rds")

context("Reproducibility")

test_that("outputs with same seed are identical regardless of parallelism", {
  registerDoSEQ()
  set.seed(123)
  output1 <- mf_sim(init = initial_cond, parameters = try_parms, times=0:1000, n_sims = 5)
  registerDoMC(cores = 2)
  set.seed(123)
  output2 <- mf_sim(init = initial_cond, parameters = try_parms, times=0:1000, n_sims = 5)
  expect_identical(output1, output2)
})

test_that("sims within run are not identical", {
  output1 <- mf_sim(init = initial_cond, parameters = try_parms, times=0:1000, n_sims = 3)
  sim1 <- output1 %>%
    filter(sim == 1)
  sim2 <- output1 %>%
    filter(sim == 2)
  sim3 <- output1 %>%
    filter(sim == 3)
  expect_false(isTRUE(identical(sim1,sim2)))
  expect_false(isTRUE(identical(sim2,sim3)))
  expect_false(isTRUE(identical(sim1,sim3)))
})

test_that("patches are not identical", {
  patch3_parms = list(
    beta = 0.004,                     #contact rate for direct transmission
    gamma = 0.167,                    #recovery rate
    mu = 0,                           #base mortality rate
    alpha = 0,                        #disease mortality rate
    phi = 1.96e-4,                    #infectiousness of environmental virions
    eta = 0.14,                       #degradation rate of environmental virions
    nu =  0.001,                      #uptake rate of environmental virion
    sigma = 0,                      #virion shedding rate
    omega = 0,                        #movement rate
    rho = 0,                          #contact rate nonlinearity 0=density-dependent, 1=frequency-dependent
    lambda = 0,                      #force of infection from external sources
    tau_crit = 0,                     #critical suveillance time
    I_crit = 0,                       #threshold for reporting
    pi_report = 0,                    #reporting probability
    pi_detect = 0,                    #detection probability
    chi = matrix(c(1,0,0,0,1,0,0,0,1), nrow=3)  #patch connectivity matrix
  )
  initial_cond_patch3 <- matrix(c(99, 1, 0, 0), nrow=3, ncol=4, byrow=TRUE)
  output1 <- mf_sim(init = initial_cond_patch3, parameters = patch3_parms, times=0:1000, n_sims = 1)
  p1 <- output1 %>%
    filter(patch == 1)
  p2 <- output1 %>%
    filter(patch == 2)
  p3 <- output1 %>%
    filter(patch == 3)
  expect_false(isTRUE(identical(p1,p2)))
  expect_false(isTRUE(identical(p2,p3)))
  expect_false(isTRUE(identical(p1,p3)))
})

test_that("non-seeded runs after same-seeded runs are not identical", {
  registerDoMC(cores = 2)
  set.seed(123)
  output1 <- mf_sim(init = initial_cond, parameters = try_parms, times=0:1000, n_sims = 2)
  output1a <- mf_sim(init = initial_cond, parameters = try_parms, times=0:1000, n_sims = 2)
  output2 <- mf_sim(init = initial_cond, parameters = try_parms, times=0:1000, n_sims = 2)
  output2a <- mf_sim(init = initial_cond, parameters = try_parms, times=0:1000, n_sims = 2)
  expect_false(isTRUE(identical(output1a,output2a)))
  #non-parallel check
  registerDoSEQ()
  set.seed(123)
  output1 <- mf_sim(init = initial_cond, parameters = try_parms, times=0:1000, n_sims = 2)
  set.seed(123)
  output1a <- mf_sim(init = initial_cond, parameters = try_parms, times=0:1000, n_sims = 2)
  output2 <- mf_sim(init = initial_cond, parameters = try_parms, times=0:1000, n_sims = 2)
  output2a <- mf_sim(init = initial_cond, parameters = try_parms, times=0:1000, n_sims = 2)
  expect_false(isTRUE(identical(output1a,output2a)))
})

test_that("RDS-saved model same as model created with same settings", {
  set.seed(123)
  output3 <- mf_sim(init = initial_cond, parameters = try_parms, times=0:1000, n_sims = 2)
  load1 <- readRDS("stored_sim_results.rds")
  expect_identical(output3, load1)
})
