suppressWarnings(suppressMessages(library(tidyverse)))

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
  chi = matrix(c(1,0,0,1), nrow=2)  #patch connectivity matrix
)

initial_cond <- matrix(c(99, 1, 0, 0,
                         99, 1, 0, 0),
                       nrow=2, byrow=TRUE)

rdsOutput <- mf_sim(init = initial_cond, parameters = try_parms, times=0:1000, n_sims = 2, seed = 2)
saveRDS(rdsOutput,"temp.rds")

context("Reproducibility")

test_that("outputs with same seed are identical regardless of parallelism", {
  options(mc.cores=2)
  output1 <- mf_sim(init = initial_cond, parameters = try_parms, times=0:1000, n_sims = 2, seed = 2, parallel = FALSE)
  output2 <- mf_sim(init = initial_cond, parameters = try_parms, times=0:1000, n_sims = 2, seed = 2, parallel = TRUE)
  expect_identical(output1, output2)
})

test_that("sims within run are not identical", {
  output1 <- mf_sim(init = initial_cond, parameters = try_parms, times=0:1000, n_sims = 3, seed = 2, parallel = FALSE)
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
    chi = matrix(c(1,0,0,0,1,0,0,0,1), nrow=3)  #patch connectivity matrix
  )
  output1 <- mf_sim(init = initial_cond, parameters = patch3_parms, times=0:1000, n_sims = 1, seed = 2, parallel = FALSE)
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
  options(mc.cores=2)
  output1 <- mf_sim(init = initial_cond, parameters = try_parms, times=0:1000, n_sims = 2, seed = 2, parallel = TRUE)
  output1a <- mf_sim(init = initial_cond, parameters = try_parms, times=0:1000, n_sims = 2, parallel = TRUE)
  output2 <- mf_sim(init = initial_cond, parameters = try_parms, times=0:1000, n_sims = 2, seed = 2, parallel = TRUE)
  output2a <- mf_sim(init = initial_cond, parameters = try_parms, times=0:1000, n_sims = 2, parallel = TRUE)
  expect_false(isTRUE(identical(output1a,output2a)))
  #non-parallel check
  output1 <- mf_sim(init = initial_cond, parameters = try_parms, times=0:1000, n_sims = 2, seed = 2, parallel = FALSE)
  output1a <- mf_sim(init = initial_cond, parameters = try_parms, times=0:1000, n_sims = 2, parallel = FALSE)
  output2 <- mf_sim(init = initial_cond, parameters = try_parms, times=0:1000, n_sims = 2, seed = 2, parallel = FALSE)
  output2a <- mf_sim(init = initial_cond, parameters = try_parms, times=0:1000, n_sims = 2, parallel = FALSE)
  expect_false(isTRUE(identical(output1a,output2a)))
})

test_that("RDS-saved model same as model created with same settings", {
  output2 <- mf_sim(init = initial_cond, parameters = try_parms, times=0:1000, n_sims = 2, seed = 2)
  load1 <- readRDS("temp.rds")
  expect_identical(output2, load1)
})
