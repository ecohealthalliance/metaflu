library(metaflu)

l <- readRDS(here::here("inst/testnets.rds"))

network <- as.matrix(l[[1]])

inits <- l[[2]]

# patch 1275 starts with infection

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
  tau_crit = 5,   #critical surveillance time
  I_crit = 5,     #threshold for reporting
  pi_report = 0, #reporting probability
  pi_detect = 0, #detection probability
  cull_time = 10,#time to culling
  network_type = "smallworld",
  network_parms = list(dim = 1, size = dim(inits)[1], nei = 2.33, p = 0.0596, multiple = FALSE, loops = FALSE),
  stochastic_network = FALSE
)

parms$chi <- network

res <- mf_sim(init = inits, parameters = parms, times = 1:28, n_sims = 10)

View(res["I",,,1])

res_abort <- mf_sim(init = inits, parameters = parms, times = 1:28, n_sims = 10, abort = TRUE)

View(res_abort["I",,,1])