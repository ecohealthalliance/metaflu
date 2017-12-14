library(metaflu)
library(igraph)
l <- readRDS(here::here("inst/testnets.rds"))

network <- as_adj(make_lattice(dim=1, length=10), sparse = FALSE)

inits <- matrix(0, nrow=10, ncol=5)
inits[, 1] <- 50
inits[1, 1] <- 0
inits[1, 2] <- 50
# patch 1275 starts with infection

parms = list(
  beta = 0.1,   #contact rate for direct transmission
  gamma = 0,  #recovery rate
  mu = 0,         #base mortality rate
  alpha = 0,      #disease mortality rate
  phi = 0,  #infectiousness of environmental virions
  eta = 0,     #degradation rate of environmental virions
  nu =  0.00,    #uptake rate of environmental virion
  sigma = 0,      #virion shedding rate
  omega = 0.005,   #movement rate
  rho = 0,        #contact  nonlinearity 0=dens-dependent, 1=freq-dependent
  lambda = 0,     #force of infection from external sources
  tau_crit = 1,   #critical surveillance time
  I_crit = 5,     #threshold for reporting
  pi_report = 0, #reporting probability
  pi_detect = 0, #detection probability
  cull_time = 10,#time to culling
#  network_type = "smallworld",
#  network_parms = list(dim = 1, size = dim(inits)[1], nei = 2.33, p = 0.0596, multiple = FALSE, loops = FALSE),
  stochastic_network = FALSE,
  abort=TRUE
)

parms$chi <- network

set.seed(0)
system.time(res <- mf_sim(init = inits, parameters = parms, times = 1:50, n_sims = 1))
res["I",,,1]

parms$abort <- FALSE

set.seed(0)
system.time(res <- mf_sim(init = inits, parameters = parms, times = 1:50, n_sims = 1))
res["I",,,1]
