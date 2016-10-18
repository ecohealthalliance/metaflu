#' Influenza Metapopulation Simulator
#'
#' The simulator uses the \strong{foreach} and \strong{doRNG} packages to ensure
#' reproducible, parallelizable results.  To run code in parallel, register
#' a \strong{foreach} back-end such as \code{\link[doMC]{registerDoMC}}.
#' @param init Intital conditions - a patches X 4 (SIRV) matrix of populations
#' @param parameters a list of parameter values.  See below for details
#' @param times a vector of output times to report
#' @param n_sims number of simulations to run
#' @examples
#'  library(doMC)
#'  registerDoMC(cores=2)
#'  parms = list(
#'    beta = 0.004,                     #contact rate for direct transmission
#'    gamma = 0.167,                    #recovery rate
#'    mu = 0,                           #base mortality rate
#'    alpha = 0,                        #disease mortality rate
#'    phi = 1.96e-4,                    #infectiousness of environmental virions
#'    eta = 0.14,                       #degradation rate of environmental virions
#'    nu =  0.001,                      #uptake rate of environmental virion
#'    sigma = 0,                        #virion shedding rate
#'    omega = 0,                        #movement rate
#'    rho = 0,                          #contact rate nonlinearity 0=density-dependent, 1=frequency-dependent
#'    lambda = 0,                      #force of infection from external sources
#'    chi = matrix(c(1,0,0,1), nrow=2)  #patch connectivity matrix
#'    )
#'  initial_cond <- matrix(c(99, 1, 0, 0), nrow=2, ncol=4, byrow=TRUE)
#'  output <- mf_sim(init = initial_cond, parameters = parms, times=0:1000, n_sims = 2)
#' @importFrom dplyr rename_ mutate_ recode select_ arrange_ bind_rows
#' @importFrom tibble as_tibble
#' @importFrom foreach foreach %dopar%
#' @importFrom doRNG %dorng%
#' @export
mf_sim <- function(init, parameters, times, n_sims=1) {

  if(!is.null(parameters[["network_type"]]) && !parameters[["stochastic_network"]]) {
    parameters[["chi"]] <- make_net(parameters[["network_type"]],
                                    parameters[["network_parms"]])
  }

  sim_fun <- function() {
    if(!is.null(parameters[["network_type"]]) && parameters[["stochastic_network"]]) {
      parameters[["chi"]] <- make_net(parameters[["network_type"]],
                                      parameters[["network_parms"]])
    }
    return(reshape2::melt(sim_gillespie(init=init, parmlist=parameters, times=times, progress=FALSE)))
  }

  suppressWarnings(suppressMessages({
    results = foreach(i=seq_len(n_sims)) %dorng% { sim_fun() }
  }))

  results = as_tibble(bind_rows(results, .id = "sim"))
  results = rename_(results, .dots = c("patch"="Var1", "class"="Var2", "time"="Var3", "population"="value"))
#  results = mutate_(results, class = ~as.character(class))
  results = mutate_(results, class =  ~recode(class, `1`="S", `2`="I", `3`="R", `4`="V"))
  results = mutate_(results, class =  ~factor(class, levels=(c("S", "I", "R", "V"))))
  results = select_(results, "sim", "time", "patch", "class", "population")
  results = arrange_(results, "sim", "time", "patch", "class")
  return(results)
}

#' @import igraph
make_net <- function(network_type, network_parms) {
  net = as_adj(do.call(paste0("sample_", network_type), network_parms),
               sparse=FALSE)
  net = net/rowSums(net)

}


# TODO: option to output all event times, not a grid
#