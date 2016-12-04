#' Influenza Metapopulation Simulator
#'
#' The simulator uses the **foreach** and **doRNG** packages to ensure
#' reproducible, parallelizable results.  To run code in parallel, register
#' a **foreach** back-end such as [doMC::registerDoMC()].
#' @param init Intital conditions - a patches X 4 (SIRV) matrix of populations
#' @param parameters a list of parameter values.  See below for details
#' @param times a vector of output times to report
#' @param n_sims number of simulations to run
#' @examples
#'  library(doMC)
#'  registerDoMC(cores=2)
#'  parms = list(
#'    beta = 0.004,   #contact rate for direct transmission
#'    gamma = 0.167,  #recovery rate
#'    mu = 0,         #base mortality rate
#'    alpha = 0,      #disease mortality rate
#'    phi = 1.96e-4,  #infectiousness of environmental virions
#'    eta = 0.14,     #degradation rate of environmental virions
#'    nu =  0.001,    #uptake rate of environmental virion
#'    sigma = 0,      #virion shedding rate
#'    omega = 0,      #movement rate
#'    rho = 0,        #contact  nonlinearity 0=dens-dependent, 1=freq-dependent
#'    lambda = 0,    #force of infection from external sources
#'    chi = matrix(c(1,0,0,1), nrow=2)  #patch connectivity matrix
#'    )
#'  initial_cond <- matrix(c(99, 1, 0, 0), nrow=2, ncol=4, byrow=TRUE)
#'  output <- mf_sim(init = initial_cond, parameters = parms, times=0:1000, n_sims = 2)
#' @importFrom dplyr as_data_frame lst_
#' @importFrom tidyr crossing_
#' @importFrom stringi stri_subset_regex
#' @importFrom foreach foreach %dopar%
#' @importFrom doRNG %dorng%
#' @importFrom purrr transpose
#' @export
mf_sim <- function(init, parameters, times, n_sims=1) {

  times <- as.double(times)
  init <- as.vector(t(init))
  if(!is.null(parameters[["network_type"]]) && !parameters[["stochastic_network"]]) {
    parameters[["chi"]] <- make_net(parameters[["network_type"]],
                                    parameters[["network_parms"]])
  }

  if(!is.null(parameters[["network_parms"]])) {
    n_patches = parameters[["network_parms"]][["size"]]
  } else {
    n_patches = nrow(parameters[["chi"]])
  }

  sim_fun <- function() {
    if(!is.null(parameters[["network_type"]]) && parameters[["stochastic_network"]]) {
      parameters[["chi"]] <- make_net(parameters[["network_type"]],
                                      parameters[["network_parms"]])
      net_gen <- parameters[["chi"]]
    } else {
      net_gen <- NULL
    }
    return(list(
      sim_gillespie(init=init, parmlist=parameters, times=times, progress=FALSE),
      net_gen))
  }

  suppressWarnings(suppressMessages({
    outputs = foreach(i=seq_len(n_sims)) %dorng% {
      sim_fun()
      }
  }))
  outputs <- purrr::transpose(outputs)
  networks <- outputs[[2]]
  names(networks) <- 1:n_sims
  results <- crossing_(list(
    sim = seq_len(n_sims),
    time = times,
    patch = seq_len(n_patches),
    class = factor(c("S", "I", "R", "V"), levels = c("S", "I", "R", "V"))
  ))
  results[["population"]] <- unlist(outputs[[1]], recursive = FALSE, use.names = FALSE)
  if(!is.null(parameters[["network_type"]]) && parameters[["stochastic_network"]]) {
    attr(results, "networks") <- networks
  } else {
    attr(results, "network") <- parameters[["chi"]]
  }
  return(results)
}

#' @importFrom igraph as_adj
make_net <- function(network_type, network_parms) {
  net_fun <- get(paste0("sample_", network_type), asNamespace("igraph"))

  net = as_adj(do.call(net_fun, network_parms),
               sparse=FALSE)
  net = net/rowSums(net)

}
