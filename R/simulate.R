#' Influenza Metapopulation Simulator
#'
#' The simulator uses the **foreach** and **doRNG** packages to ensure
#' reproducible, parallelizable results.  To run code in parallel, register
#' a **foreach** back-end such as [doMC::registerDoMC()].
#' @param init Initial conditions - a patches X 4 (SIRV) matrix of populations
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
#'    cull_time = 0,    #average time, in days, it takes for a reported patch to be culled
#'    chi = matrix(c(1,0,0,1), nrow=2)  #patch connectivity matrix
#'    )
#'  initial_cond <- matrix(c(99, 1, 0, 0), nrow=2, ncol=4, byrow=TRUE)
#'  output <- mf_sim(init = initial_cond, parameters = parms, times=0:1000, n_sims = 2)
#' @importFrom dplyr as_data_frame lst_
#' @importFrom tidyr crossing_
#' @importFrom stringi stri_subset_regex
#' @importFrom foreach foreach
#' @importFrom doRNG %dorng%
#' @importFrom purrr transpose
#' @export
mf_sim <- function(init, parameters, times, n_sims=1, return_array=TRUE) {

  times <- as.double(times)
  init <- as.vector(t(init))

  #Get number of patches
  if(!is.null(parameters[["network_parms"]])) {
    n_patches = parameters[["network_parms"]][["size"]]
  } else {
    n_patches = nrow(parameters[["chi"]])
  }

  #Convert scalar parameters to vectors
  parameters <-  repeat_parms(parameters, n_patches, c("lambda", "tau_crit", "I_crit", "pi_report", "pi_detect", "cull_time"))

  if(!is.null(parameters[["network_type"]]) && !parameters[["stochastic_network"]]) {
    parameters[["chi"]] <- make_net(parameters[["network_type"]],
                                    parameters[["network_parms"]])
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
    outputs = foreach(i=seq_len(n_sims),
                      .options.multicore=list(preschedule= FALSE)) %dorng% {
      sim_fun()
      }
  }))
  outputs <- purrr::transpose(outputs)
  networks <- outputs[[2]]
  names(networks) <- seq_len(n_sims)

  if(return_array) {
    results <- array(data = unlist(outputs[[1]], recursive = FALSE, use.names = FALSE),
                     dim = c(5, n_patches, length(times), n_sims),
                     dimnames =  list(c("S", "I", "R", "V", "C"),
                                      seq_len(n_patches),
                                      times,
                                      seq_len(n_sims)))
  } else {
  results <- crossing_(list(
    sim = seq_len(n_sims),
    time = times,
    patch = seq_len(n_patches),
    class = factor(c("S", "I", "R", "V"), levels = c("S", "I", "R", "V", "C"))
  ))
  results[["population"]] <- unlist(outputs[[1]], recursive = FALSE, use.names = FALSE)
  }
  if(!is.null(parameters[["network_type"]]) && parameters[["stochastic_network"]]) {
    attr(results, "networks") <- networks
  } else {
    attr(results, "network") <- parameters[["chi"]]
  }
  return(results)
}

#' @importFrom igraph as_adj
#' @export
make_net <- function(network_type, network_parms) {
  net_fun <- get(paste0("sample_", network_type), asNamespace("igraph"))

  net = as_adj(do.call(net_fun, network_parms),
               sparse=FALSE)
  net = net/rowSums(net)

}


repeat_parms <- function(parameters, n_patches, parms_to_extend) {
  for(parm in parms_to_extend) {
    parameters[[parm]] = rep_len(parameters[[parm]], n_patches)
  }
  return(parameters)
}
