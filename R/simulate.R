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
#' @importFrom dplyr mutate_ select_ arrange_ bind_rows as_data_frame
#' @importFrom tidyr separate_ gather_
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
      tibble::as_data_frame(sim_gillespie(init=init, parmlist=parameters, times=times, progress=FALSE)),
      net_gen))
  }

  suppressWarnings(suppressMessages({
    results = foreach(i=seq_len(n_sims)) %dorng% {
      sim_fun()
      }
  }))
  results <- purrr::transpose(results)
  networks <- results[[2]]
  names(networks) <- 1:n_sims
  results <- results[[1]]
  results = bind_rows(results, .id = "sim")
  names(results)[-1] <- c("time", paste(rep(1:n_patches, each=4), c("S", "I", "R", "V"),  sep="_"))
  results <- gather_(results, "class", "population", gather_cols = stri_subset_regex(names(results), "\\d_\\w"))
  results <- separate_(results, "class", into=c("patch", "class"))

  results = mutate_(results, class =  ~factor(class, levels=(c("S", "I", "R", "V"))))
  results = select_(results, "sim", "time", "patch", "class", "population")

  results = arrange_(results, "sim", "time", "patch", "class")
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
