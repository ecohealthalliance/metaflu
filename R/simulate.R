#' Influenza Metapopulation Simulator
#' @param init Intital conditions - a patches X 4 (SIRV) matrix of populations
#' @param parameters a list of parameter values.  See below for details
#' @param times a vector of output times to report
#' @param n_sims number of simulations to run
#' @param parallel run the simulations in parallel?
#' @param progress TRUE/FALSE show a progress bar?
#' @importFrom dplyr rename_ mutate_ recode select_ arrange_
#' @importFrom tibble as_tibble
#' @export
mf_sim <- function(init, parameters, times, n_sims=1, seed=NULL, parallel=FALSE, progress=TRUE) {
  progress = progress & !parallel  #progress bars don't work in parallel

  old_rng_kind = RNGkind("L'Ecuyer-CMRG")  # set a parallel-friendly RNG generator
  on.exit(RNGkind(old_rng_kind[1], old_rng_kind[2])) # change it back when done

  if(!is.null(seed)) set.seed(seed)

  results = sim_gillespie(init=init, parmlist=parameters, times=times, progress=progress)

  results = tibble::as_tibble(reshape2::melt(results))
  results = rename_(results, .dots = c("patch"="Var1", "class"="Var2", "time"="Var3", "population"="value"))
  results = mutate_(results, class = ~as.character(class))
  results = mutate_(results, class =  ~recode(class, `1`="S", `2`="I", `3`="R", `4`="V"))
  results = mutate_(results, class =  ~factor(class, levels=(c("S", "I", "R", "V"))))
  results = select_(results, "time", "patch", "class", "population")
  results = arrange_(results, "time", "patch", "class")
  return(results)
}

# TODO: option to output all event times, not a grid