#' Influenza Metapopulation Simulator
#' @param init Intital conditions - a patches X 4 (SIRV) matrix of populations
#' @param parameters a list of parameter values.  See below for details
#' @param times a vector of output times to report
#' @param n_sims number of simulations to run
#' @param parallel run the simulations in parallel?
#' @param progress TRUE/FALSE show a progress bar?
#' @importFrom dplyr rename_ mutate_ recode select_ arrange_ bind_rows
#' @importFrom tibble as_tibble
#' @importFrom parallel mclapply
#' @export
mf_sim <- function(init, parameters, times, n_sims=1, seed=NULL, parallel=FALSE, progress=FALSE) {
  progress = progress & !parallel  #progress bars don't work in parallel

  old_rng_kind = RNGkind("L'Ecuyer-CMRG")  # set a parallel-friendly RNG generator
  on.exit(RNGkind(old_rng_kind[1], old_rng_kind[2]), add=TRUE)

  if(!is.null(seed)) {
    old_seed <- .Random.seed
    set.seed(seed)
    on.exit(assign(".Random.seed", old_seed, envir=globalenv()), add=TRUE)
  }


  if(parallel) n_cores = getOption("mc.cores", 2L) else n_cores = 1

  results = mclapply(seq_len(n_sims), function(sim) {
    reshape2::melt(sim_gillespie(init=init, parmlist=parameters, times=times, progress=progress))
  }, mc.set.seed=TRUE, mc.silent = TRUE)

  results = as_tibble(bind_rows(results, .id = "sim"))
  results = rename_(results, .dots = c("patch"="Var1", "class"="Var2", "time"="Var3", "population"="value"))
#  results = mutate_(results, class = ~as.character(class))
  results = mutate_(results, class =  ~recode(class, `1`="S", `2`="I", `3`="R", `4`="V"))
  results = mutate_(results, class =  ~factor(class, levels=(c("S", "I", "R", "V"))))
  results = select_(results, "sim", "time", "patch", "class", "population")
  results = arrange_(results, "sim", "time", "patch", "class")
  return(results)
}

# TODO: option to output all event times, not a grid