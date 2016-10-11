#' Influenza Metapopulation Simulator
#' @export
mf_sim <- function(init, parameters, times=NULL, n_sims=1, seed=NULL, .parallel=FALSE, .progress=TRUE) {
  progress = progress & !parallel

  old_rng_kind = RNGkind("L'Ecuyer-CMRG")
  on.exit(RNGkind(old_rng_kind[1], old_rng_kind[2]))
  if(!is.null(seed)) set.seed(seed)

  results = sim_gillespie(init, parameters, times, progress)

  return(results)
}
