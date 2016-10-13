#ifndef INCLUDES
#define INCLUDES

#include <RcppArmadillo.h>
#include <progress.hpp>
#include "data_structures.h"

void update_rates(arma::mat &rates, const arma::mat &state, const mf_parmlist &parms);
arma::Row<arma::sword> act(arma::uword action);
void update_state(arma::mat &state, const mf_parmlist &parms, arma::uword action, arma::uword patch);

#endif