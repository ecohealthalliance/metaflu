#include "includes.h"
using namespace Rcpp;
using namespace arma;

// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(RcppProgress)]]
// [[Rcpp::export]]
arma::cube sim_gillespie(const arma::mat &init, const List parmlist, const arma::vec &times, const bool &progress) {

  mf_parmlist parms(parmlist);
  arma::mat state = init;
  arma::cube results = zeros<arma::cube>(init.n_rows, init.n_cols, times.n_elem);

  double time = times(1);
  double time_next;
  arma::vec::const_iterator next_record_time = times.begin();
  uword record = 0;
  double time_max = as_scalar(times.tail(1));

  arma::mat rates = arma::zeros<arma::mat>(state.n_rows, parms.n_actions);
  uword event;
  uword action;
  uword patch;

  results.slice(record) = state;

  Progress p(times.n_elem, progress);

  while(*(next_record_time) < time_max) {
    update_rates(rates, state, parms);
    double sum_rates = accu(rates);
    time_next = time + R::rexp(sum_rates);

    while(time_next > *(next_record_time)) {  //only store times on the grid
      Progress::check_abort();
      p.update(time);
      results.slice(record) = state;
      ++(next_record_time);
      ++record;
    }

    event = as_scalar(arma::find(R::runif(0,1) < cumsum(vectorise(rates) / sum_rates), 1, "last"));
    patch = event % parms.K;
    action = event / parms.K;
    update_state(state, parms, action, patch);
    time = time_next;
  }
  return results;
};