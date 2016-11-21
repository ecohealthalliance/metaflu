// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(RcppProgress)]]
// [[Rcpp::plugins(cpp11)]]

#include "includes.h"
using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
arma::mat sim_gillespie(const arma::vec &init, const List parmlist, const arma::vec &times, const bool &progress) {
  mf_parmlist parms(parmlist);
  mf_vals vals;
  arma::vec state = init;
  arma::vec rates = set_rates(state, parms, vals);

  arma::vec cumrates = zeros<arma::vec>(rates.n_elem);

  arma::mat results = zeros<arma::mat>(init.n_elem, times.n_elem);
  double time = times(0);
  double time_next;
  arma::vec::const_iterator next_record_time = times.begin();
  uword record = 0;
  double time_max = as_scalar(times.tail(1));

  uword event;

  results.col(record) = state;

  Progress p(times.n_elem, progress);

  while(*(next_record_time) < time_max) {
    double sum_rates = accu(rates);
    if(sum_rates == 0) {
      time_next = time_max + 1;
    } else {
      time_next = time + R::rexp(1/sum_rates);
    }

    while(time_next > *(next_record_time) && *(next_record_time) < time_max) {  //only store times on the grid
      Progress::check_abort();
      p.update(time);
      ++(next_record_time);
      ++record;
      results.col(record) = state;
    }

    if(sum_rates == 0) break;

    cumrates = cumsum(rates) / sum_rates;
    event = as_scalar(find(cumrates >= R::runif(0,1), 1));

    update_state(state, parms, vals, event);
    update_rates(rates, state, parms, vals, event);
    //Rcout << time << "; " << event << "; " << vals.action << "; " << vals.patch << "; " << cumrates.t();
    //Rcout << state.t();

    time = time_next;
  }

  return join_rows(times, trans(results));
};
