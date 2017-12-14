// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(RcppProgress)]]
// [[Rcpp::plugins(cpp11)]]

#include "includes.h"
using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
arma::mat sim_gillespie(arma::vec &init, const List parmlist, const arma::vec &times, const bool &progress, const bool &abort) {
  mf_parmlist parms(parmlist);

  mf_vals vals;
  arma::vec state = init;
  //Rcout << state << std::endl;
  arma::vec rates = set_rates(state, parms, vals);
  arma::vec cumrates = cumsum(rates);
  double sum_rates = cumrates(cumrates.n_elem - 1);
  arma::mat results = zeros<arma::mat>(init.n_elem, times.n_elem);
  arma::vec results_vec = arma::vec(results.memptr(), results.n_elem, false, true);

  double time = times(0);
  double time_next;
  arma::vec::const_iterator next_record_time = times.begin();
  uword record = 0;
  double time_max = as_scalar(times.tail(1));
  bool abort_yes = false;
  uword event;

  results.col(record) = state;

  Progress p(times.n_elem, progress);

  while(*(next_record_time) < time_max) {
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
      //Rcpp::Rcout << abort << std::endl;
      if(abort) {
        abort_yes = abort_criterion(init, state, parms, vals);
      }
      //Rcpp::Rcout << abort_yes << std::endl;
    }

    if(sum_rates == 0 || abort_yes) break;
    if(!is_finite(time)) {
      std::stringstream err_msg;
      err_msg<< "Time has an invalid value: " << time << std::endl;
      stop(err_msg.str());
    }

    sample_cumrates(event, cumrates, sum_rates);
    update_state(state, parms, time, vals, event);
    update_rates(rates, cumrates, state, parms, vals, event, time);
    sum_rates = cumrates(cumrates.n_elem - 1);

    //Rcout << std::endl;
    time = time_next;
  }

    return results_vec;
};
