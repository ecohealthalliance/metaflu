// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(RcppProgress)]]
// [[Rcpp::plugins(cpp11)]]

#include "includes.h"
using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
arma::cube sim_gillespie(const arma::mat &init, const List parmlist, const arma::vec &times, const bool &progress) {
  mf_parmlist parms(parmlist);
  arma::mat state = init;
  arma::cube results = zeros<arma::cube>(init.n_rows, init.n_cols, times.n_elem);
  double time = times(0);
  double time_next;
  arma::vec::const_iterator next_record_time = times.begin();
  uword record = 0;
  double time_max = as_scalar(times.tail(1));
  arma::mat rates = arma::zeros<arma::mat>(state.n_rows, action_matrix.n_rows);
  uword event;
  uword action;
  uword patch;

 // double randnum;
//  arma::vec vvec;
  results.slice(record) = state;
  Progress p(times.n_elem, progress);
  while(*(next_record_time) < time_max) {

    update_rates(rates, state, parms);
    double sum_rates = accu(rates);
    if(sum_rates == 0) {
      time_next = time_max + 1;
    } else {
      time_next = time + R::rexp(1/sum_rates);
    }
//    Rcout << time << "; " << time_next  << "; " << *(next_record_time) << "; " << time_max  << "; " << sum_rates << std::endl;
//    Rcout << rates << std::endl;
//    Rcout << sum_rates << ";  " << time << "; " << time_next << std::endl;
      while(time_next > *(next_record_time) && *(next_record_time) < time_max) {  //only store times on the grid
      Progress::check_abort();
      p.update(time);
      ++(next_record_time);
      ++record;
      results.slice(record) = state;
    }
    if(sum_rates == 0) break;
    event = as_scalar(arma::find(cumsum(vectorise(rates) / sum_rates) > R::runif(0,1), 1));
    patch = event % parms.K;
    action = event / parms.K;
    update_state(state, parms, action, patch);
    time = time_next;
  }
  return results;
};
