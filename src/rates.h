#ifndef RATES
#define RATES


#include "includes.h"
using namespace Rcpp;
using namespace arma;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]

arma::vec set_rates(arma::vec &state, const mf_parmlist &parms, mf_vals &vals) {
    arma::vec rates(parms.K * parms.n_actions);
    vals.rates_mat = arma::mat(rates.memptr(), parms.n_actions, parms.K, false, false);
    vals.state_mat = arma::mat(state.memptr(), parms.n_pstates, parms.K, false, false);

    vals.rates_mat.row(0) = vals.state_mat.row(0) % (parms.beta *  vals.state_mat.row(1) / pow(sum(vals.state_mat.rows(0, 2), 0), parms.rho) + //infection events (contact)
                                   parms.nu * ( 1 - exp(-parms.phi * vals.state_mat.row(3))) + //virion uptake
                                   parms.lambda); //external infection

    vals.rates_mat.row(1) = vals.state_mat.row(1) * parms.gamma; //recovery events
    vals.rates_mat.row(2) = vals.state_mat.row(0) * parms.mu;    //susceptible mortality
    vals.rates_mat.row(3) = vals.state_mat.row(1) * (parms.mu + parms.alpha); //infected mortality
    vals.rates_mat.row(4) = vals.state_mat.row(2) * parms.mu;     //recovered mortality
    vals.rates_mat.row(5) = vals.state_mat.row(1) * parms.sigma;   //virion shedding
    vals.rates_mat.row(6) = vals.state_mat.row(3) * parms.eta;     //viron degradation
    vals.rates_mat.rows(7, 9) = vals.state_mat.rows(0, 2) * parms.omega; // emigration of all classes

    //Rcout << rates.t();
    return rates;

};

void update_state(arma::vec &state, const mf_parmlist &parms, mf_vals &vals, arma::uword event) {

  vals.patch = event / parms.n_actions;
  vals.action = event % parms.n_actions;

  state.subvec(vals.patch * 4, vals.patch * 4 + 3) = state.subvec(vals.patch * 4, vals.patch * 4 + 3) + parms.action_matrix.col(vals.action);  //Apply the given action to the patch
  //Rcout << state.row(patch)<< std::endl;
  if(vals.action >= 7) {
    double rand_value = R::runif(0,1);
    vals.target_patch = as_scalar(find(parms.chi_cum.row(vals.patch) > rand_value, 1, "first"));
 //   Rcout << vals.patch << " " << vals.target_patch <<  " " << vals.action << " " << std::endl;
    state.subvec(vals.target_patch * 4, vals.target_patch * 4 + 3) = state.subvec(vals.target_patch * 4, vals.target_patch * 4 + 3) - parms.action_matrix.col(vals.action);
  }
};

void update_rates(arma::vec &rates, const arma::vec &state, const mf_parmlist &parms,  mf_vals &vals, arma::uword event) {

 // Rcout << vals.patch << ";  " << vals.action<< std::endl;
  vals.rates_mat(0, vals.patch) = vals.state_mat(0, vals.patch) *
    (parms.beta *  vals.state_mat(1, vals.patch) / pow(accu(vals.state_mat.submat(0, vals.patch, 2, vals.patch)), parms.rho)) + //infection events (contact)
    parms.nu * ( 1 - exp(-parms.phi * vals.state_mat(3, vals.patch))) + //virion uptake
    parms.lambda; //external infection
  vals.rates_mat(1, vals.patch) = vals.state_mat(1, vals.patch) * parms.gamma; //recovery events
  vals.rates_mat(2, vals.patch) = vals.state_mat(0, vals.patch) * parms.mu;    //susceptible mortality
  vals.rates_mat(3, vals.patch) = vals.state_mat(1, vals.patch) * (parms.mu + parms.alpha); //infected mortality
  vals.rates_mat(4, vals.patch) = vals.state_mat(2, vals.patch) * parms.mu;     //recovered mortality
  vals.rates_mat(5, vals.patch) = vals.state_mat(1, vals.patch) * parms.sigma;   //virion shedding
  vals.rates_mat(6, vals.patch) = vals.state_mat(3, vals.patch) * parms.eta;     //viron degradation
  vals.rates_mat.submat(7, vals.patch, 9, vals.patch) = vals.state_mat.submat(0, vals.patch, 2, vals.patch) * parms.omega; // emigration of all classes

  if(vals.action >= 7) {    //If an emigration action, send the individual to another patch
    vals.rates_mat(0, vals.target_patch) = vals.state_mat(0, vals.target_patch) *
      (parms.beta *  vals.state_mat(1, vals.target_patch) / pow(accu(vals.state_mat.submat(0, vals.target_patch, 2, vals.target_patch)), parms.rho)) + //infection events (contact)
      parms.nu * ( 1 - exp(-parms.phi * vals.state_mat(3, vals.target_patch))) + //virion uptake
      parms.lambda; //external infection
    vals.rates_mat(1, vals.target_patch) = vals.state_mat(1, vals.target_patch) * parms.gamma; //recovery events
    vals.rates_mat(2, vals.target_patch) = vals.state_mat(0, vals.target_patch) * parms.mu;    //susceptible mortality
    vals.rates_mat(3, vals.target_patch) = vals.state_mat(1, vals.target_patch) * (parms.mu + parms.alpha); //infected mortality
    vals.rates_mat(4, vals.target_patch) = vals.state_mat(2, vals.target_patch) * parms.mu;     //recovered mortality
    vals.rates_mat(5, vals.target_patch) = vals.state_mat(1, vals.target_patch) * parms.sigma;   //virion shedding
    vals.rates_mat(6, vals.target_patch) = vals.state_mat(3, vals.target_patch) * parms.eta;     //viron degradation
    vals.rates_mat.submat(7, vals.target_patch, 9, vals.target_patch) = vals.state_mat.submat(0, vals.target_patch, 2, vals.target_patch) * parms.omega; // emigration of all classes
  }



}
#endif
