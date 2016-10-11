#include "includes.h"
using namespace Rcpp;
using namespace arma;

// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::depends(RcppArmadillo)]]

void update_rates(arma::mat &rates, const arma::mat &state, const mf_parmlist &parms) {
    //arma::mat rates = arma::zeros<arma::mat>(state.n_rows, parms.n_actions);
    rates.col(1) = state.col(1) % (parms.beta *  state.col(2) + parms.nu * ( 1 - exp(-parms.phi * state.col(4)))); //infection events
    rates.col(2) = state.col(2) * parms.gamma; //recovery events
    rates.col(3) = state.col(1) * parms.mu;    //susceptible mortality
    rates.col(4) = state.col(2) * (parms.mu + parms.alpha); //infected mortality
    rates.col(5) = state.col(3) * parms.mu;     //recovered mortality
    rates.col(6) = state.col(2) * parms.sigma;   //virion shedding
    rates.col(7) = state.col(4) * parms.eta;     //viron degradation
    rates.cols(8, 10) = state.cols(1, 3) * parms.omega; // emigration of all classes
    //return(rates);
};

arma::Row<sword> act(uword action) {
  arma::Mat<sword> action_matrix = {
                                    {-1,  1,  0,  0},  //1 infection
                                    { 0, -1, +1,  0},  //2 recovery
                                    {-1,  0,  0,  0},  //3 S mortality
                                    { 0, -1,  0,  0},  //4 I mortality
                                    { 0,  0, -1,  0},  //5 R mortality
                                    { 0,  0,  0, +1},  //6 shedding
                                    { 0,  0,  0, -1},  //7 degradation
                                    {-1,  0,  0,  0},  //8 S emigration
                                    { 0, -1,  0,  0},  //9 I emigration
                                    { 0,  0, -1,  0},  //10 R emigration
  };
  return action_matrix.row(action);
};

void update_state(arma::mat &state, const mf_parmlist &parms, arma::uword action, arma::uword patch) {
  state.row(patch) += state.row(patch) + act(action);  //Apply the given action to the patch
  if(action >= 8) {    //If an emigration action, send the individual to another patch
    double rand_value = R::runif(0,1);
    uword target_patch = as_scalar(find(parms.chi_cum.row(patch) > rand_value, 1, "last"));
    state.row(target_patch) = state.row(target_patch) - act(action);
  }
};