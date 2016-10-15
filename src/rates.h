#ifndef RATES
#define RATES


#include "includes.h"
using namespace Rcpp;
using namespace arma;

// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::depends(RcppArmadillo)]]

void update_rates(arma::mat &rates, const arma::mat &state, const mf_parmlist &parms) {
    rates.col(0) = state.col(0) % (parms.beta *  state.col(1) + parms.nu * ( 1 - exp(-parms.phi * state.col(3)))); //infection events
    rates.col(1) = state.col(1) * parms.gamma; //recovery events
    rates.col(2) = state.col(0) * parms.mu;    //susceptible mortality
    rates.col(3) = state.col(1) * (parms.mu + parms.alpha); //infected mortality
    rates.col(4) = state.col(2) * parms.mu;     //recovered mortality
    rates.col(5) = state.col(1) * parms.sigma;   //virion shedding
    rates.col(6) = state.col(3) * parms.eta;     //viron degradation
    rates.cols(7, 9) = state.cols(0, 2) * parms.omega; // emigration of all classes
};

arma::imat action_matrix = {
                                    {-1,  1,  0,  0},  //0 infection
                                    { 0, -1, +1,  0},  //1 recovery
                                    {-1,  0,  0,  0},  //2 S mortality
                                    { 0, -1,  0,  0},  //3 I mortality
                                    { 0,  0, -1,  0},  //4 R mortality
                                    { 0,  0,  0, +1},  //5 shedding
                                    { 0,  0,  0, -1},  //6 degradation
                                    {-1,  0,  0,  0},  //7 S emigration
                                    { 0, -1,  0,  0},  //8 I emigration
                                    { 0,  0, -1,  0},  //9 R emigration
};


void update_state(arma::mat &state, const mf_parmlist &parms, arma::uword action, arma::uword patch) {
  //Rcout << action << " " << patch << " " << act(action) << " " << state.row(patch)<< std::endl;
  state.row(patch) = state.row(patch) + action_matrix.row(action);  //Apply the given action to the patch
  //Rcout << state.row(patch)<< std::endl;
  if(action >= 7) {    //If an emigration action, send the individual to another patch
    double rand_value = R::runif(0,1);
    uword target_patch = as_scalar(find(parms.chi_cum.row(patch) > rand_value, 1, "last"));
    state.row(target_patch) = state.row(target_patch) - action_matrix.row(action);
  }
};

#endif
