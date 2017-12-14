#ifndef RATES
#define RATES


#include "includes.h"
#include <sstream>
using namespace Rcpp;
using namespace arma;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]
void sample_cumrates(uword &event, const arma::vec &cumrates, const double &sum_rates) {

  double rnum = R::runif(0, sum_rates);
  event = 0;
  while(cumrates[event] < rnum) {
    event++;
  }
};

void sample_cumrates_row(uword &event, const arma::rowvec &cumrates, const double &sum_rates) {
  double rnum = R::runif(0, sum_rates);
  event = 0;
  while(cumrates[event] < rnum) {
    event++;
  }
};


arma::vec set_rates(arma::vec &state, const mf_parmlist &parms, mf_vals &vals) {
  arma::vec rates(parms.K * parms.n_actions);
  vals.rates_mat = arma::mat(rates.memptr(), parms.n_actions, parms.K, false, false);
  vals.state_mat = arma::mat(state.memptr(), parms.n_pstates, parms.K, false, false);
  vals.init_I_ind = vals.state_mat.row(1) > 0;
  //Rcpp::Rcout << vals.init_I_ind << std::endl;
  vals.init_I2_ind = (sum(parms.chi.rows(find(vals.init_I_ind)), 0) + vals.init_I_ind) > 0;
  //Rcpp::Rcout << vals.init_I2_ind << std::endl;
  vals.deathtimes = std::vector<std::priority_queue<double, std::vector<double>, std::greater<double>>>(parms.K);
  vals.cull = false;

  vals.rates_mat.row(0) = vals.state_mat.row(0) % (parms.beta *  vals.state_mat.row(1) / pow(sum(vals.state_mat.rows(0, 2), 0) + (vals.state_mat.row(0) == 0), parms.rho) + //infection events (contact)
    parms.nu * ( 1 - exp(-parms.phi * vals.state_mat.row(3))) + //virion uptake
    parms.lambda); //external infection

  vals.rates_mat.row(1) = vals.state_mat.row(1) * parms.gamma; //recovery events
  vals.rates_mat.row(2) = vals.state_mat.row(0) * parms.mu;    //susceptible mortality
  vals.rates_mat.row(3) = vals.state_mat.row(1) * (parms.mu + parms.alpha); //infected mortality
  vals.rates_mat.row(4) = vals.state_mat.row(2) * parms.mu;     //recovered mortality
  vals.rates_mat.row(5) = vals.state_mat.row(1) * parms.sigma;   //virion shedding
  vals.rates_mat.row(6) = vals.state_mat.row(3) * parms.eta;     //viron degradation
  vals.rates_mat.rows(7, 9) = vals.state_mat.rows(0, 2) * parms.omega; // emigration of all classes
  vals.rates_mat.row(10) = (vals.state_mat.row(4) == 1) % parms.cull_rate; // culling

  //Rcout << rates.t();
  return rates;

};

void update_state(arma::vec &state, mf_parmlist &parms, const double &time, mf_vals &vals, arma::uword event) {

  vals.patch = event / parms.n_actions;
  vals.action = event % parms.n_actions;
  //Rcout << vals.action << "; " << vals.patch << "; " << state.t() << "; " << time << std::endl;
  if(vals.action == 2 || vals.action == 3 || vals.action ==4) {  //in the case of death,
    //  Rcout << "test: " << vals.deathtimes[vals.patch].size() << std::endl;
    vals.deathtimes[vals.patch].push(time);                          //record the time
    while(vals.deathtimes[vals.patch].top() < (time - parms.tau_crit[vals.patch])) {
      vals.deathtimes[vals.patch].pop();
    }  //drop any deaths further back than tau_crit

    //Rcout << "test: " << vals.deathtimes[vals.patch].size() << "; " << vals.deathtimes[vals.patch].top() << std::endl;
    //Rcout << vals.patch << "; " << parms.tau_crit[vals.patch] << "; " << vals.deathtimes[vals.patch].size() << "; " << parms.I_crit[vals.patch] << std::endl;
    if(vals.deathtimes[vals.patch].size() >= parms.I_crit[vals.patch] &&   //if I_crit deaths have occurred in the past tau_crit time...
       R::runif(0,1) < parms.pi_report[vals.patch]                 &&   //and the deaths are reported
       R::runif(0,1) < (1 - pow((1 - parms.pi_detect[vals.patch]), state[vals.patch * 4 + 1]))) {   //and disease is detected
      state(vals.patch * parms.n_pstates + 4) = 1;  //set the state to reported
      //state.subvec(vals.patch * parms.n_pstates, vals.patch * parms.n_pstates + 3) = zeros<vec>(4);  //Cull
      //vals.cull = true;
      while(vals.deathtimes[vals.patch].size() > 0) {
        vals.deathtimes[vals.patch].pop();
      }
    } else {
      for(uword i = 0; i < parms.n_pstates; i++) {
        state[vals.patch * parms.n_pstates + i] += parms.action_matrix.at(i, vals.action);  //Apply the given action to the patch
      }
    }
  } else {
    for(uword i = 0; i < parms.n_pstates; i++) {
      state[vals.patch * parms.n_pstates + i] += parms.action_matrix.at(i, vals.action);  //Apply the given action to the patch
    }
    if(vals.action == 10) { //for culling, in addition to updating the cull state (done via action matrix)
      state.subvec(vals.patch * parms.n_pstates, vals.patch * parms.n_pstates + 3) = zeros<vec>(4);  //kill 'em all
      parms.chi.col(vals.patch) = zeros<vec>(parms.chi.n_rows);  //remove patch from the marix by setting its column to zero.
      parms.chi.each_col() /= sum(parms.chi, 1); //renormalize the rows - divide each column element-wise by the row sums
      parms.chi_cum = cumsum(parms.chi, 1); //recalculate cumulative probabilities
      //Rcpp::Rcout << "Culled farm #" << vals.patch << " at time " << time << std::endl;
    }

    //Rcout << state.row(patch)<< std::endl;
    if(vals.action == 7 | vals.action == 8 | vals.action == 9)  {
      //double rand_value = R::runif(0,1);
      //vals.target_patch = as_scalar(find(parms.chi_cum.row(vals.patch) > rand_value, 1, "first"));
      sample_cumrates_row(vals.target_patch, parms.chi_cum.row(vals.patch), 1);
      for(uword i = 0; i < parms.n_pstates; i++) {
        state[vals.target_patch * parms.n_pstates + i] -= parms.action_matrix.at(i, vals.action);  //Apply the given action to the patch
      }
      if(vals.action == 8){
        //std::cout << "An infected chicken has moved!" << std::endl;
      }

      }
  }


  // if(any(vals.rates_mat.col(vals.patch) < 0) || any(vals.state_mat.col(vals.patch) < 0)) {
  //   std::stringstream err_string;
  //   err_string << "Negative states or rates " << time << "; " << event << "; " << vals.patch <<  "; " << vals.action << "; " <<  std::endl << std::endl <<
  //     vals.state_mat << std::endl <<
  //       vals.rates_mat << std::endl;
  //   err_string << vals.state_mat << std::endl;
  //   Rcpp::stop(err_string.str());
  // }
};

void update_rates(arma::vec &rates, arma::vec &cumrates, const arma::vec &state, const mf_parmlist &parms, mf_vals &vals, arma::uword event, const double &time) {

  // if(time < 2) {
  //   Rcout << time << "; " << event << "; " << vals.patch <<  "; " << vals.action << "; " <<  std::endl << vals.rates_mat << vals.state_mat << std::endl;
  // }

  // Rcout << vals.patch << ";  " << vals.action<< std::endl;
  vals.min_patch = vals.patch*parms.n_actions;
  // if(vals.cull) {
  //   vals.rates_mat.col(vals.patch) = zeros<vec>(10);
  //   vals.cull = false;
  // } else {
    vals.rates_mat.at(0, vals.patch) =
      vals.state_mat.at(0, vals.patch) == 0 ? 0 :
    vals.state_mat.at(0, vals.patch) * (
        (parms.beta *  vals.state_mat.at(1, vals.patch) / pow(accu(vals.state_mat.submat(0, vals.patch, 2, vals.patch)), parms.rho)) + //infection events (contact)
          parms.nu * ( 1 - exp(-parms.phi * vals.state_mat.at(3, vals.patch))) + //virion uptake
          parms.lambda(vals.patch)); //external infection
    vals.rates_mat.at(1, vals.patch) = vals.state_mat.at(1, vals.patch) * parms.gamma; //recovery events
    vals.rates_mat.at(2, vals.patch) = vals.state_mat.at(0, vals.patch) * parms.mu;    //susceptible mortality
    vals.rates_mat.at(3, vals.patch) = vals.state_mat.at(1, vals.patch) * (parms.mu + parms.alpha); //infected mortality
    vals.rates_mat.at(4, vals.patch) = vals.state_mat.at(2, vals.patch) * parms.mu;     //recovered mortality
    vals.rates_mat.at(5, vals.patch) = vals.state_mat.at(1, vals.patch) * parms.sigma;   //virion shedding
    vals.rates_mat.at(6, vals.patch) = vals.state_mat.at(3, vals.patch) * parms.eta;     //viron degradation
    vals.rates_mat.submat(7, vals.patch, 9, vals.patch) = vals.state_mat.submat(0, vals.patch, 2, vals.patch) * parms.omega; // emigration of all classes
    vals.rates_mat.at(10, vals.patch) = (as_scalar(vals.state_mat.at(4, vals.patch)) == 1) ? parms.cull_rate[vals.patch] : 0.0;     //culling

    if(vals.action == 7 | vals.action == 8 | vals.action == 9) {    //If an emigration action, send the individual to another patch
      vals.rates_mat.at(0, vals.target_patch) = vals.state_mat.at(0, vals.target_patch) * (
        (parms.beta *  vals.state_mat.at(1, vals.target_patch) / pow(accu(vals.state_mat.submat(0, vals.target_patch, 2, vals.target_patch)), parms.rho)) + //infection events (contact)
          parms.nu * ( 1 - exp(-parms.phi * vals.state_mat.at(3, vals.target_patch))) + //virion uptake
          parms.lambda(vals.patch)); //external infection
      vals.rates_mat.at(1, vals.target_patch) = vals.state_mat.at(1, vals.target_patch) * parms.gamma; //recovery events
      vals.rates_mat.at(2, vals.target_patch) = vals.state_mat.at(0, vals.target_patch) * parms.mu;    //susceptible mortality
      vals.rates_mat.at(3, vals.target_patch) = vals.state_mat.at(1, vals.target_patch) * (parms.mu + parms.alpha); //infected mortality
      vals.rates_mat.at(4, vals.target_patch) = vals.state_mat.at(2, vals.target_patch) * parms.mu;     //recovered mortality
      vals.rates_mat.at(5, vals.target_patch) = vals.state_mat.at(1, vals.target_patch) * parms.sigma;   //virion shedding
      vals.rates_mat.at(6, vals.target_patch) = vals.state_mat.at(3, vals.target_patch) * parms.eta;     //viron degradation
      vals.rates_mat.submat(7, vals.target_patch, 9, vals.target_patch) = vals.state_mat.submat(0, vals.target_patch, 2, vals.target_patch) * parms.omega; // emigration of all classes
      vals.rates_mat.at(10, vals.target_patch) = (as_scalar(vals.state_mat.at(4, vals.target_patch)) == 1) ? parms.cull_rate[vals.target_patch] : 0.0;     //culling
      vals.min_patch = std::min(vals.target_patch*parms.n_actions, vals.min_patch);
    }
 // }

  //Rcout << vals.state_mat.col(vals.patch).t() << std::endl;

  if(vals.min_patch == 0) {
    cumrates[0] = rates[0];
    vals.min_patch++;
  }
  for( ; vals.min_patch < cumrates.size(); vals.min_patch++) {
    cumrates[vals.min_patch] = cumrates[vals.min_patch-1] + rates[vals.min_patch];
  }
}

bool abort_criterion(arma::vec &inits, arma::vec &state, const mf_parmlist &parms, mf_vals &vals) {
  //bool out = sum((vals.state_mat.row(1) > 0) + (1 - vals.init_I2_ind)) >= 2;
  //Rcpp::Rcout << out << std::endl;
  return false;
}

#endif
