#ifndef STRUCTURES  //compliler flags prevent this from being re-run if included
#define STRUCTURES  // in multiple files

#include <RcppArmadillo.h>
#include <vector>
#include <queue>

using namespace Rcpp;
using namespace arma;

// A structure to hold a single list of parameters
// [[Rcpp::depends(RcppArmadillo)]]
struct mf_parmlist {
  double beta;              //contact rate for direct transmission
  double gamma;             //recovery rate
  double mu;                //base mortality rate
  double alpha;             //disease mortality rate
  double phi;               //infectiousness of environmental virions
  double eta;               //degradation rate of environmental virions
  double nu;                //uptake rate of environmental virion
  double sigma;             //virion shedding rate
  double omega;             //movement rate
  double rho;               //contact rate nonlinearity 0=density-dependent, 1=frequency-dependent
  arma::rowvec lambda;           //force of infection from external sources
  arma::mat chi;            //patch connectivity matrix
  arma::mat chi_cum;        //patch connectivity matrix, row-wise cumulative
  int K;                    //total number of patches
  arma::rowvec tau_crit;          //passive surveillance period
  arma::irowvec I_crit;               //number of deaths that would trigger reporting
  arma::rowvec pi_report;         //probability of reporting deaths
  arma::rowvec pi_detect;         //probability of detecting AI
  arma::rowvec cull_rate;      //rate at which culling occurs at a reported farm
  bool abort;

  int n_actions;            // number of possible actions
  int n_pstates;            // number of possible classes
  arma::imat action_matrix; // action-class matrix

  mf_parmlist() ;  //default constructor

  mf_parmlist(List parmlist) {  //Creates this structure from an R list
    beta = as<double>(parmlist["beta"]);
    gamma = as<double>(parmlist["gamma"]);
    mu = as<double>(parmlist["mu"]);
    alpha = as<double>(parmlist["alpha"]);
    phi = as<double>(parmlist["phi"]);
    eta = as<double>(parmlist["eta"]);
    nu = as<double>(parmlist["nu"]);
    sigma = as<double>(parmlist["sigma"]);
    omega = as<double>(parmlist["omega"]);
    rho = as<double>(parmlist["rho"]);
    lambda = as<arma::rowvec>(parmlist["lambda"]);
    chi = as<arma::mat>(parmlist["chi"]);
    chi.each_col() /= sum(chi, 1);
    chi_cum = cumsum(chi, 1);
    K = chi.n_rows;
    tau_crit = as<arma::rowvec>(parmlist["tau_crit"]);
    I_crit = as<arma::irowvec>(parmlist["I_crit"]);
    pi_report = as<arma::rowvec>(parmlist["pi_report"]);
    pi_detect = as<arma::rowvec>(parmlist["pi_detect"]);
    cull_rate = 1 / as<arma::rowvec>(parmlist["cull_time"]);
    abort = as<bool>(parmlist["abort"]);
    n_actions = 11;
    n_pstates = 5;
    action_matrix = {
      {-1,  1,  0,  0,  0},  //0 infection
      { 0, -1, +1,  0,  0},  //1 recovery
      {-1,  0,  0,  0,  0},  //2 S mortality
      { 0, -1,  0,  0,  0},  //3 I mortality
      { 0,  0, -1,  0,  0},  //4 R mortality
      { 0,  0,  0, +1,  0},  //5 shedding
      { 0,  0,  0, -1,  0},  //6 degradation
      {-1,  0,  0,  0,  0},  //7 S emigration
      { 0, -1,  0,  0,  0},  //8 I emigration
      { 0,  0, -1,  0,  0},  //9 R emigration
      { 0,  0,  0,  0, +1},  //10 culling
    };
    action_matrix = trans(action_matrix);

  }
};

struct mf_vals {  //holds interim objects

  arma::mat rates_mat;
  arma::mat state_mat;
  arma::urowvec init_I_ind;
  arma::urowvec init_I2_ind;
  uword target_patch;
  uword patch;
  uword action;
  uword min_patch;
  std::vector<std::priority_queue<double, std::vector<double>, std::greater<double>>> deathtimes;
  bool cull;

  mf_vals() {
    patch = 0;
    target_patch = 0;
    min_patch = 0;
  };
};

#endif
