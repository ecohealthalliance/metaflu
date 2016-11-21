#ifndef STRUCTURES  //compliler flags prevent this from being re-run if included
#define STRUCTURES  // in multiple files

#include <RcppArmadillo.h>

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
  double lambda;           //force of infection from external sources
  arma::mat chi;            //patch connectivity matrix
  arma::mat chi_cum;        //patch connectivity matrix, row-wise cumulative
  int K;                    //total number of patches

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
    lambda = as<double>(parmlist["lambda"]);
    chi = as<arma::mat>(parmlist["chi"]);
    chi_cum = cumsum(chi, 1);
    K = chi.n_rows;

    n_actions = 10;
    n_pstates = 4;
    action_matrix = {
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
    action_matrix = trans(action_matrix);

  }
};

struct mf_vals {  //holds interim objects

  arma::mat rates_mat;
  arma::mat state_mat;
  uword target_patch;
  uword patch;
  uword action;
  mf_vals() {};

};

#endif
