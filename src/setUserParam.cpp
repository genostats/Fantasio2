#include <Rcpp.h>
#include <RcppEigen.h>
#include "LBFGSB.h"
#include "getUserParam.h"

static userParam<double> userParamDouble;
static userParam<float>  userParamFloat;

// voir la description des paramètres 
// https://lbfgspp.statr.me/doc/classLBFGSpp_1_1LBFGSBParam.html

//[[Rcpp::export]]
void setUserParam(int m, double epsilon, int past, double delta, int max_iterations, int max_submin, 
                  int max_linesearch, double min_step, double max_step, double ftol, double wolfe, int max_retries,
                  Rcpp::NumericVector lower, Rcpp::NumericVector upper, int n_threads, bool use_float, int debug, 
                  bool use_froh, int minNbSNPs, double minROHlength, double minDistHet, double maxGapLength) {

  userParamDouble.set(m, epsilon, past, delta, max_iterations, max_submin, max_linesearch, min_step, max_step,
                      ftol, wolfe, max_retries, lower, upper, n_threads, use_float, debug,
                      use_froh, minNbSNPs, minROHlength, minDistHet, maxGapLength);

  userParamFloat.set(m, epsilon, past, delta, max_iterations, max_submin, max_linesearch, min_step, max_step,
                     ftol, wolfe, max_retries, lower, upper, n_threads, use_float, debug,
                      use_froh, minNbSNPs, minROHlength, minDistHet, maxGapLength);
}

template<>
userParam<float> getUserParam() {
  return userParamFloat;
}

template<>
userParam<double> getUserParam() {
  return userParamDouble;
}

//[[Rcpp::export]]
Rcpp::List getUserParam() {
  return userParamDouble.toList();
}

