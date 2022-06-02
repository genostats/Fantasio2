#include <Rcpp.h>
#include <RcppEigen.h>
#include "LBFGSB.h"
#include "getUserParam.h"

static LBFGSpp::LBFGSBParam<double> userParametersDouble;
static LBFGSpp::LBFGSBParam<float>  userParametersFloat;

static double lb_a = 0;
static double ub_a = INFINITY;
static double lb_f = 0;
static double ub_f = 1;

static bool _debug_ = false;

// voir la description des paramètres 
// https://lbfgspp.statr.me/doc/classLBFGSpp_1_1LBFGSBParam.html

//[[Rcpp::export]]
void setUserParam(int m, double epsilon, int past, double delta, int max_iterations, int max_submin, 
                  int max_linesearch, double min_step, double max_step, double ftol, double wolfe,
                  Rcpp::NumericVector lower, Rcpp::NumericVector upper, bool debug_) {

  userParametersDouble.m              = m;
  userParametersDouble.epsilon        = epsilon;
  userParametersDouble.past           = past;
  userParametersDouble.delta          = delta;
  userParametersDouble.max_iterations = max_iterations;
  userParametersDouble.max_submin     = max_submin;
  userParametersDouble.max_linesearch = max_linesearch;
  userParametersDouble.min_step       = min_step;
  userParametersDouble.max_step       = max_step;
  userParametersDouble.ftol           = ftol;
  userParametersDouble.wolfe          = wolfe;

  userParametersFloat.m              = m;
  userParametersFloat.epsilon        = (float) epsilon;
  userParametersFloat.past           = past;
  userParametersFloat.delta          = (float) delta;
  userParametersFloat.max_iterations = max_iterations;
  userParametersFloat.max_submin     = max_submin;
  userParametersFloat.max_linesearch = max_linesearch;
  userParametersFloat.min_step       = (float) min_step;
  userParametersFloat.max_step       = (float) max_step;
  userParametersFloat.ftol           = (float) ftol;
  userParametersFloat.wolfe          = (float) wolfe;

  lb_a = lower[0]; lb_f = lower[1];
  ub_a = upper[0]; ub_f = upper[1];
 
  _debug_ = debug_;
}

template<>
LBFGSpp::LBFGSBParam<float> getUserParam() {
  return userParametersFloat;
}

template<>
LBFGSpp::LBFGSBParam<double> getUserParam() {
  return userParametersDouble;
}

template<>
VECTOR<float> getLb() {
  VECTOR<float> lb(2);
  lb << (float) lb_a, (float) lb_f;
  return lb;
}

template<>
VECTOR<double> getLb() {
  VECTOR<double> lb(2);
  lb << lb_a, lb_f;
  return lb;
}

template<>
VECTOR<float> getUb() {
  VECTOR<float> ub(2);
  ub << (float) ub_a, (float) ub_f;
  return ub;
}

template<>
VECTOR<double> getUb() {
  VECTOR<double> ub(2);
  ub << ub_a, ub_f;
  return ub;
}

bool debug() {
  return _debug_;
}

//[[Rcpp::export]]
Rcpp::List getUserParam() {
  Rcpp::List L;
  L["m"] = userParametersDouble.m;
  L["epsilon"] = userParametersDouble.epsilon;
  L["past"] = userParametersDouble.past;
  L["delta"] = userParametersDouble.delta;
  L["max_iterations"] = userParametersDouble.max_iterations;
  L["max_submin"] = userParametersDouble.max_submin;
  L["max_linesearch"] = userParametersDouble.max_linesearch;
  L["min_step"] = userParametersDouble.min_step;
  L["max_step"] = userParametersDouble.max_step;
  L["ftol"] = userParametersDouble.ftol;
  L["wolfe"] = userParametersDouble.wolfe;
  L["lower"] = Rcpp::NumericVector::create(Rcpp::_["lb.a"] = lb_a, Rcpp::_["lb.f"] = lb_f);
  L["upper"] = Rcpp::NumericVector::create(Rcpp::_["ub.a"] = ub_a, Rcpp::_["ub.f"] = ub_f);
  L["debug"] = _debug_;
  return L;
}

