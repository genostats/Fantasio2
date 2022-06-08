#include <Rcpp.h>
#include <RcppEigen.h>
#include "LBFGSB.h"

#ifndef userParams
#define userParams

template<typename scalar_t>
using VECTOR = Eigen::Matrix<scalar_t, Eigen::Dynamic, 1>;

template<typename scalar_t>
class userParam {
public:
  LBFGSpp::LBFGSBParam<scalar_t> BFGSparam;
  VECTOR<scalar_t> lb;
  VECTOR<scalar_t> ub;
  int max_retries;
  
  int n_threads;
  bool use_float;
  int debug;

  userParam() : BFGSparam(), lb(2), ub(2), max_retries(5), n_threads(1), use_float(false), debug(0) {
    lb << 0, 0;
    ub << INFINITY, 1;
  }


  void set(int m, double epsilon, int past, double delta, int max_iterations, int max_submin, 
           int max_linesearch, double min_step, double max_step, double ftol, double wolfe, int max_retries_,
           Rcpp::NumericVector lower, Rcpp::NumericVector upper, int n_threads_, bool use_float_, int debug_) {
  
    BFGSparam.m              = m;
    BFGSparam.epsilon        = epsilon;
    BFGSparam.past           = past;
    BFGSparam.delta          = delta;
    BFGSparam.max_iterations = max_iterations;
    BFGSparam.max_submin     = max_submin;
    BFGSparam.max_linesearch = max_linesearch;
    BFGSparam.min_step       = min_step;
    BFGSparam.max_step       = max_step;
    BFGSparam.ftol           = ftol;
    BFGSparam.wolfe          = wolfe;
  
    n_threads = n_threads_;
    max_retries = max_retries_;
    use_float = use_float_;
    debug = debug_;
  
    lb[0] = lower[0];
    lb[1] = lower[1];
    ub[0] = upper[0];
    ub[1] = upper[1];
  }
 
  Rcpp::List toList() {
    Rcpp::List L;
    L["m"] = BFGSparam.m;
    L["epsilon"] = BFGSparam.epsilon;
    L["past"] = BFGSparam.past;
    L["delta"] = BFGSparam.delta;
    L["max_iterations"] = BFGSparam.max_iterations;
    L["max_submin"] = BFGSparam.max_submin;
    L["max_linesearch"] = BFGSparam.max_linesearch;
    L["min_step"] = BFGSparam.min_step;
    L["max_step"] = BFGSparam.max_step;
    L["ftol"] = BFGSparam.ftol;
    L["wolfe"] = BFGSparam.wolfe;
    L["max_retries"] = max_retries;
  
    L["lower"] = Rcpp::NumericVector::create(Rcpp::_["lb.a"] = lb[0], Rcpp::_["lb.f"] = lb[1]);
    L["upper"] = Rcpp::NumericVector::create(Rcpp::_["ub.a"] = ub[0], Rcpp::_["ub.f"] = ub[1]);

    L["n_threads"] = n_threads;
    L["use_float"] = use_float;
    L["debug"] = debug;
    return L;
  }

};

template<typename scalar_t>
userParam<scalar_t> getUserParam();

#endif
