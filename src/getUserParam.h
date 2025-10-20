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

  /// parameters for BFGS
  LBFGSpp::LBFGSBParam<scalar_t> BFGSparam;
  VECTOR<scalar_t> lb;
  VECTOR<scalar_t> ub;
  int max_retries;
  
  // parameters for Fantasio
  int n_threads;   // in festim, pHBD and logitModel
  bool use_float;  // in festim, pHBD and logitModel
  int debug;       // for optimisation procedure only at this moment
  bool verbose;    // in R code (to be generalised)

  // parameter for fROH
  bool use_froh;
  int minNbSNPs; 
  double minROHlength; 
  double minDistHet; 
  double maxGapLength;

  userParam() : BFGSparam(), lb(2), ub(2), max_retries(5), n_threads(1), use_float(false), debug(0), verbose(false), 
                use_froh(true), minNbSNPs(400), minROHlength(2), minDistHet(1), maxGapLength(1) {
    lb << 0, 0;
    ub << INFINITY, 1;
  }


  void set(int m, double epsilon, int past, double delta, int max_iterations, int max_submin, 
           int max_linesearch, double min_step, double max_step, double ftol, double wolfe, int max_retries_,
           Rcpp::NumericVector lower, Rcpp::NumericVector upper, int n_threads_, bool use_float_, int debug_,
           bool verbose_, bool use_froh_, int minNbSNPs_, double minROHlength_, double minDistHet_, double maxGapLength_) {
  
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
    verbose = verbose_;

    use_froh = use_froh_;
    minNbSNPs = minNbSNPs_;
    minROHlength = minROHlength_;
    minDistHet = minDistHet_;
    maxGapLength = maxGapLength_;
  
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
    L["verbose"] = verbose;
    
    L["use_froh"]     = use_froh;
    L["minNbSNPs"]    = minNbSNPs;
    L["minROHlength"] = minROHlength;
    L["minDistHet"]   = minDistHet;
    L["maxGapLength"] = maxGapLength;

    return L;
  }

};

template<typename scalar_t>
userParam<scalar_t> getUserParam();

#endif
