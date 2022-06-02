#include <Rcpp.h>
#include <RcppEigen.h>
#include <iostream>
#include "gaston/matrix4.h"
#include "emiss.h"
#include "likelihoodGradient.h"
#include "LBFGSB.h"
#include "getUserParam.h"

using namespace Rcpp;

//[[Rcpp::export]]
NumericVector testOptimLikelihood(XPtr<matrix4> p_A, NumericVector p, IntegerVector map, NumericVector deltaDist, double epsilon, int i) {
  std::vector<double> dDist;
  for(double a : deltaDist) 
    dDist.push_back(a);

  emiss<double> EM(p_A, p, map, epsilon);
  likelihoodGradient<double> LG( EM.getLogEmiss(i) , dDist, -1);

  LBFGSpp::LBFGSBParam<double> param = getUserParam<double>();
  LBFGSpp::LBFGSBSolver<double> solver(param);

  VECTOR<double> lb = getLb<double>();
  VECTOR<double> ub = getUb<double>();

  VECTOR<double> x(2);
  x << 0.05, 0.05;
  double fx;

  int niter = 0;
  int max_tries = 30;
  int count = 1;
  while(true) {
    try {
      niter += solver.minimize(LG, x, fx, lb, ub);
      break;
    } catch(std::runtime_error const & e) {
      warning(e.what());
      niter++;
      if(count++ > max_tries) 
        break;
    }
  }

  NumericVector R = NumericVector::create( _["n.iter"] = niter, _["a"] = x[0], _["f"] = x[1], _["fx"] = fx);
  return R;
}
