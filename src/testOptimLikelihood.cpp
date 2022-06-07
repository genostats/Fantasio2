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
NumericVector testOptimLikelihood(XPtr<matrix4> p_A, NumericVector p_, IntegerVector map_, NumericVector deltaDist, double epsilon, int i) {
  RVector<double> p(p_);
  RVector<int> map(map_);
  matrix4 * PA(p_A);

  std::vector<double> dDist;
  for(double a : deltaDist) 
    dDist.push_back(a);

  emiss<double> EM(PA, p, map, epsilon);
  likelihoodGradient<double> LG( EM.getLogEmiss(i) , dDist, -1);

  userParam<double> pars = getUserParam<double>();

  LBFGSpp::LBFGSBParam<double> param = pars.BFGSparam;
  LBFGSpp::LBFGSBSolver<double> solver(param);

  VECTOR<double> lb = pars.lb;
  VECTOR<double> ub = pars.ub;

  VECTOR<double> x(2);
  x << 0.05, 0.05;
  double fx;

  int niter = 0;
  int count = 1;
  while(true) {
    try {
      niter += solver.minimize(LG, x, fx, lb, ub);
      break;
    } catch(std::runtime_error const & e) {
      niter++;
      if(count++ > pars.max_retries) {
        warning(e.what());
        break;
      }
    }
  }

  NumericVector R = NumericVector::create( _["n.iter"] = niter, _["a"] = x[0], _["f"] = x[1], _["fx"] = fx);
  return R;
}
