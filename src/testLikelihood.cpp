#include <Rcpp.h>
#include <iostream>
#include "gaston/matrix4.h"
#include "emiss.h"
#include "likelihoodGradient.h"
#include "LBFGSB.h"

using namespace Rcpp;

//[[Rcpp::export]]
NumericVector testLikelihood(XPtr<matrix4> p_A, NumericVector p, IntegerVector map, NumericVector deltaDist, double epsilon, int i, double a, double f) {
  emiss<double> EM(p_A, p, map, epsilon);
  likelihoodGradient<double> LG( EM.getLogEmiss(i) , deltaDist);
  Eigen::VectorXd grad(2);
  Eigen::VectorXd x(2);
  x << a, f;
  double lik = LG( x, grad );
  NumericVector R(3);
  R[0] = lik; R[1] = grad[0]; R[2] = grad[1];
  return R;
}

//[[Rcpp::export]]
NumericVector testOptimLikelihood(XPtr<matrix4> p_A, NumericVector p, IntegerVector map, NumericVector deltaDist, double epsilon, int i) {
  emiss<double> EM(p_A, p, map, epsilon);
  likelihoodGradient<double> LG( EM.getLogEmiss(i) , deltaDist, -1);

  LBFGSpp::LBFGSBParam<double> param;
  LBFGSpp::LBFGSBSolver<double> solver(param);

  Eigen::VectorXd lb(2); 
  lb << 0, 0;

  Eigen::VectorXd ub(2); 
  ub << std::numeric_limits<double>::infinity(), 1;

  Eigen::VectorXd x(2);
  x << 0.1, 0.1;
  double fx;
  int niter = solver.minimize(LG, x, fx, lb, ub);

  NumericVector R(3);
  R[0] = niter; R[1] = x[0]; R[2] = x[1];
  return R;
}
