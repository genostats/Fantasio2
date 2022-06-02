#include <Rcpp.h>
#include <RcppEigen.h>
#include <iostream>
#include "gaston/matrix4.h"
#include "emiss.h"
#include "likelihoodGradient.h"

using namespace Rcpp;

//[[Rcpp::export]]
NumericVector testLikelihood(XPtr<matrix4> p_A, NumericVector p, IntegerVector map, NumericVector deltaDist, 
                             double epsilon, int i, double a, double f) {
  std::vector<double> dDist;
  for(double a : deltaDist) 
    dDist.push_back(a);

  emiss<double> EM(p_A, p, map, epsilon);
  likelihoodGradient<double> LG( EM.getLogEmiss(i) , dDist);
  Eigen::VectorXd grad(2);
  Eigen::VectorXd x(2);
  x << a, f;
  double lik = LG( x, grad );
  NumericVector R(3);
  R[0] = lik; R[1] = grad[0]; R[2] = grad[1];
  return R;
}

