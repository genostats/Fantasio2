#include <Rcpp.h>
#include <RcppEigen.h>
#include <iostream>
#include "gaston/matrix4.h"
#include "emiss.h"
#include "likelihoodGradient.h"

using namespace Rcpp;

//[[Rcpp::export]]
NumericVector testLikelihood(XPtr<matrix4> p_A, NumericVector p_, IntegerVector map_, NumericVector deltaDist, 
                             double epsilon, int i, double a, double f) {

  RVector<double> p(p_);
  RVector<int> map(map_);
  matrix4 * PA(p_A);

  std::vector<double> dDist;
  for(double a : deltaDist) 
    dDist.push_back(a);

  emiss<double> EM(PA, p, map, epsilon);
  likelihoodGradient<double> LG( EM.getLogEmiss(i) , dDist);
  Eigen::VectorXd grad(2);
  Eigen::VectorXd x(2);
  x << a, f;
  double lik = LG( x, grad );
  NumericVector R(3);
  R[0] = lik; R[1] = grad[0]; R[2] = grad[1];
  return R;
}

