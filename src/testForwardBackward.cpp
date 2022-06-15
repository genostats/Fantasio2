#include <Rcpp.h>
#include <RcppEigen.h>
#include <iostream>
#include "gaston/matrix4.h"
#include "emiss.h"
#include "forwardBackward.h"

using namespace Rcpp;

//[[Rcpp::export]]
NumericVector testForwardBackward(XPtr<matrix4> p_A, NumericVector p_, IntegerVector map_, NumericVector deltaDist, 
                             double epsilon, int i, double a, double f) {

  RVector<double> p(p_);
  RVector<int> map(map_);
  matrix4 * PA(p_A);

  std::vector<double> dDist;
  for(double a : deltaDist) 
    dDist.push_back(a);

  emiss<double> EM(PA, p, map, epsilon);
  return forwardBackward<double>( EM.getLogEmiss(i) , dDist, a, f);
}

