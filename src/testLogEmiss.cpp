#include <Rcpp.h>
#include <iostream>
#include "gaston/matrix4.h"
#include "emiss.h"
using namespace Rcpp;

//[[Rcpp::export]]
NumericVector testLogEmiss(XPtr<matrix4> p_A, NumericVector p, IntegerVector map, double epsilon, int i) {
  emiss<double> EM(p_A, p, map, epsilon);
  std::vector<double> logEmiss( EM.getLogEmiss(i) );
  return wrap( logEmiss );
  NumericVector R(5);
  return R;
}
  
