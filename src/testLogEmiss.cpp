#include <Rcpp.h>
#include <iostream>
#include "gaston/matrix4.h"
#include "emiss.h"
#include "RVector.h"
using namespace Rcpp;

//[[Rcpp::export]]
NumericVector testLogEmiss(XPtr<matrix4> p_A, NumericVector p_, IntegerVector map_, double epsilon, int i) {
  RVector<double> p(p_);
  RVector<int> map(map_);
  matrix4 * PA(p_A);

  emiss<double> EM(PA, p, map, epsilon);
  std::vector<double> logEmiss( EM.getLogEmiss(i) );
  return wrap( logEmiss );
}
  
