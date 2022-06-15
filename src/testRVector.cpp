#include <Rcpp.h>
#include <iostream>
#include "RVector.h"

using namespace Rcpp;

//[[Rcpp::export]]
void testRVector(NumericVector x, IntegerVector y, LogicalVector z) {
  RVector<double> a(x);
  RVector<int> b(y);
  RVector<int> c(z);

  for(auto & e : a) std::cout << e << "\n";
  for(auto & e : b) std::cout << e << "\n";
  for(auto & e : c) std::cout << e << "\n";

}


