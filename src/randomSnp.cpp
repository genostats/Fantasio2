#include <Rcpp.h>
#include <iostream>
#include <random>

using namespace Rcpp;

double mt_runif();

// [[Rcpp::export]]
IntegerVector randomSnp(List L) {
  IntegerVector beg = L["beg"] ;
  IntegerVector end = L["end"] ;
  int n = beg.length();
  IntegerVector submap = no_init(n);
  for(int i = 0; i < n; i++) {
    submap[i] = beg[i] + std::floor( mt_runif() * (end[i] - beg[i] + 1) );
  }
  return submap;
}
