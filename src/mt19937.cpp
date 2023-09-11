#include <Rcpp.h>
#include <iostream>
#include <random>

using namespace Rcpp;
//[[Rcpp::plugins(cpp11)]]

static std::mt19937 mt;

// [[Rcpp::export]]
IntegerVector getSeed() {
  std::stringstream sx;
  sx << mt;
  uint32_t n;
  std::vector<int> seed;
  while(sx >> n) seed.push_back((int) n);
  return wrap(seed);
}

// [[Rcpp::export]]
void setSeed(IntegerVector seed) {
  std::stringstream sx;
  for(int n : seed) sx << (uint32_t) n << " "; 
  sx >> mt;
}

// [[Rcpp::export]]
double mt_runif() {
  return mt() / 4294967296.;
}
