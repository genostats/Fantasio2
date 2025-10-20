#include <Rcpp.h>
#include <iostream>
using namespace Rcpp;

// [[Rcpp::export]]
double acfd(double d, NumericVector z, IntegerVector Chr, NumericVector Dist) {
  unsigned int n = z.size();
  double S = 0;
  unsigned int i = 0, j = 1, k = 0;
  int chri = Chr[i];
  while(i < n && j < n - 1) {
    bool endOfChr = false;
    double dij = Dist[j] - Dist[i];
    while(j < n - 1) {
      if(Chr[j+1] != chri) {
        endOfChr = true;
        break;
      }
      double dij1 = Dist[j+1] - Dist[i];
      if( std::abs(dij1 - d) > std::abs(dij - d) ) break;
      j++;
      dij = dij1;
    }
    if(endOfChr) {
      i = j+1;
      j = i+1;
      chri = Chr[i];
    } else {
      S = S + z[i]*z[j];
      k++;
      i++;
    }
  }
  return S / (double) k;
}
