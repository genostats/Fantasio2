#include <Rcpp.h>
#include <iostream>
#include "gaston/matrix4.h"
#include <math.h>
#include "LSE.h"

using namespace Rcpp;

//[[Rcpp::export]]
NumericMatrix m4_logEmiss(XPtr<matrix4> p_A, NumericVector p, IntegerVector map, double epsilon) {
  int n = p_A->nrow; // nb snps
  int m = p_A->ncol; // nb inds
  int true_ncol = p_A->true_ncol; //
  double logeps = log(epsilon);
  int nb_snps = map.length();

  NumericMatrix logEmiss(2*m,nb_snps);
  double gg[4], hh[4]; // pour les log probas d'émission à HBD 0 et 1
  gg[3] = hh[3] = 0; // manquant -> proba 1
  for(int i = 0; i < nb_snps; i++) { // snp loop
    if(map[i] < 0 || map[i] > n) stop("Out of range SNP index");
    uint8_t * data = p_A->data[map[i]-1];
    double p_ = p[map[i]-1];
    gg[0] = 2*log(1-p_);  // p^2
    hh[0] = log( (1-epsilon)*(1-p_) + epsilon*(1-p_)*(1-p_) );
    gg[1] = log(2.) + log(1-p_) + log(p_); // 2pq
    hh[1] = logeps + gg[1];  // epsilon*2pq
    gg[2] = 2*log(p_);  // q^2
    hh[2] = log( (1-epsilon)*p_ + epsilon*p_*p_ );
    size_t k = 0;
    for(size_t j = 0; j < true_ncol; j++) {
      uint8_t x = data[j];
       for(int ss = 0; ss < 4 && (4*j + ss) < m; ss++) {
         logEmiss(2*k,i)     = gg[x&3];
         logEmiss(1+2*k++,i) = hh[x&3];
         x >>= 2;
       }
    }
  }
  return logEmiss;
}

