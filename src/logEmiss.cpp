#include <Rcpp.h>
#include <iostream>
#include <ctime>
#include <math.h>
#include "LSE.h"

// --------------- this should disappear, kept for comparison purposes -----------------

using namespace Rcpp;

// epsilon = proba d'erreur de génotypage
// logFreq = log des fréquences : (max allele) x N
// output : matrice 2 x N des probas d'emission (ligne i : conditionnellement à S = i)
// [[Rcpp::export]]
NumericMatrix logEmiss(NumericVector Y1, NumericVector Y2, NumericMatrix logFreq, double epsilon) {
  int N = Y1.length();
  double logeps = log(epsilon);
  double logumeps = log1p(-epsilon);
  NumericMatrix logEmiss(2,N);
  for(int n = 0; n < N; n++) {
    if(Y1[n] == 0 && Y2[n] == 0) { // manquant -> proba 1
      logEmiss(0, n) = logEmiss(1, n) = 0;
      continue;
    } else {
      if(Y1[n] == Y2[n]) { // homozygote
        logEmiss(0,n) = 2*logFreq(n,Y1[n]-1);  // p^2
        logEmiss(1,n) = LSE(logumeps+logFreq(n,Y1[n]-1),logeps+2*logFreq(n,Y1[n]-1)); // (1-epsilon)*p + epsilon*p^2 [?]
      } else { // hétérozygote
        logEmiss(0,n) = log(2.) + logFreq(n,Y1[n]-1) + logFreq(n,Y2[n]-1); // 2pq
        logEmiss(1,n) = logeps + log(2.) + logFreq(n,Y1[n]-1) + logFreq(n,Y2[n]-1); // epsilon*2pq [?]
      }
    }
  }
  return logEmiss;
}

