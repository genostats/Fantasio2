// [[Rcpp::depends(RcppParallel)]]
#include <Rcpp.h>
#include <iostream>
#include <ctime>
#include <math.h>
#include "LSE.h"
#define SHOW(x) Rcpp::Rcout << #x << " = " << (x) << std::endl;

using namespace Rcpp;

inline void logTrans4(double d, double a, double f, double lf, double lumf, double & lt00, double & lt01, double & lt10, double & lt11) {
  if(d < 0) { // changement de chromosomoe
    lt00 = lumf;
    lt01 = lf;
    lt10 = lumf;
    lt11 = lf;
    return;
  }
  // attention à la façon de calculer log(1 - exp(-a*d))
  // quand a petit, log1p(-exp(-a*d)) marche moins bien que log(-expm1(-a*d))
  // [enfin ça change pas grand chose]
  lt00 = log1p(f*expm1(-a*d));
  lt01 = lf + log(-expm1(-a*d));
  lt10 = lumf + log(-expm1(-a*d));
  lt11 = log1p((1-f)*expm1(-a*d));
  return;
}

/**** algorithme forward - backward pour le calcul des probas a posteriori ****/

// [[Rcpp::export]]
NumericMatrix forward_backward(NumericMatrix logEmiss, NumericVector Dist, double a, double f) {
  double logf   = log(f);
  double logumf = log(1-f);
  double lt00, lt01, lt10, lt11; // log proba transition
  int N = logEmiss.ncol();
  // stocker les alpha
  NumericMatrix Alpha(2,N);
  // état initial
  double alpha0 = logumf;
  double alpha1 = logf;
  double alpha0_;
  Alpha(0,0) = alpha0;
  Alpha(1,0) = alpha1;
  // forward
  for(int n = 1; n < N; n++) {
    // calcul des lt
    logTrans4(Dist[n-1], a, f, logf, logumf, lt00, lt01, lt10, lt11);
    alpha0_ = LSE( lt00 + logEmiss(0,n-1) + alpha0, lt10 + logEmiss(1,n-1) + alpha1 );
    alpha1  = LSE( lt01 + logEmiss(0,n-1) + alpha0, lt11 + logEmiss(1,n-1) + alpha1 );
    alpha0  = alpha0_;
    Alpha(0,n) = alpha0;
    Alpha(1,n) = alpha1;
  }
  // backward
  // [on calcule les deux probas, car sinon je crains qu'on puisse avoir un souci d'arrondi quand on a une proba proche de 1]
  NumericMatrix Beta(2,N);
  double beta0 = 1/(1 + exp( alpha1 + logEmiss(1,N-1) - alpha0 - logEmiss(0,N-1)));
  double beta1 = 1/(1 + exp( alpha0 + logEmiss(0,N-1) - alpha1 - logEmiss(1,N-1)));
  double beta0_;
  Beta(0,N-1) = beta0; Beta(1,N-1) = beta1;
  for(int n = N-2; n >=0; n--) {
    logTrans4(Dist[n], a, f, logf, logumf, lt00, lt01, lt10, lt11);
    beta0_ = exp(logEmiss(0,n))*( beta0 * exp(lt00 + Alpha(0,n) - Alpha(0,n+1)) + beta1 * exp(lt01 + Alpha(0,n) - Alpha(1,n+1)) );
    beta1  = exp(logEmiss(1,n))*( beta0 * exp(lt10 + Alpha(1,n) - Alpha(0,n+1)) + beta1 * exp(lt11 + Alpha(1,n) - Alpha(1,n+1)) );
    beta0 = beta0_;
    Beta(0,n) = beta0; Beta(1,n) = beta1;
  }
  return Beta;
}



