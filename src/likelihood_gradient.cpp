// [[Rcpp::depends(RcppParallel)]]
#include <Rcpp.h>
#include <iostream>
#include <ctime>
#include <math.h>
#include "LSE.h"
#ifndef SHOW
#define SHOW(x) Rcpp::Rcout << #x << " = " << (x) << std::endl;
#endif

// --------------- this should disappear, kept for comparison purposes -----------------

using namespace Rcpp;

// pour 0 < f < 1
NumericVector logLikelihood_gradientf(NumericMatrix logEmiss, NumericVector Dist, double a, double f) {
  double lt00, lt01, lt10, lt11; // log proba transition
  double df_lt00, df_lt01, df_lt10, df_lt11;
  double da_lt00, da_lt01, da_lt10, da_lt11;
  int N = logEmiss.ncol();
  // état initial
  double alpha0 = log1p(-f);
  double alpha1 = log(f);

  double da_alpha0 = 0, da_alpha1 = 0;
  double df_alpha0 = -1/(1-f);
  double df_alpha1 = 1/f;  

  // ceux là ne changent pas de valeur !!  
  df_lt01 = 1/f;
  df_lt10 = -1/(1-f);

  double da_alpha0_, df_alpha0_;

  for(int n = 1; n < N; n++) {
    // calcul des lt -------------------------
    double d = Dist[n-1];
    if(d < 0) { // changement de chromosomoe
      lt00 = log(1-f);
      lt01 = log(f);
      lt10 = log(1-f);
      lt11 = log(f);
      df_lt00 = -1/(1-f);
      df_lt11 = 1/f;
      da_lt00 = da_lt01 = da_lt10 = da_lt11 = 0.;
    } else {
      double ex = expm1(-a*d);
      double t00 = 1 + f*ex;
      double t01 = -f*ex;
      double t10 = -(1-f)*ex;
      double t11 = 1 + (1-f)*ex;

      lt00 = log(t00);
      lt01 = log(t01);
      lt10 = log(t10);
      lt11 = log(t11);

      df_lt00 = ex/t00; 
      df_lt11 = -ex/t11; 

      da_lt00 = -d*f*(ex+1)/t00;
      // gérer le cas d = 0 
      // (qui ne devrait normalement pas être rencontré avec des données parfaitement propres...)
     if(d > 0) 
        da_lt01 = da_lt10 = d*(ex+1)/(-ex);
      else 
        da_lt01 = da_lt10 = 1/a;
      da_lt11 = -d*(1-f)*(ex+1)/t11;
    }
    // ---------------------------

    double u0 = lt00 + logEmiss(0,n-1) + alpha0;
    double v0 = lt10 + logEmiss(1,n-1) + alpha1;
    double u1 = lt01 + logEmiss(0,n-1) + alpha0;
    double v1 = lt11 + logEmiss(1,n-1) + alpha1;

    alpha0 = LSE( u0, v0 );
    alpha1 = LSE( u1, v1 );

    df_alpha0_ = (df_lt00 + df_alpha0)/(1+exp(v0-u0)) + (df_lt10 + df_alpha1)/(1 + 1/exp(v0-u0));    
    df_alpha1  = (df_lt01 + df_alpha0)/(1+exp(v1-u1)) + (df_lt11 + df_alpha1)/(1 + 1/exp(v1-u1));
    df_alpha0 = df_alpha0_; 
   
    if(a == 0 && d >= 0) {
      da_alpha0_ = da_alpha0 + Dist[n-1]*f*expm1(alpha1-alpha0);
      da_alpha1  = da_alpha1 + Dist[n-1]*(1-f)*expm1(alpha0-alpha1);
      da_alpha0 = da_alpha0_;
    } else {
      da_alpha0_ = (da_lt00 + da_alpha0)/(1+exp(v0-u0)) + (da_lt10 + da_alpha1)/(1 + 1/exp(v0-u0));
      da_alpha1  = (da_lt01 + da_alpha0)/(1+exp(v1-u1)) + (da_lt11 + da_alpha1)/(1 + 1/exp(v1-u1));
      da_alpha0 = da_alpha0_;
    }
  }
  double u = alpha0 + logEmiss(0,N-1);
  double v = alpha1 + logEmiss(1,N-1);

  double lik = LSE(u, v);
  double df_lik = df_alpha0/(1+exp(v-u)) + df_alpha1/(1+1/exp(v-u));
  double da_lik = da_alpha0/(1+exp(v-u)) + da_alpha1/(1+1/exp(v-u));

  return NumericVector::create(lik, da_lik, df_lik);
}

// f == 0
NumericVector logLikelihood_gradientf0(NumericMatrix logEmiss, NumericVector Dist, double a) {
  int N = logEmiss.ncol();
  // état initial
  double alpha0 = 0;
  double df_alpha0 = -1;
  double df_alpha1 = 1;     // ici, df_alpha1 est df a1 / a0 et non df a1 / a1
  
  for(int n = 1; n < N; n++) {
    double d = Dist[n-1];
    alpha0 += logEmiss(0,n-1);

    if(d >=0 ) {
      df_alpha0 += expm1(-a*d)*( 1 - df_alpha1*exp( logEmiss(1,n-1) - logEmiss(0,n-1) ) );
      df_alpha1 = -expm1(-a*d) + df_alpha1 * exp(logEmiss(1,n-1) - logEmiss(0,n-1) - a*d );
    } else { // correspond à d = +inf
      df_alpha0 += -( 1 - df_alpha1*exp( logEmiss(1,n-1) - logEmiss(0,n-1) ) );
      df_alpha1 = 1;
    }
  }

  double lik = alpha0 + logEmiss(0,N-1);
  double df_lik = df_alpha0 + df_alpha1 * exp(logEmiss(1,N-1)-logEmiss(0,N-1));
  double da_lik = 0;

  return NumericVector::create(lik, da_lik, df_lik);
}


// f == 1
NumericVector logLikelihood_gradientf1(NumericMatrix logEmiss, NumericVector Dist, double a) {
  int N = logEmiss.ncol();
  // état initial
  double alpha1 = 0;
  double df_alpha0 = -1;   // ici, df_alpha0 est df a0 / a1 et non df a0 / a0
  double df_alpha1 = 1;  
  
  for(int n = 1; n < N; n++) {
    double d = Dist[n-1];
    alpha1 += logEmiss(1, n-1);

    if(d >= 0) {
      df_alpha1 += -expm1(-a*d)*(1 + df_alpha0*exp(logEmiss(0,n-1) - logEmiss(1,n-1)));
      df_alpha0 = expm1(-a*d) + df_alpha0 * exp(logEmiss(0,n-1) - logEmiss(1, n-1) - a*d);
    } else {
      df_alpha1 += (1 + df_alpha0*exp(logEmiss(0,n-1) - logEmiss(1,n-1)));
      df_alpha0 = -1;
    }
  }

  double lik = alpha1 + logEmiss(1, N-1);
  double df_lik = df_alpha1 + df_alpha0 * exp(logEmiss(0,N-1)-logEmiss(1,N-1));
  double da_lik = 0;

  return NumericVector::create(lik, da_lik, df_lik);
}



// [[Rcpp::export]]
NumericVector logLikelihood_gradient(NumericMatrix logEmiss, NumericVector Dist, double a, double f) {
  if(f == 0)
    return logLikelihood_gradientf0(logEmiss, Dist, a);
  else if(f == 1)
    return logLikelihood_gradientf1(logEmiss, Dist, a);
  else
    return logLikelihood_gradientf(logEmiss, Dist, a, f);
}

