#include <RcppEigen.h>
#include <Rcpp.h>
#include <iostream>
#include <ctime>
#include <math.h>
#include "LSE.h"
#define SHOW(x) Rcpp::Rcout << #x << " = " << (x) << std::endl;

using namespace Rcpp;

// pour 0 < f < 1
template<typename scalar_t>
class likelihoodGradient {
private:
  std::vector<scalar_t> & logEmiss;
  NumericVector deltaDist;
  scalar_t scale;
  
public:
  likelihoodGradient(std::vector<scalar_t> & logEmiss_, NumericVector deltaDist_) : logEmiss(logEmiss_), deltaDist(deltaDist_), scale(1) {}
  likelihoodGradient(std::vector<scalar_t> & logEmiss_, NumericVector deltaDist_, scalar_t scale_) : logEmiss(logEmiss_), deltaDist(deltaDist_), scale(scale_) {}
  
  scalar_t operator()(const Eigen::VectorXd & x, Eigen::VectorXd & grad) {
    scalar_t a = x[0];
    scalar_t f = x[1];
    if(f == 0)
      return f0(a, grad);
    else if(f == 1)
      return f1(a, grad);
    else
      return ff(a, f, grad);
  }
 
private:
 
  scalar_t ff(scalar_t a, scalar_t f, Eigen::VectorXd & grad) {
SHOW(a)
SHOW(f)
    scalar_t lt00, lt01, lt10, lt11; // log proba transition
    scalar_t df_lt00, df_lt01, df_lt10, df_lt11;
    scalar_t da_lt00, da_lt01, da_lt10, da_lt11;

    // int N = logEmiss.ncol();
    int N = deltaDist.size() + 1;
    // état initial
    scalar_t alpha0 = log1p(-f);
    scalar_t alpha1 = log(f);

    scalar_t da_alpha0 = 0, da_alpha1 = 0;
    scalar_t df_alpha0 = -1/(1-f);
    scalar_t df_alpha1 = 1/f;  

    // ceux là ne changent pas de valeur !!  
    df_lt01 = 1/f;
    df_lt10 = -1/(1-f);

    scalar_t da_alpha0_, df_alpha0_;

    for(int n = 0; n < N-1; n++) {
      // calcul des lt -------------------------
      scalar_t d = deltaDist[n];
      if(d < 0) { // changement de chromosomoe
        lt00 = log(1-f);
        lt01 = log(f);
        lt10 = log(1-f);
        lt11 = log(f);
        df_lt00 = -1/(1-f);
        df_lt11 = 1/f;
        da_lt00 = da_lt01 = da_lt10 = da_lt11 = 0.;
      } else {
        scalar_t ex = expm1(-a*d);
        scalar_t t00 = 1 + f*ex;
        scalar_t t01 = -f*ex;
        scalar_t t10 = -(1-f)*ex;
        scalar_t t11 = 1 + (1-f)*ex;
  
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
      // -------- fin calcul des lt ---------------

      scalar_t u0 = lt00 + logEmiss[2*n]   + alpha0;
      scalar_t v0 = lt10 + logEmiss[2*n+1] + alpha1;
      scalar_t u1 = lt01 + logEmiss[2*n]   + alpha0;
      scalar_t v1 = lt11 + logEmiss[2*n+1] + alpha1;

      alpha0 = LSE( u0, v0 );
      alpha1 = LSE( u1, v1 );

      df_alpha0_ = (df_lt00 + df_alpha0)/(1+exp(v0-u0)) + (df_lt10 + df_alpha1)/(1 + 1/exp(v0-u0));    
      df_alpha1  = (df_lt01 + df_alpha0)/(1+exp(v1-u1)) + (df_lt11 + df_alpha1)/(1 + 1/exp(v1-u1));
      df_alpha0 = df_alpha0_; 
   
      if(a == 0 && d >= 0) {
        da_alpha0_ = da_alpha0 + deltaDist[n-1]*f*expm1(alpha1-alpha0);
        da_alpha1  = da_alpha1 + deltaDist[n-1]*(1-f)*expm1(alpha0-alpha1);
        da_alpha0 = da_alpha0_;
      } else {
        da_alpha0_ = (da_lt00 + da_alpha0)/(1+exp(v0-u0)) + (da_lt10 + da_alpha1)/(1 + 1/exp(v0-u0));
        da_alpha1  = (da_lt01 + da_alpha0)/(1+exp(v1-u1)) + (da_lt11 + da_alpha1)/(1 + 1/exp(v1-u1));
        da_alpha0 = da_alpha0_;
      }
    }
    scalar_t u = alpha0 + logEmiss[2*N-2];
    scalar_t v = alpha1 + logEmiss[2*N-1];

    scalar_t lik = LSE(u, v);
    scalar_t df_lik = df_alpha0/(1+exp(v-u)) + df_alpha1/(1+1/exp(v-u));
    scalar_t da_lik = da_alpha0/(1+exp(v-u)) + da_alpha1/(1+1/exp(v-u));

    grad[0] = scale*da_lik;
    grad[1] = scale*df_lik;
    return scale*lik;
  }

  // quand f == 0
  scalar_t f0(scalar_t a, Eigen::VectorXd & grad) {
    int N = deltaDist.size() + 1;
    // état initial
    scalar_t alpha0 = 0;
    scalar_t df_alpha0 = -1;
    scalar_t df_alpha1 = 1;     // ici, df_alpha1 est df a1 / a0 et non df a1 / a1
  
    for(int n = 0; n < N-1; n++) {
      scalar_t d = deltaDist[n];
      alpha0 += logEmiss[2*n];

      if(d >=0 ) {
        df_alpha0 += expm1(-a*d)*( 1 - df_alpha1*exp( logEmiss[2*n+1] - logEmiss[2*n] ) );
        df_alpha1 = -expm1(-a*d) + df_alpha1 * exp(logEmiss[2*n+1] - logEmiss[2*n] - a*d );
      } else { // correspond à d = +inf
        df_alpha0 += -( 1 - df_alpha1*exp( logEmiss[2*n+1] - logEmiss[2*n] ) );
        df_alpha1 = 1;
      }
    }

    scalar_t lik = alpha0 + logEmiss[2*N-2];
    scalar_t df_lik = df_alpha0 + df_alpha1 * exp(logEmiss[2*N-1]-logEmiss[2*N-2]);
    scalar_t da_lik = 0;

    grad[0] = scale*da_lik;
    grad[1] = scale*df_lik;
    return scale*lik;
  }


  // f == 1
  scalar_t f1(scalar_t a, Eigen::VectorXd & grad) {
    int N = deltaDist.size() + 1;
    // état initial
    scalar_t alpha1 = 0;
    scalar_t df_alpha0 = -1;   // ici, df_alpha0 est df a0 / a1 et non df a0 / a0
    scalar_t df_alpha1 = 1;  
  
    for(int n = 0; n < N-1; n++) {
      scalar_t d = deltaDist[n];
      alpha1 += logEmiss[2*n+1];

      if(d >= 0) {
        df_alpha1 += -expm1(-a*d)*(1 + df_alpha0*exp(logEmiss[2*n] - logEmiss[2*n+1]));
        df_alpha0 = expm1(-a*d) + df_alpha0 * exp(logEmiss[2*n] - logEmiss[2*n+1] - a*d);
      } else {
        df_alpha1 += (1 + df_alpha0*exp(logEmiss[2*n] - logEmiss[2*n+1]));
        df_alpha0 = -1;
      }
    }

    scalar_t lik = alpha1 + logEmiss[2*N-1];
    scalar_t df_lik = df_alpha1 + df_alpha0 * exp(logEmiss[2*N-2]-logEmiss[2*N-1]);
    scalar_t da_lik = 0;
  
    grad[0] = scale*da_lik;
    grad[1] = scale*df_lik;
    return scale*lik;
  }

};
