#include <Rcpp.h>
#include <RcppEigen.h>
#include <iostream>
#include <ctime>
#include <math.h>
#include "LSE.h"
#include "getUserParam.h"

#ifndef __likelihood__
#define __likelihood__

#ifndef SHOW
#define SHOW(x) Rcpp::Rcout << #x << " = " << (x) << std::endl;
#endif

using namespace Rcpp;

template<typename scalar_t>
using VECTOR = Eigen::Matrix<scalar_t, Eigen::Dynamic, 1>;


// compute the log likelihood for 1 individual, given the log emissions...
template<typename scalar_t>
class likelihood {
private:
  std::vector<scalar_t> & logEmiss;
  std::vector<scalar_t> & deltaDist;
  // NumericVector deltaDist;
public:
  
  likelihood(std::vector<scalar_t> & logEmiss_, std::vector<scalar_t> & deltaDist_) : logEmiss(logEmiss_), deltaDist(deltaDist_) {}
  
  scalar_t operator()(scalar_t a, scalar_t f) {
    userParam<scalar_t> pars = getUserParam<scalar_t>();
    
    if(pars.debug > 1) {
      std::cout << "a = " << a << ", f = " << f;
    }
    
    double fx;
    if(f == 0) 
      fx = f0(a);
    else if(f == 1)
      fx = f1(a);
    else
      fx = ff(a, f);

    if(pars.debug > 1) {
      std::cout << ", value = " << fx << "\n";
    }

    return fx;
  }
 
private:
 
  // pour 0 < f < 1
  scalar_t ff(scalar_t a, scalar_t f) {
    scalar_t lt00, lt01, lt10, lt11; // log proba transition

    // int N = logEmiss.ncol();
    int N = deltaDist.size() + 1;
    // état initial
    scalar_t alpha0 = log1p(-f);
    scalar_t alpha1 = log(f);

    for(int n = 0; n < N-1; n++) {
      // calcul des lt -------------------------
      scalar_t d = deltaDist[n];
      if(d < 0) { // changement de chromosomoe
        lt00 = log(1-f);
        lt01 = log(f);
        lt10 = log(1-f);
        lt11 = log(f);
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
      }
      // -------- fin calcul des lt ---------------

      scalar_t u0 = lt00 + logEmiss[2*n]   + alpha0;
      scalar_t v0 = lt10 + logEmiss[2*n+1] + alpha1;
      scalar_t u1 = lt01 + logEmiss[2*n]   + alpha0;
      scalar_t v1 = lt11 + logEmiss[2*n+1] + alpha1;

      alpha0 = LSE( u0, v0 );
      alpha1 = LSE( u1, v1 );

    }
    scalar_t u = alpha0 + logEmiss[2*N-2];
    scalar_t v = alpha1 + logEmiss[2*N-1];

    scalar_t lik = LSE(u, v);

    return lik;
  }

  // quand f == 0
  scalar_t f0(scalar_t a) {
    int N = deltaDist.size() + 1;
    // état initial
    scalar_t alpha0 = 0;
  
    for(int n = 0; n < N-1; n++) {
      scalar_t d = deltaDist[n];
      alpha0 += logEmiss[2*n];
    }
    scalar_t lik = alpha0 + logEmiss[2*N-2];

    return lik;
  }


  // f == 1
  scalar_t f1(scalar_t a) {
    int N = deltaDist.size() + 1;
    // état initial
    scalar_t alpha1 = 0;
  
    for(int n = 0; n < N-1; n++) {
      scalar_t d = deltaDist[n];
      alpha1 += logEmiss[2*n+1];
    }
    scalar_t lik = alpha1 + logEmiss[2*N-1];

    return lik;
  }
};

#endif
