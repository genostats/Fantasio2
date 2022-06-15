#include <Rcpp.h>
#include <RcppEigen.h>
#include <iostream>
#include <ctime>
#include "gaston/matrix4.h"
#include "emiss.h"
#include "likelihoodGradient.h"
#include "LBFGSB.h"
#include "getUserParam.h"

#ifndef __festim__
#define __festim__

template<typename scalar_t>
List festim(XPtr<matrix4> p_A, NumericVector p_, IntegerVector map_, NumericVector deltaDist, scalar_t epsilon) {
  RVector<double> p(p_);
  RVector<int> map(map_);
  matrix4 * PA(p_A);

  std::vector<scalar_t> dDist;
  for(double a : deltaDist) 
    dDist.push_back(a);

  userParam<scalar_t> pars = getUserParam<scalar_t>();

  LBFGSpp::LBFGSBParam<scalar_t> param = pars.BFGSparam;

  VECTOR<scalar_t> lb = pars.lb;
  VECTOR<scalar_t> ub = pars.ub;

  std::vector<scalar_t> A(p_A->ncol);
  std::vector<scalar_t> F(p_A->ncol);
  std::vector<scalar_t> LIK0(p_A->ncol);
  std::vector<scalar_t> LIK1(p_A->ncol);

  emiss<scalar_t> EM(PA, p, map, epsilon);
  LBFGSpp::LBFGSBSolver<scalar_t> solver(param);
  clock_t beg = clock();

#pragma omp parallel num_threads(pars.n_threads)
#pragma omp for firstprivate(EM, solver)
  for(int i0 = 0; i0 < p_A->ncol; i0 += 4) { 
    for(int i1 = 0; i1 < 4 & i0 + i1 < p_A->ncol; i1++) {
      int i = i0 + i1;   
      likelihoodGradient<scalar_t> LG( EM.getLogEmiss(i) , dDist, -1); // (-1) is the scale parameter
      VECTOR<scalar_t> x(2);
      VECTOR<scalar_t> grad(2);
   
      x << 0, 0;
      
      LIK0[i] = -LG(x, grad);

      x << 0.05, 0.05;
      scalar_t fx;
      int count = 1;
    
      while(true) {
        try {
          solver.minimize(LG, x, fx, lb, ub);
          break;
        } catch(std::runtime_error const & e) {
          if(pars.debug > 0) 
            std::cout << e.what() << " -- x = " << x.transpose() << " [" << count << "/" << pars.max_retries << "]\n";
          // on ne ré essaie que pour ce cas ci
          if( e.what() != std::string("the line search routine reached the maximum number of iterations") || ++count > pars.max_retries) {
            break;
          }
        }
      }
      A[i] = x[0];
      F[i] = x[1];
      LIK1[i] = -fx;  // likelihood = - f(x) 
    }
  }
  if(pars.debug) {
    std::cout << (float) (clock() - beg) / CLOCKS_PER_SEC << " spent in likelihood maximization\n";
  }
  List L;
  L["a"] = wrap(A);
  L["f"] = wrap(F);
  L["likelihood.0"] = wrap(LIK0);
  L["likelihood.1"] = wrap(LIK1);
  return L;
}

#endif
