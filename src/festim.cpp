#include <Rcpp.h>
#include <RcppEigen.h>
#include <iostream>
#include "gaston/matrix4.h"
#include "emiss.h"
#include "likelihoodGradient.h"
#include "LBFGSB.h"
#include "getUserParam.h"

//[[Rcpp::export]]
List festim(XPtr<matrix4> p_A, NumericVector p_, IntegerVector map_, NumericVector deltaDist, double epsilon) {
  RVector<double> p(p_);
  RVector<int> map(map_);
  matrix4 * PA(p_A);

  std::vector<double> dDist;
  for(double a : deltaDist) 
    dDist.push_back(a);

  userParam<double> pars = getUserParam<double>();

  LBFGSpp::LBFGSBParam<double> param = pars.BFGSparam;

  VECTOR<double> lb = pars.lb;
  VECTOR<double> ub = pars.ub;

  std::vector<int> NITER(p_A->ncol);
  std::vector<double> A(p_A->ncol);
  std::vector<double> F(p_A->ncol);
  std::vector<double> FX(p_A->ncol);

    emiss<double> EM(PA, p, map, epsilon);
    LBFGSpp::LBFGSBSolver<double> solver(param);
#pragma omp parallel num_threads(pars.n_threads)
#pragma omp for firstprivate(EM, solver)
  for(int i0 = 0; i0 < p_A->ncol; i0 += 4) { 
    for(int i1 = 0; i1 < 4 & i0 + i1 < p_A->ncol; i1++) {
      int i = i0 + i1;   
      likelihoodGradient<double> LG( EM.getLogEmiss(i) , dDist, -1);
      VECTOR<double> x(2);
      x << 0.05, 0.05;
      double fx;
  
      int niter = 0;
      int count = 1;
    
      while(true) {
        try {
          niter += solver.minimize(LG, x, fx, lb, ub);
          break;
        } catch(std::runtime_error const & e) {
          niter++;
          if(pars.debug) 
            std::cout << e.what() << " -- x = " << x.transpose() << " [" << count << "/" << pars.max_retries << "]\n";
          if(count++ > pars.max_retries) {
            // std::cerr << e.what() << "\n";
            break;
          }
        }
      }
      NITER[i] = niter;
      A[i] = x[0];
      F[i] = x[1];
      FX[i] = fx;
    }
  }
  List L;
  L["n.iter"] = wrap(NITER);
  L["a"] = wrap(A);
  L["f"] = wrap(F);
  L["fx"] = wrap(FX);
  return L;
}
