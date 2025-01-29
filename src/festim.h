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
List festim(XPtr<matrix4> p_A, NumericVector p_, IntegerVector submap_, NumericVector deltaDist, scalar_t epsilon, 
            NumericVector f_start, NumericVector a_start) {
  RVector<double> p(p_);
  RVector<int> submap(submap_);
  matrix4 * PA(p_A);

  std::vector<scalar_t> dDist;
  for(double a : deltaDist) 
    dDist.push_back(a);

  // parametres
  userParam<scalar_t> pars = getUserParam<scalar_t>();

  LBFGSpp::LBFGSBParam<scalar_t> param = pars.BFGSparam;
  VECTOR<scalar_t> lb = pars.lb;
  VECTOR<scalar_t> ub = pars.ub;

  // le solver
  LBFGSpp::LBFGSBSolver<scalar_t> solver(param);

  // stockage des résultats
  std::vector<scalar_t> A(p_A->ncol);
  std::vector<scalar_t> F(p_A->ncol);
  std::vector<scalar_t> LIK0(p_A->ncol);
  std::vector<scalar_t> LIK1(p_A->ncol);

  // l'objet qui calcule les log emissions
  emiss<scalar_t> EM(PA, p, submap, epsilon);

  clock_t beg = clock();

#pragma omp parallel num_threads(pars.n_threads)
#pragma omp for firstprivate(EM, solver)
  for(int i0 = 0; i0 < p_A->ncol; i0 += 4) {  // boucle sur les individus, de 4 en 4. Chaque thread a une copie de 'EM' -> 4 individus précalculés
    for(int i1 = 0; i1 < 4 & i0 + i1 < p_A->ncol; i1++) { // boucle de 0 à 3.
      int i = i0 + i1;   
      likelihoodGradient<scalar_t> LG( EM.getLogEmiss(i) , dDist, -1); // (-1) is the scale parameter
      VECTOR<scalar_t> x(2);
      VECTOR<scalar_t> grad(2);
   
      x << 0, 0;
      
      LIK0[i] = -LG(x, grad);

      // si a_start[i] est NAN on met a = 0, f = 0 sans faire l'EMV
      // si a_start[i] n'est pas NAN on maximise
      if(a_start[i] != a_start[i]) {
        A[i] = 0;
        F[i] = 0;
        LIK1[i] = LIK0[i];
      } else {
        // valeurs initiales
        x[0] = a_start[i];
        x[1] = f_start[i];
        scalar_t fx;
        int count = 1;
     
        // une boucle qui reprend la maximisation dans le cas où l'erreur est dans la line search...
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
