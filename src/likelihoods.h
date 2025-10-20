#include <Rcpp.h>
#include <RcppEigen.h>
#include <iostream>
#include <ctime>
#include "gaston/matrix4.h"
#include "emiss.h"
#include "likelihood.h"
#include "LBFGSB.h"
#include "getUserParam.h"

#ifndef __likelihoods__
#define __likelihoods__

// calcule les vraisemblances en parallèle...
template<typename scalar_t>
NumericVector likelihoods(XPtr<matrix4> p_A, NumericVector p_, IntegerVector submap_, NumericVector deltaDist, scalar_t epsilon, NumericVector a_, NumericVector f_) {
  RVector<double> p(p_);
  RVector<int> submap(submap_);
  RVector<double> a(a_);
  RVector<double> f(f_);
  matrix4 * PA(p_A);

  std::vector<scalar_t> dDist;
  for(double a : deltaDist) 
    dDist.push_back(a);

  // parametres
  userParam<scalar_t> pars = getUserParam<scalar_t>();

  // stockage des résultats
  std::vector<scalar_t> LIK(p_A->ncol);

  // l'objet qui calcule les log emissions
  emiss<scalar_t> EM(PA, p, submap, epsilon);

  clock_t beg = clock();

#pragma omp parallel num_threads(pars.n_threads)
#pragma omp for firstprivate(EM)
  for(int i0 = 0; i0 < p_A->ncol; i0 += 4) {  // boucle sur les individus, de 4 en 4. Chaque thread a un copie de 'EM' -> 4 individus précalculés
    for(int i1 = 0; i1 < 4 & i0 + i1 < p_A->ncol; i1++) { // boucle de 0 à 3.
      int i = i0 + i1;   
      // objet qui calcule la vraisemblance p
      likelihood<scalar_t> LG( EM.getLogEmiss(i) , dDist); 
      LIK[i] = LG(a[i], f[i]);
    }
  }
  if(pars.debug) {
    std::cout << (float) (clock() - beg) / CLOCKS_PER_SEC << " spent in likelihood computation\n";
  }
  return wrap(LIK);
}

#endif
