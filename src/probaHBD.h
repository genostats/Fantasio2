#include "HBDmatrix.h"
#include <iostream>
#include <ctime>
#include "gaston/matrix4.h"
#include "emiss.h"
#include "forwardBackward.h"
#include "getUserParam.h"

#ifndef __festim__
#define __festim__

template<typename scalar_t>
HBDMATRIX<scalar_t> probaHBD(XPtr<matrix4> p_A, NumericVector p_, IntegerVector submap_, NumericVector deltaDist, LogicalVector whichInds_, 
                       NumericVector a, NumericVector f, double epsilon) {
  RVector<double> p(p_);
  RVector<int> submap(submap_);
  RVector<int> whichInds(whichInds_);
  matrix4 * PA(p_A);

  std::vector<scalar_t> dDist;
  for(double a : deltaDist) 
    dDist.push_back(a);

  // parametres
  userParam<scalar_t> pars = getUserParam<scalar_t>();


  // l'objet qui calcule les log emissions
  emiss<scalar_t> EM(PA, p, submap, whichInds, epsilon);

  clock_t beg = clock();

  // stockage des résultats
  HBDMATRIX<scalar_t> PHBD(p_A->ncol);

#pragma omp parallel for num_threads(pars.n_threads) firstprivate(EM) 
  for(int i0 = 0; i0 < p_A->ncol; i0 += 4) {  // boucle sur les individus, de 4 en 4. Chaque thread a un copie de 'EM' -> 4 individus précalculés [au cas où...]
    for(int i1 = 0; i1 < 4 & i0 + i1 < p_A->ncol; i1++) {
      int i = i0 + i1;   
      if(!whichInds[i])
        continue;
      PHBD[i] = std::move( forwardBackward<scalar_t>( EM.getLogEmiss(i), dDist, a[i], f[i]) );
    }
  }

  if(pars.debug) {
    std::cout << "computed proba HBD (" << PHBD.size() << " inds) in ";
    std::cout << (float) (clock() - beg) / CLOCKS_PER_SEC << " secs\n";
  }

  return PHBD;
}

#endif
