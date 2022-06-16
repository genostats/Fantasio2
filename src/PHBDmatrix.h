#include <vector>
#include <unordered_map>
#include <RcppEigen.h>
#include "RVector.h"

#ifndef _phbdmatrix_
#define _phbdmatrix_

template<typename scalar_t>
using MATRIX = Eigen::Matrix<scalar_t, Eigen::Dynamic, Eigen::Dynamic>;

// template<typename scalar_t, template<typename> class matrix_t = MATRIX>
template<typename scalar_t>
class PHBDmatrix {
private:
  unsigned int nbSNPs; // nb de SNPs (= nb de probas HBD pour chaque individu)
  std::unordered_map<int, int> index; // individu i (conservé pour le calcul de HBD) -> indice de colonne j (dans la matrice qui va bien)
  MATRIX<scalar_t> PHBD;

public:
  // constructeur 
  // whichInds : le vecteur des indices i des individus 
  template<typename vec>
  PHBDmatrix(vec whichInds, unsigned int nbSNPs_) : nbSNPs(nbSNPs_) {
    // compter les individus / leur assigner une colonne de la matrice
    unsigned int ncol = 0;
    for(unsigned int i = 0; i < whichInds.size(); i++) {
      if(whichInds[i]) {
        index.insert( std::make_pair(i, ncol++) );
      }
      PHBD = MATRIX<scalar_t>(nbSNPs, ncol);
    }
  }

  unsigned int ncol() {
    return PHBD.cols();
  } 

  unsigned int nrow() {
    return PHBD.rows();
  } 
 
  unsigned int getColIndex(unsigned int i) {
    auto x = index.find(i);
    if(x == index.end())
      stop("Individual not in matrix");
    return x->second;
  }
 
  RVector<scalar_t> getCol(unsigned int i) {
    unsigned int j = getColIndex(i);
    return RVector<scalar_t>( &PHBD(0,0) + j*nbSNPs, &PHBD(0,0) + (j+1)*nbSNPs );
  }

  MATRIX<scalar_t> getMatrix() {
    return PHBD;
  }
};

#endif
