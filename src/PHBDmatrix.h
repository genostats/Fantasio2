#include <vector>
#include <unordered_map>
#include <RcppEigen.h>
#include "RVector.h"
#include "houba/MMatrix.h"

#ifndef _phbdmatrix_
#define _phbdmatrix_

// template<typename scalar_t>
// using MATRIX = Eigen::Matrix<scalar_t, Eigen::Dynamic, Eigen::Dynamic>;

template<typename scalar_t>
using MMatrix = houba::MMatrix<scalar_t>;

template<typename scalar_t>
class PHBDmatrix {
private:
  unsigned int nbSNPs; // nb de SNPs (= nb de probas HBD pour chaque individu)
  std::unordered_map<int, int> index; // individu i (conservé pour le calcul de HBD) -> indice de colonne j (dans la matrice qui va bien)
  MMatrix<scalar_t> * PHBD;

public:
  // constructeur 
  // whichInds : le vecteur des indices i des individus 
  template<typename vec>
  PHBDmatrix(vec whichInds, unsigned int nbSNPs_, std::string file) : nbSNPs(nbSNPs_), 
                         PHBD(new MMatrix<scalar_t>(file, nbSNPs, std::count(whichInds.begin(), whichInds.end(), true))) {
    // compter les individus / leur assigner une colonne de la matrice
    unsigned int ncol = 0;
    for(unsigned int i = 0; i < whichInds.size(); i++) {
      if(whichInds[i]) {
        index.insert( std::make_pair(i, ncol++) );
      }
      // PHBD = MMatrix<scalar_t>(file, nbSNPs, ncol);
    }
  }

  unsigned int ncol() {
    return PHBD->ncol();
  } 

  unsigned int nrow() {
    return PHBD->nrow();
  } 
 
  unsigned int getColIndex(unsigned int i) {
    auto x = index.find(i);
    if(x == index.end())
      stop("Individual not in matrix");
    return x->second;
  }
 
  RVector<scalar_t> getCol(unsigned int i) {
    unsigned int j = getColIndex(i);
    return RVector<scalar_t>( PHBD->data() + j*nbSNPs, PHBD->data() + (j+1)*nbSNPs );
  }

  MMatrix<scalar_t> * getMatrix() {
    return PHBD;
  }
};

#endif
