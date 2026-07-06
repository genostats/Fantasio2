#include <vector>
#include <unordered_map>
#include <RcppEigen.h>
#include "RVector.h"
#include "debug.h"

#ifndef _phbdmatrix_
#define _phbdmatrix_

template<typename matrixType>
class PHBDmatrix {
  using scalar_t = typename matrixType::value_type;

private:
  matrixType & PHBD;
  unsigned int nbSNPs; // nb de SNPs (= nb de probas HBD pour chaque individu)
  std::unordered_map<int, int> index; // individu i (conservé pour le calcul de HBD) -> indice de colonne j (dans la matrice qui va bien)

public:
  // constructeur 
  // whichInds : le vecteur des indices i des individus 
  template<typename vec>
  PHBDmatrix(matrixType & target, vec whichInds, unsigned int nbSNPs_) : PHBD(target), nbSNPs(nbSNPs_) {
    // compter les individus / leur assigner une colonne de la matrice
    unsigned int ncol = 0;
    for(unsigned int i = 0; i < whichInds.size(); i++) {
      if(whichInds[i]) {
        index.insert( std::make_pair(i, ncol++) );
      }
    }
    // verifier les dimensions
    if(PHBD.nrow() != nbSNPs || PHBD.ncol() != ncol) {
      throw std::runtime_error("Dimensions mismatch");
    }
  }

  unsigned int ncol() {
    return PHBD.ncol();
  } 

  unsigned int nrow() {
    return PHBD.nrow();
  } 
 
  unsigned int getColIndex(unsigned int i) {
    auto x = index.find(i);
    if(x == index.end())
      throw std::runtime_error("Individual not in matrix");
    return x->second;
  }
 
  RVector<scalar_t> getCol(unsigned int i) {
    unsigned int j = getColIndex(i);
    return RVector<scalar_t>( &PHBD(0,0) + j*nbSNPs, &PHBD(0,0) + (j+1)*nbSNPs );
  }
};

#endif
