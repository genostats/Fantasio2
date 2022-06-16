#include <Rcpp.h>
#include <RcppEigen.h>
#include <iostream>
#include "gaston/matrix4.h"
#include "PHBDmatrix.h"

using namespace Rcpp;

//[[Rcpp::export]]
NumericMatrix testPHBDmatrix(LogicalVector z, int nbSNPs, int i) {
  PHBDmatrix<double> A(z, nbSNPs);
  std::cout << A.ncol() << " cols " << A.nrow() << " rows\n";
  std::cout << "Col for individual " << i << " = " << A.getColIndex(i) << "\n";

  RVector<double> B = A.getCol(i);
  for(int u = 0; u < nbSNPs; u++) B[u] = R::unif_rand();

  return wrap(A.getMatrix());
}
