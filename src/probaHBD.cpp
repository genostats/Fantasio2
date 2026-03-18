#include "wrap_vecvec.h"
#include <Rcpp.h>
#include <RcppEigen.h>
#include <iostream>
#include "gaston/matrix4.h"
#include "getUserParam.h"
#include "PHBDmatrix.h"
#include "probaHBD.h"
#include "houba/MMatrix.h"

//[[Rcpp::export]]
SEXP probaHBD(XPtr<matrix4> p_A, NumericVector p, IntegerVector submap, NumericVector deltaDist, LogicalVector whichInds, NumericVector a, NumericVector f, double epsilon, std::string file) {
  if(getUserParam<double>().use_float) {
    PHBDmatrix<float> R = probaHBD<float>(p_A, p, submap, deltaDist, whichInds, a, f, epsilon, file);
    return( XPtr<houba::MMatrix<float>>(R.getMatrix()) );
  }
  else {
    PHBDmatrix<double> R = probaHBD<double>(p_A, p, submap, deltaDist, whichInds, a, f, epsilon, file);
    return( XPtr<houba::MMatrix<double>>(R.getMatrix()) );
  }
}
