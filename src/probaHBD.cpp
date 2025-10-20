#include "wrap_vecvec.h"
#include <Rcpp.h>
#include <RcppEigen.h>
#include <iostream>
#include "gaston/matrix4.h"
#include "getUserParam.h"
#include "probaHBD.h"
#include "PHBDmatrix.h"

//[[Rcpp::export]]
NumericMatrix probaHBD(XPtr<matrix4> p_A, NumericVector p, IntegerVector submap, NumericVector deltaDist, LogicalVector whichInds, NumericVector a, NumericVector f, double epsilon) {
  if(getUserParam<double>().use_float) {
    PHBDmatrix<float> R = probaHBD<float>(p_A, p, submap, deltaDist, whichInds, a, f, epsilon);
    return( wrap(R.getMatrix()) );
  }
  else {
    PHBDmatrix<double> R = probaHBD<double>(p_A, p, submap, deltaDist, whichInds, a, f, epsilon);
    return wrap(R.getMatrix());
  }
}
