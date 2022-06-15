#include "HBDmatrix.h"
#include <Rcpp.h>
#include <RcppEigen.h>
#include <iostream>
#include "gaston/matrix4.h"
#include "getUserParam.h"
#include "probaHBD.h"

//[[Rcpp::export]]
NumericMatrix probaHBD(XPtr<matrix4> p_A, NumericVector p, IntegerVector submap, NumericVector deltaDist, LogicalVector whichInds, NumericVector a, NumericVector f, double epsilon) {
  if(getUserParam<double>().use_float)   
    return wrap(probaHBD<float>(p_A, p, submap, deltaDist, whichInds, a, f, epsilon));
  else {
    return wrap(probaHBD<double>(p_A, p, submap, deltaDist, whichInds, a, f, epsilon));
  }
}
