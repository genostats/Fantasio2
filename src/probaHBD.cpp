#include "wrap_vecvec.h"
#include <Rcpp.h>
#include <RcppEigen.h>
#include <iostream>
#include "gaston/matrix4.h"
#include "getUserParam.h"
#include "probaHBD.h"
#include "PHBDmatrix.h"
#include "RMatrix.h"
//[[Rcpp::export]]
void probaHBD_matrix(XPtr<matrix4> p_A, NumericMatrix PHBD_, NumericVector p, IntegerVector submap, NumericVector deltaDist, 
                     LogicalVector whichInds, NumericVector a, NumericVector f, double epsilon) {

  RMatrix<double> PHBD(PHBD_);
  PHBDmatrix<RMatrix<double>> M(PHBD, whichInds, submap.size());
  probaHBD(p_A, M, p, submap, deltaDist, whichInds, a, f, epsilon);
}

/*
  if(getUserParam<double>().use_float) {
    PHBDmatrix<float> R = probaHBD<float>(p_A, PHBD, p, submap, deltaDist, whichInds, a, f, epsilon);
    return( wrap(R.getMatrix()) );
  } else {
    PHBDmatrix<double> R = probaHBD<double>(p_A, PHBD, p, submap, deltaDist, whichInds, a, f, epsilon);
    return wrap(R.getMatrix());
  }
}
*/

/*
//[[Rcpp::export]]
void probaHBD_matrix(XPtr<matrix4> p_A, S4 PHBD, NumericVector p, IntegerVector submap, NumericVector deltaDist, 
                     LogicalVector whichInds, NumericVector a, NumericVector f, double epsilon) {
  std::string datatype(PHBD.slot("datatype"));
  
  if(datatype == "float") {
    Rcpp::XPtr<houba::MMatrix<float>> mmPHBD(PHBD.slot("ptr"));
    PHBDmatrix<float> R = probaHBD<float>(p_A, mmPHBD, p, submap, deltaDist, whichInds, a, f, epsilon);
    return( wrap(R.getMatrix()) );
  } else if(datatype = "double") {
    Rcpp::XPtr<houba::MMatrix<double>> mmPHBD(PHBD.slot("ptr"));
    PHBDmatrix<double> R = probaHBD<double>(p_A, mmPHBD, p, submap, deltaDist, whichInds, a, f, epsilon);
    return wrap(R.getMatrix());
  } else {
    stop("datatype must be double of float")
  }
}
*/
