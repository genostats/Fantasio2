#include <Rcpp.h>
#include <RcppEigen.h>
#include <iostream>
#include "gaston/matrix4.h"
#include "getUserParam.h"
#include "festim.h"

//[[Rcpp::export]]
List festim(XPtr<matrix4> p_A, NumericVector p, IntegerVector submap, NumericVector deltaDist, double epsilon, 
            NumericVector f_start, NumericVector a_start) {
  if(getUserParam<double>().use_float)   
    return festim<float>(p_A, p, submap, deltaDist, epsilon, f_start, a_start);
  else
    return festim<double>(p_A, p, submap, deltaDist, epsilon, f_start, a_start);
}
