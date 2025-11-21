#include <Rcpp.h>
#include <RcppEigen.h>
#include <iostream>
#include "gaston/matrix4.h"
#include "getUserParam.h"
#include "likelihoods.h"

//[[Rcpp::export]]
NumericVector likelihoods_(XPtr<matrix4> p_A, NumericVector p, IntegerVector submap, NumericVector deltaDist, double epsilon, NumericVector a, NumericVector f) {
  if(getUserParam<double>().use_float)   
    return likelihoods<float>(p_A, p, submap, deltaDist, epsilon, a, f);
  else
    return likelihoods<double>(p_A, p, submap, deltaDist, epsilon, a, f);
}
