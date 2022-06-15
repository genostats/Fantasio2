#include <Rcpp.h>
#include <RcppEigen.h>
#include <iostream>
#include "gaston/matrix4.h"
#include "emiss.h"
#include "forwardBackward.h"

using namespace Rcpp;

//[[Rcpp::export]]
NumericVector testForwardBackward(XPtr<matrix4> p_A, NumericVector p_, IntegerVector submap_, NumericVector deltaDist, 
                             double epsilon, int i, double a, double f) {

  RVector<double> p(p_);
  RVector<int> submap(submap_);

  LogicalVector whichInd(p_A->ncol);
  whichInd[i] = true;
  RVector<int> wI(whichInd);

  matrix4 * PA(p_A);

  std::vector<double> dDist;
  for(double a : deltaDist) 
    dDist.push_back(a);

  emiss<double> EM(PA, p, submap, wI, epsilon);
  return wrap(forwardBackward<double>( EM.getLogEmiss(i) , dDist, a, f));
}

