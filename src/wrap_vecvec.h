#ifdef Rcpp_hpp
#error Don t include Rcpp.h before this file
#endif

#ifndef _wrap_vecvec_
#define _wrap_vecvec_

#include <RcppCommon.h>
// [[Rcpp::plugins(cpp11)]]

namespace Rcpp {
  template <typename scalar_t> SEXP wrap(const std::vector<std::vector<scalar_t>> & x);
}

#include <Rcpp.h>
namespace Rcpp {
  template <typename scalar_t> SEXP wrap(const std::vector<std::vector<scalar_t>> & x) {
    if(x.size() == 0)
      return Rcpp::NumericMatrix(0,0);

    unsigned int ncol = x.size();
    unsigned int nrow = x[0].size();
    Rcpp::NumericMatrix X = no_init(nrow, ncol);
    for(unsigned int i = 0; i < ncol; i++) {
      if(x[i].size() != nrow) 
        stop("Can't wrap in a matrix");
      for(unsigned int j = 0; j < nrow; j++)
        X(j,i) = (double) x[i][j];
    }
    return X;
  }
}

#endif
