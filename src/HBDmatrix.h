#ifndef _wrap_vecvec_
#define _wrap_vecvec_

#ifdef Rcpp_hpp
#error Don t include Rcpp.h before this file
#endif

#include <RcppCommon.h>

template<typename scalar_t>
using HBDMATRIX = std::vector<std::vector<scalar_t>>; 

namespace Rcpp {
  template <typename scalar_t> SEXP wrap(const std::vector<std::vector<scalar_t>> & x);
}

#include <Rcpp.h>
namespace Rcpp {
  template <typename scalar_t> SEXP wrap(HBDMATRIX<scalar_t> & x) {
    size_t ncol = 0;
    size_t nrow = 0;
    for(size_t i = 0; i < x.size(); i++) {
      if(x[i].size() > 0) ncol++;
      if(nrow == 0) 
        nrow = x[i].size();
      if(x[i].size() != nrow) 
        stop("Can't wrap this in a matrix");
    }
    if(ncol == 0)
      return Rcpp::NumericMatrix(0,0);

    Rcpp::NumericMatrix X = no_init(nrow, ncol);
    size_t k = 0;
    for(size_t i = 0; i < x.size(); i++) {
      if(x[i].size() == 0) 
        continue;
      for(size_t j = 0; j < nrow; j++)
        X(j,k) = (double) x[i][j];
      k++;
    }
    return X;
  }
}

#endif
