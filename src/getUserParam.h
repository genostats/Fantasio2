#include <Rcpp.h>
#include <RcppEigen.h>
#include "LBFGSB.h"

#ifndef userParams
#define userParams

template<typename scalar_t>
LBFGSpp::LBFGSBParam<scalar_t> getUserParam();

template<typename scalar_t>
using VECTOR = Eigen::Matrix<scalar_t, Eigen::Dynamic, 1>;

template<typename scalar_t>
VECTOR<scalar_t> getLb();

template<typename scalar_t>
VECTOR<scalar_t> getUb();

bool debug();

#endif
