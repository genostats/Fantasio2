#include "logitModelScore.h"
#include "getUserParam.h"

template<typename scalar_t>
using MATRIX = Eigen::Matrix<scalar_t, Eigen::Dynamic, Eigen::Dynamic>;

template<typename scalar_t>
using VECTOR = Eigen::Matrix<scalar_t, Eigen::Dynamic, 1>;

using namespace Rcpp;
using namespace Eigen;

//[[Rcpp::export]]
List logitModelScore(NumericVector Y1, NumericVector W, NumericMatrix A, NumericMatrix H, unsigned int beg, unsigned int end) {
  if(getUserParam<double>().use_float) 
    return logitModelScore<float>(Y1, W, A, H, beg, end);
  else 
    return logitModelScore<double>(Y1, W, A, H, beg, end);
}

//[[Rcpp::export]]
List logitModelScore_nocovar(NumericVector Y1, double w, NumericMatrix H, unsigned int beg, unsigned int end, bool compute_variance = true) {
  if(getUserParam<double>().use_float) 
    return logitModelScore_nocovar<float>(Y1, w, H, beg, end, compute_variance);
  else 
    return logitModelScore_nocovar<double>(Y1, w, H, beg, end, compute_variance);
}
