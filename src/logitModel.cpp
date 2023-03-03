#include "logitModel.h"
#include "getUserParam.h"

template<typename scalar_t>
using MATRIX = Eigen::Matrix<scalar_t, Eigen::Dynamic, Eigen::Dynamic>;

template<typename scalar_t>
using VECTOR = Eigen::Matrix<scalar_t, Eigen::Dynamic, 1>;

using namespace Rcpp;
using namespace Eigen;

// la dernière colonne de X doit être 'vide' (elle va servir à copier 
// les colonnes de H une à une)
//[[Rcpp::export]]
List logitModel(NumericVector Y, NumericMatrix X, NumericMatrix H, unsigned int beg, unsigned int end) {
  if(getUserParam<double>().use_float) 
    return logitModel<float>(Y, X, H, beg, end);
  else 
    return logitModel<double>(Y, X, H, beg, end);
}
