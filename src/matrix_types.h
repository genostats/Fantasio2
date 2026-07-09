#include <Rcpp.h>
#include <RcppEigen.h>

#ifndef _gaston_matrix_types_
#define _gaston_matrix_types_

template<typename scalar_t>
using MATRIX = Eigen::Matrix<scalar_t, Eigen::Dynamic, Eigen::Dynamic>;

template<typename scalar_t>
using VECTOR = Eigen::Matrix<scalar_t, Eigen::Dynamic, 1>;

// types pour des matrices / vecteurs construits à partir d'une NumericMatrix
// si scalar_t est double : Map Matrix ou Map Vector
// si scalar_t est float (ou autre) : Eigen Matrix
// le nom est pourri mais en pratique dans les template on met 'auto'
template<typename scalar_t>
using MATRIX_ = std::conditional_t<std::is_same_v<scalar_t, double>, Eigen::Map<MATRIX<double>>, MATRIX<scalar_t>>;

template<typename scalar_t>
using VECTOR_ = std::conditional_t<std::is_same_v<scalar_t, double>, Eigen::Map<VECTOR<double>>, VECTOR<scalar_t>>;

// le "getter" pour récuper une MATRIX_ à partir d'une NumericMatrix
template<typename scalar_t>
MATRIX_<scalar_t> get_matrix(Rcpp::NumericMatrix X);

template<>
inline MATRIX_<float> get_matrix<float>(Rcpp::NumericMatrix X) {
  // l'opérateur .cast permet de faire la copie à partir d'une Map Matrix, pas besoin d'écrire de boucle
  // pas facile à relire mais on peut compter sur Eigen pour que ça soit optimisé
  return MATRIX<float>(Eigen::Map<MATRIX<double>>(X.begin(), X.nrow(), X.ncol()).template cast<float>());
}

template<>
inline MATRIX_<double> get_matrix<double>(Rcpp::NumericMatrix X) {
  // on renvoie une map matrix
  return Eigen::Map<MATRIX<double>>(X.begin(), X.nrow(), X.ncol());
}

// le "getter" pour récuper un VECTOR_ à partir d'un NumericVector
template<typename scalar_t>
VECTOR_<scalar_t> get_vector(Rcpp::NumericVector V);

template<>
inline VECTOR_<float> get_vector<float>(Rcpp::NumericVector V) {
   return VECTOR<float>(Eigen::Map<VECTOR<double>>(V.begin(), V.size()).template cast<float>());
}

template<>
inline VECTOR_<double> get_vector<double>(Rcpp::NumericVector V) {
   return Eigen::Map<VECTOR<double>>(V.begin(), V.size());
}

#endif
