#include "logitModel.h"
#include "getUserParam.h"
#include "RMatrix.h"
#include "houba/MMatrix.h"

template<typename scalar_t>
using MATRIX = Eigen::Matrix<scalar_t, Eigen::Dynamic, Eigen::Dynamic>;

template<typename scalar_t>
using VECTOR = Eigen::Matrix<scalar_t, Eigen::Dynamic, 1>;

using namespace Rcpp;
using namespace Eigen;

// la dernière colonne de X doit être 'vide' (elle va servir à copier 
// les colonnes de H une à une)
//[[Rcpp::export]]
List logitModel_matrix(NumericVector Y, NumericMatrix X, NumericMatrix H, unsigned int beg, unsigned int end, 
                       bool centered, NumericVector rowMeansH) {
  if(getUserParam<double>().use_float) 
    return logitModel<float, Rcpp::NumericMatrix>(Y, X, H, beg, end, centered, rowMeansH);
  else 
    return logitModel<double, Rcpp::NumericMatrix>(Y, X, H, beg, end, centered, rowMeansH);
}

//[[Rcpp::export]]
List logitModel_mmatrix(NumericVector Y, NumericMatrix X, S4 H, unsigned int beg, unsigned int end,
                        bool centered, NumericVector rowMeansH) {
  
  std::string datatype = Rcpp::as<std::string>(H.slot("datatype"));

  if(datatype == "float") {
    Rcpp::XPtr<houba::MMatrix<float>> mmH(H.slot("ptr"));
    return logitModel<float>(Y, X, *mmH, beg, end, centered, rowMeansH);
  } else if(datatype == "double") {
    Rcpp::XPtr<houba::MMatrix<double>> mmH(H.slot("ptr"));
    return logitModel<double>(Y, X, *mmH, beg, end, centered, rowMeansH);
  } else {
    stop("datatype must be double of float");
  }
}
