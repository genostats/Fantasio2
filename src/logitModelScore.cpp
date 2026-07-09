#include "logitModelScore.h"
#include "getUserParam.h"
#include "houba/MMatrix.h"

template<typename scalar_t>
using MATRIX = Eigen::Matrix<scalar_t, Eigen::Dynamic, Eigen::Dynamic>;

template<typename scalar_t>
using VECTOR = Eigen::Matrix<scalar_t, Eigen::Dynamic, 1>;

using namespace Rcpp;
using namespace Eigen;

//[[Rcpp::export]]
List logitModelScore_matrix(NumericVector Y1, NumericVector W, NumericMatrix A, NumericMatrix H, unsigned int beg, unsigned int end, 
                            bool centered, NumericVector rowMeansH) {
  return logitModelScore<double>(Y1, W, A, H, beg, end, centered, rowMeansH);
}

//[[Rcpp::export]]
List logitModelScore_nocovar_matrix(NumericVector Y1, double w, NumericMatrix H, unsigned int beg, unsigned int end, bool compute_variance, 
                                    bool centered, NumericVector rowMeansH) {
  return logitModelScore_nocovar<double>(Y1, w, H, beg, end, compute_variance, centered, rowMeansH);
}

//[[Rcpp::export]]
List logitModelScore_mmatrix(NumericVector Y1, NumericVector W, NumericMatrix A, S4 H, unsigned int beg, unsigned int end,
                             bool centered, NumericVector rowMeansH) {
     
  std::string datatype = Rcpp::as<std::string>(H.slot("datatype"));

  if(datatype == "float") {
    Rcpp::XPtr<houba::MMatrix<float>> mmH(H.slot("ptr"));
    return logitModelScore<float>(Y1, W, A, *mmH, beg, end, centered, rowMeansH);
  } else if(datatype == "double") {
    Rcpp::XPtr<houba::MMatrix<double>> mmH(H.slot("ptr"));
    return logitModelScore<double>(Y1, W, A, *mmH, beg, end, centered, rowMeansH);
  } else {
    stop("datatype must be double of float");
  }
}

//[[Rcpp::export]]
List logitModelScore_nocovar_mmatrix(NumericVector Y1, double w, S4 H, unsigned int beg, unsigned int end, bool compute_variance,
                                     bool centered, NumericVector rowMeansH) {
     
  std::string datatype = Rcpp::as<std::string>(H.slot("datatype"));

  if(datatype == "float") {
    Rcpp::XPtr<houba::MMatrix<float>> mmH(H.slot("ptr"));
    return logitModelScore_nocovar<float>(Y1, w, *mmH, beg, end, compute_variance, centered, rowMeansH);
  } else if(datatype == "double") {
    Rcpp::XPtr<houba::MMatrix<double>> mmH(H.slot("ptr"));
    return logitModelScore_nocovar<double>(Y1, w, *mmH, beg, end, compute_variance, centered, rowMeansH);
  } else {
    stop("datatype must be double of float");
  }
}
