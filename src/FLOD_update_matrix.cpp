#include "RMatrix.h"
#include "RVector.h"
#include "houba/MMatrix.h"
#include <stdexcept>
#include <cmath>		
#include <Rcpp.h>

template<typename matrixType, typename vectorType> 
void FLOD_update_matrix(matrixType & HBD, matrixType & FLOD, vectorType & f, double q_) {
  using scalar_t = typename matrixType::value_type;

  scalar_t q = (scalar_t) q_;
  unsigned int n = HBD.nrow();
  unsigned int m = HBD.ncol();
  if(n != FLOD.nrow() || m != FLOD.ncol() || m != f.size()) {
    throw std::runtime_error("Dimensions mismatch in HFLOD()");
  }

  for(unsigned int j = 0; j < m; j++) {
    scalar_t fj = (scalar_t) f[j];
    scalar_t log_denom = std::log10(fj + q * (1 - fj));
    for(unsigned int i = 0; i < n; i++) {
      FLOD(i, j) += std::log10( HBD(i, j) + q * ( 1. - HBD(i, j) ) ) - log_denom;
    }
  }
}

// [[Rcpp::export]]
void FLOD_update_matrix(Rcpp::NumericMatrix HBD, Rcpp::NumericMatrix FLOD, Rcpp::NumericVector f, double q) {
  RMatrix<double> HBD_(HBD);
  RMatrix<double> FLOD_(FLOD);
  RVector<double> f_(f);
  FLOD_update_matrix(HBD_, FLOD_, f_, q);
}

// [[Rcpp::export]]
void FLOD_update_mmatrix(Rcpp::S4 HBD, Rcpp::S4 FLOD, Rcpp::NumericVector f, double q) {
  std::string datatype = Rcpp::as<std::string>(HBD.slot("datatype"));
  RVector<double> f_(f);
  if(datatype == "float") {
    Rcpp::XPtr<houba::MMatrix<float>> mmHBD(HBD.slot("ptr"));
    Rcpp::XPtr<houba::MMatrix<float>> mmFLOD(FLOD.slot("ptr"));
    FLOD_update_matrix(*mmHBD, *mmFLOD, f_, q);
  } else if(datatype == "double") {
    Rcpp::XPtr<houba::MMatrix<double>> mmHBD(HBD.slot("ptr"));
    Rcpp::XPtr<houba::MMatrix<double>> mmFLOD(FLOD.slot("ptr"));
    FLOD_update_matrix(*mmHBD, *mmFLOD, f_, q);
  } else {
    Rcpp::stop("datatype must be double of float");
  }
}
