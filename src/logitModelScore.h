#include "milorGWAS/logit_model_score.h"
#include <Rcpp.h>
#include <RcppEigen.h>
#include "getUserParam.h"
using namespace Rcpp;

template<typename scalar_t>
using MATRIX = Eigen::Matrix<scalar_t, Eigen::Dynamic, Eigen::Dynamic>;

template<typename scalar_t>
using VECTOR = Eigen::Matrix<scalar_t, Eigen::Dynamic, 1>;

// type des matrices / vecteurs construits à partir d'une NumericMatrix
// si scalar_t est double : Map Matrix ou Map Vector
// si scalar_t est float (ou autre) : Eigen Matrix
// le nom est pourri mais en pratique dans les template on met 'auto'
template<typename scalar_t>
using MATRIX_ = std::conditional_t<std::is_same_v<scalar_t, double>, 
                                  Eigen::Map<MATRIX<double>>, 
                                  MATRIX<scalar_t>>;

template<typename scalar_t>
using VECTOR_ = std::conditional_t<std::is_same_v<scalar_t, double>,
                                  Eigen::Map<VECTOR<double>>,
                                  VECTOR<scalar_t>>;

// le "getter" pour récuper une MATRIX_ à partir d'une NumericMatrix
template<typename scalar_t>
MATRIX_<scalar_t> get_matrix(Rcpp::NumericMatrix A, int p, int n);

template<>
MATRIX_<float> get_matrix<float>(Rcpp::NumericMatrix A, int p, int n) {
  // l'opérateur .cast permet de faire la copie à partir d'une Map Matrix, pas besoin d'écrire de boucle
  return MATRIX<float>(Eigen::Map<MATRIX<double>>(A.begin(), p, n).template cast<float>());
}

template<>
MATRIX_<double> get_matrix<double>(Rcpp::NumericMatrix A, int p, int n) {
  // on renvoie une map matrix
  return Eigen::Map<MATRIX<double>>(A.begin(), p, n);
}

// le "getter" pour récuper un VECTOR_ à partir d'un NumericVector
template<typename scalar_t>
VECTOR_<scalar_t> get_vector(Rcpp::NumericVector V, int n);

template<>
VECTOR_<float> get_vector<float>(Rcpp::NumericVector V, int n) {
   return VECTOR<float>(Eigen::Map<VECTOR<double>>(V.begin(), n).template cast<float>());
}

template<>
VECTOR_<double> get_vector<double>(Rcpp::NumericVector V, int n) {
   return Eigen::Map<VECTOR<double>>(V.begin(), n);
}
// ----------------------------------------------------------------------

// Y1, W, A : cf logit_model_score.h
// H la matrice dont on va tester les colonnes (de beg à end) une à une
template<typename scalar_t, typename matrixType>
List logitModelScore(NumericVector Y1, NumericVector W, NumericMatrix A, matrixType & H, unsigned int beg, unsigned int end) {
  int n = Y1.size();
  int p = A.nrow();
  if(n != A.ncol() | n != W.size() | n != H.nrow()) stop("Dimensions mismatch");

  // paramètres
  userParam<scalar_t> pars = getUserParam<scalar_t>();

  auto y1 = get_vector<scalar_t>(Y1, n);
  auto w  = get_vector<scalar_t>(W, n);
  auto a  = get_matrix<scalar_t>(A, p, n);

  // pour les résultats [thread safe vectors!]
  // on met des double parce que ça finit par un wrap()
  VECTOR<double> SCORE(end-beg+1);
  VECTOR<double> VARIANCE(end-beg+1);

  bool printed = true;
#pragma omp parallel for firstprivate(printed) num_threads(pars.n_threads)
  for(unsigned int i = beg; i <= end; i++) {
    if(!printed) {
      std::cout << "thread " << omp_get_thread_num() << "\n";
      printed = true;
    }
    scalar_t score, variance;

    // et encore une copie
    VECTOR<scalar_t> G(n);
    for(unsigned int k = 0; k < n; k++) G[k] = (scalar_t) H(k, i);

    logistic_model_score<scalar_t>(y1, G, w, a, score, variance);
    SCORE(i-beg) = (double) score;
    VARIANCE(i-beg) = (double) variance;
  }

  // on renvoie ça.
  List R;
  R["score"] = wrap(SCORE);
  R["variance"] = wrap(VARIANCE);
  return R;
}

/*********************************************************************************
 *                                                                               *
 *                              INTERCEPT ONLY                                   *
 *                                                                               *
 *********************************************************************************/

template<typename scalar_t, typename matrixType>
List logitModelScore_nocovar(NumericVector Y1, scalar_t w, matrixType & H, unsigned int beg, unsigned int end, bool compute_variance) {
  int n = Y1.size();
  if(n != H.nrow()) stop("Dimensions mismatch");

  // paramètres
  userParam<scalar_t> pars = getUserParam<scalar_t>();

  auto y1 = get_vector<scalar_t>(Y1, n);

  // pour les résultats [thread safe vectors!]
  // on met des double parce que ça finit par un wrap()
  VECTOR<double> SCORE(end-beg+1);
  VECTOR<double> VARIANCE(end-beg+1);

  bool printed = true;
#pragma omp parallel for firstprivate(printed) num_threads(pars.n_threads) 
  for(unsigned int i = beg; i <= end; i++) {
    if(!printed) {
       std::cout << "thread " << omp_get_thread_num() << "\n";
       printed = true;
    }
    scalar_t score, variance = 0;

    // et encore une copie
    VECTOR<scalar_t> G(n);
    for(unsigned int k = 0; k < n; k++) G[k] = (scalar_t) H(k, i);

    logistic_model_score_nocovar<scalar_t>(y1, G, w, score, variance, compute_variance);
    SCORE(i-beg) = (double) score;
    if(compute_variance) VARIANCE(i-beg) = (double) variance;
  }

  // on renvoie ça.
  List R;
  R["score"] = wrap(SCORE);
  if(compute_variance) R["variance"] = wrap(VARIANCE);
  return R;
}
