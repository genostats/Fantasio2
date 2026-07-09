#include "milorGWAS/logit_model_score.h"
#include <Rcpp.h>
#include <RcppEigen.h>
#include "getUserParam.h"
#include "matrix_types.h"

using namespace Rcpp;

// Y1, W, A : cf logit_model_score.h
// H la matrice dont on va tester les colonnes (de beg à end) une à une
template<typename scalar_t, typename matrixType>
List logitModelScore(NumericVector Y1, NumericVector W, NumericMatrix A, matrixType & H, unsigned int beg, unsigned int end) {
  int n = Y1.size();
  if(n != A.ncol() | n != W.size() | n != H.nrow()) stop("Dimensions mismatch");

  // paramètres
  userParam<scalar_t> pars = getUserParam<scalar_t>();

  auto y1 = get_vector<scalar_t>(Y1);
  auto w  = get_vector<scalar_t>(W);
  auto a  = get_matrix<scalar_t>(A);

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

  auto y1 = get_vector<scalar_t>(Y1);

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
