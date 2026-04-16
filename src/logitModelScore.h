#include "milorGWAS/logit_model_score.h"
#include <Rcpp.h>
#include <RcppEigen.h>
#include "getUserParam.h"
using namespace Rcpp;

template<typename scalar_t>
using MATRIX = Eigen::Matrix<scalar_t, Eigen::Dynamic, Eigen::Dynamic>;

template<typename scalar_t>
using VECTOR = Eigen::Matrix<scalar_t, Eigen::Dynamic, 1>;

// Y1, W, A : cf logit_model_score.h
// H la matrice dont on va tester les colonnes (de beg à end) une à une
template<typename scalar_t>
List logitModelScore(NumericVector Y1, NumericVector W, NumericMatrix A, NumericMatrix H, unsigned int beg, unsigned int end);

// float -------------------------------------------------------------------
template<>
List logitModelScore<float>(NumericVector Y1, NumericVector W, NumericMatrix A, NumericMatrix H, unsigned int beg, unsigned int end) {
  int n = Y1.size();
  int p = A.nrow();
  if(n != A.ncol() | n != W.size() | n != H.nrow()) stop("Dimensions mismatch");

  // paramètres
  userParam<float> pars = getUserParam<float>();
  
  // recopiage des matrices... nécessaire en float 
  VECTOR<float> y1(n);
  VECTOR<float> w(n);
  MATRIX<float> a(p, n);
  for(int i = 0; i < n; i++) y1(i) = (float) Y1[i];
  for(int i = 0; i < n; i++) w(i) = (float) W[i];

  for(int i = 0; i < p; i++)
    for(int j = 0; j < n; j++)
      a(i,j) = (float) A(i,j);

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
    float score, variance;

    // et encore une copie
    VECTOR<float> G(n);
    for(unsigned int k = 0; k < n; k++) G[k] = (float) H(k, i);

    logistic_model_score<float>(y1, G, w, a, score, variance);
    SCORE(i-beg) = (double) score;
    VARIANCE(i-beg) = (double) variance;
  }

  // on renvoie ça.
  List R;
  R["score"] = wrap(SCORE);
  R["variance"] = wrap(VARIANCE);
  return R;
}

// double ------------------------------------------------------------------------
template<>
List logitModelScore<double>(NumericVector Y1, NumericVector W, NumericMatrix A, NumericMatrix H, unsigned int beg, unsigned int end) {
  int n = Y1.size();
  int p = A.nrow();
  if(n != A.ncol() | n != W.size() | n != H.nrow()) stop("Dimensions mismatch");

  // paramètres
  userParam<double> pars = getUserParam<double>();
  
  // pas de recopiage, on peut faire des map
  Eigen::Map<VECTOR<double>> y1(&Y1[0], n);
  Eigen::Map<VECTOR<double>> w(&W[0], n);
  Eigen::Map<MATRIX<double>> a(&A(0, 0), p, n);

  // pour les résultats [thread safe vectors!]
  VECTOR<double> SCORE(end-beg+1);
  VECTOR<double> VARIANCE(end-beg+1);

  bool printed = true;
#pragma omp parallel for firstprivate(printed) num_threads(pars.n_threads)
  for(unsigned int i = beg; i <= end; i++) {
    if(!printed) {
       std::cout << "thread " << omp_get_thread_num() << "\n";
       printed = true;
    }
    double score, variance;
    Eigen::Map<VECTOR<double>> G(&H(0,i), n);
    logistic_model_score<double>(y1, G, w, a, score, variance);
    SCORE(i-beg) = score;
    VARIANCE(i-beg) = variance;
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

template<typename scalar_t>
List logitModelScore_nocovar(NumericVector Y1, scalar_t w, NumericMatrix H, unsigned int beg, unsigned int end, bool compute_variance);

// float -------------------------------------------------------------------
template<>
List logitModelScore_nocovar<float>(NumericVector Y1, float w, NumericMatrix H, unsigned int beg, unsigned int end, bool compute_variance) {
  int n = Y1.size();
  if(n != H.nrow()) stop("Dimensions mismatch");

  // paramètres
  userParam<float> pars = getUserParam<float>();
  
  // recopiage des matrices... nécessaire en float 
  VECTOR<float> y1(n);
  for(int i = 0; i < n; i++) y1(i) = (float) Y1[i];

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
    float score, variance = 0;

    // et encore une copie
    VECTOR<float> G(n);
    for(unsigned int k = 0; k < n; k++) G[k] = (float) H(k, i);

    logistic_model_score_nocovar<float>(y1, G, w, score, variance, compute_variance);
    SCORE(i-beg) = (double) score;
    if(compute_variance) VARIANCE(i-beg) = (double) variance;
  }

  // on renvoie ça.
  List R;
  R["score"] = wrap(SCORE);
  if(compute_variance) R["variance"] = wrap(VARIANCE);
  return R;
}

// double ------------------------------------------------------------------------
template<>
List logitModelScore_nocovar<double>(NumericVector Y1, double w, NumericMatrix H, unsigned int beg, unsigned int end, bool compute_variance) {
  int n = Y1.size();
  if(n != H.nrow()) stop("Dimensions mismatch");

  // paramètres
  userParam<double> pars = getUserParam<double>();
  
  // pas de recopiage, on peut faire des map
  Eigen::Map<VECTOR<double>> y1(&Y1[0], n);

  // pour les résultats [thread safe vectors!]
  VECTOR<double> SCORE(end-beg+1);
  VECTOR<double> VARIANCE(end-beg+1);

  bool printed = true;
#pragma omp parallel for firstprivate(printed) num_threads(pars.n_threads)  
  for(unsigned int i = beg; i <= end; i++) {
     if(!printed) {
       std::cout << "thread " << omp_get_thread_num() << "\n";
       printed = true;
    }
    double score, variance = 0;
    Eigen::Map<VECTOR<double>> G(&H(0,i), n);
    logistic_model_score_nocovar<double>(y1, G, w, score, variance, compute_variance);
    SCORE(i-beg) = score;
    if(compute_variance) VARIANCE(i-beg) = variance;
  }

  // on renvoie ça.
  List R;
  R["score"] = wrap(SCORE);
  if(compute_variance) R["variance"] = wrap(VARIANCE);
  return R;
}


