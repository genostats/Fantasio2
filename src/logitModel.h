#include "milorGWAS/logit_model.h"
#include <Rcpp.h>
#include <RcppEigen.h>
#include "getUserParam.h"
using namespace Rcpp;

template<typename scalar_t>
using MATRIX = Eigen::Matrix<scalar_t, Eigen::Dynamic, Eigen::Dynamic>;

template<typename scalar_t>
using VECTOR = Eigen::Matrix<scalar_t, Eigen::Dynamic, 1>;


// la dernière colonne de X doit être 'vide' (elle va servir à copier 
// les colonnes de H une à une)
template<typename scalar_t>
List logitModel(NumericVector Y, NumericMatrix X, NumericMatrix H, unsigned int beg, unsigned int end) {
  int n = Y.size();
  int r = X.ncol();
 
  // paramètres
  userParam<scalar_t> pars = getUserParam<scalar_t>();
 
  // recopiage des matrices... nécessaire en scalar_t [spécialiser template ?]
  MATRIX<scalar_t> y(n,1);
  MATRIX<scalar_t> x(n,r);
  for(int i = 0; i < n; i++) y(i,0) = (scalar_t) Y[i];

  for(int i = 0; i < n; i++)
    for(int j = 0; j < r; j++)
      x(i,j) = (scalar_t) X(i,j);

 
  // pour les résultats [thread safe vectors!]
  // on met des double parce que ça finit par un wrap()
  VECTOR<double> BETA(end-beg+1);
  VECTOR<double> SDBETA(end-beg+1);

  // pour la régression
  VECTOR<scalar_t> beta(r);
  beta.setZero();
  MATRIX<scalar_t> varbeta(r,r);

  // les copies de beta et varbeta sont nécessaires pour qu'ils soient bien dimensionnés
  // il faut évidemment une copie privée de x pour chaque thread !
#pragma omp parallel for firstprivate(beta) firstprivate(varbeta) firstprivate(x) num_threads(pars.n_threads)
  for(unsigned int i = beg; i <= end; i++) {
    // copie de la colonne i de H dans la dernière colonne de X
    for(unsigned int k = 0; k < n; k++) 
      x(k, r-1) = H(k, i);

    logistic_model2<scalar_t>(y, x, beta, varbeta);
    BETA(i-beg) = beta(r-1);
    SDBETA(i-beg) = sqrt(varbeta(r-1,r-1));
  }

  // on renvoie ça.
  List R;
  R["beta"] = wrap(BETA);
  R["sd.beta"] = wrap(SDBETA);
  return R;
}
