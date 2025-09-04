#include <Rcpp.h>
#include <iostream>
#include <random>
#include <omp.h>

using namespace Rcpp;
// [[Rcpp::plugins(openmp)]]

// R::qnorm should be thread safe
// Ceci semble plus rapide que std::normal_distribution ...
// !!! ATTENTION à bien prendre une référence au rng sinon on le copie et c'est la cata...!
template<typename TRNG>
inline double rnormal(TRNG & rng, double mu, double sigma) {
  double U = ((double) rng() ) / ((double) rng.max()) ;
  return R::qnorm(U, mu, sigma, 1, 0);
}

// nSimus = nombre de génomes simulés
// lambda, sigma, mu : paramètres du processus gaussien
// dist = distance entre les points simulés
// len = longueur des chromosomes
// seed = une graine pour le générateur aléatoire. Prendre seed = runif(1) * 4294967296 par exemple
// nThreads = nb de threads
//
// [[Rcpp::export]]
NumericVector maxGaussianGenome(int nSimus, double lambda, double sigma, double dist, NumericVector len, unsigned int seed, int nThreads = 4) {
  double Phi = exp(-lambda*dist);
  double sdEps = sqrt(1 - Phi*Phi)*sigma;

  NumericVector Max(nSimus);
#pragma omp parallel num_threads(nThreads) 
  {
    // après tout on veut que ça aille vite (en plus c'est provisoire, normalement)
    // donc on prend un LCG
    std::minstd_rand rng;
    // graine différente pour chaque thread. Pas top. On pourrait avancer la longueur du thread, c'est facile car c'est un Lehmer
    // on pourrait d'ailleurs le ré-implémenter (je ne le trouve pas très rapide par rapport au MT de R)
    rng.seed( seed + 17 * omp_get_thread_num() ); 
   
    for(int s = 0; s < nSimus; s++) {
      double Z = rnormal(rng, 0, sigma);
      double M = Z;
      for(double le : len) {
        int N = (int) (le / dist);
        for(int i = 0; i < N; i++) {
          Z = Phi * Z + rnormal(rng, 0, sdEps);
          if(Z > M) M = Z;
        }
        // fin de chromosome : nouveau départ 
        Z = rnormal(rng, 0, sigma);
      }
      Max(s) = M;
    }
  } // fin d'omp parallel
  return Max;
}


