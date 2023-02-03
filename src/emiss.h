#include <Rcpp.h>
#include <iostream>
#include "gaston/matrix4.h"
#include <math.h>
#include "LSE.h"
#include "getUserParam.h"
#include "RVector.h"

#ifndef __emiss__
#define __emiss__

#define SHOW(x) Rcpp::Rcout << #x << " = " << (x) << std::endl;

/* 
   Il y a deux constructeurs
   Un où on passe le pointeur de la bed matrix, le vecteur p des frequences de l'allele A2,
   la submap (vecteur des indices des SNPs qu'on conserve) et epsilon

   Un autre où on passe whichInds = les indices des indidus à garder (C BASED)
   dans ce cas partial est mis à true et on ne stocke que les log emissions des
   individus concernés... on ne fait pas de vérification dans getLogEmiss(i) !
*/

template<typename scalar_t>
class emiss {
private:
  matrix4 * bm;          // bed matrix
  RVector<double> p;     // vecteur de fréq de l'allele alt
  RVector<int> submap;   // indice des SNPs à considérer
  // !! R indices, ie à partir de 1, cf dans le code submap[i] - 1 !!
  double epsilon;        // proba d'erreur

  bool partial;
  RVector<int> whichInds; // les LogicalVectors sont des entiers
  // !! celui là est C based, ie à partir de 0... !!

  int nbInds;
  int nbSnps;
  int nbTotalSnps;
  int nbPrecomp;
  int preComputedFrom;

  std::vector< std::vector<scalar_t> > LEMISS;

  scalar_t gg[4], hh[4]; // pour les log probas d'émission à HBD 0 et 1

public:

  // preComputedFrom est initialisé à une valeur trop grande 
  // -> au premier appel de getLogEmiss il y aura appel de preCompute

  emiss(matrix4 * bm_, RVector<double> p_, RVector<int> submap_, double epsilon_) :
      bm(bm_), p(p_), submap(submap_), epsilon(epsilon_), partial(false), nbInds(bm->ncol), nbSnps(submap.size()),
      nbTotalSnps(bm->nrow), nbPrecomp(4), preComputedFrom(bm->ncol + 1) 
   {
     init();
   }

  // quand on n'est interessé que par certains individus
  emiss(matrix4 * bm_, RVector<double> p_, RVector<int> submap_, RVector<int> whichInds_, double epsilon_) :
      bm(bm_), p(p_), submap(submap_), epsilon(epsilon_), partial(true), whichInds(whichInds_), nbInds(bm->ncol), nbSnps(submap.size()),
      nbTotalSnps(bm->nrow), nbPrecomp(4), preComputedFrom(bm->ncol + 1) 
   {
     init();
   }


  // commun à tous les constructeurs
  void init() {
    LEMISS.resize(nbPrecomp);
    for(auto & le : LEMISS) {
      le.reserve(nbSnps);
    }
    // une fois pour toutes
    gg[3] = hh[3] = 0; // manquant -> proba 1
  }

  void preCompute(size_t from) {

    // clear buffer
    for(auto & le : LEMISS) {
      le.clear();
    }

    // remet from à un multiple de 4...
    from = 4 * (from/4);

    for(int i = 0; i < nbSnps; i++) { // snp loop
      if(submap[i] < 1 || submap[i] > nbTotalSnps) 
        stop("SNP index out of range");

      uint8_t * data = bm->data[submap[i]-1];
      scalar_t p_ = (scalar_t) p[submap[i]-1];

      if(std::isnan(p_) || p_ < 0 || p_ > 1) { 
        // no allelic frequency information
        // will be used in the f/b algorithm to insert SNPs in a submap without 
        // using the genotypes at these new positions
        // -> emission proba = 1 ; log emission proba = 0
        gg[0] = gg[1] = gg[2] = 0;
        hh[0] = hh[1] = hh[2] = 0;
      } else {
        // HBD 0
        gg[0] = 2*log(1-p_);                   // p^2
        gg[1] = log(2.) + log(1-p_) + log(p_); // 2pq
        gg[2] = 2*log(p_);                     // q^2

        // HBD 1
        hh[0] = log( (1-epsilon)*(1-p_) + epsilon*(1-p_)*(1-p_) );
        hh[1] = log(epsilon) + gg[1];  // epsilon*2pq
        hh[2] = log( (1-epsilon)*p_ + epsilon*p_*p_ );
      }

      size_t k = 0;
      for(size_t j = 0; j < nbPrecomp/4; j++) {
        uint8_t x = data[from/4 + j];
         for(int ss = 0; ss < 4 && (4*j + ss) < nbInds; ss++) {
           if(!partial || whichInds[from + k]) {
             LEMISS[k].push_back( gg[x&3] );
             LEMISS[k].push_back( hh[x&3] );
           }
           k++;
           x >>= 2;
         }
      }
    }
    preComputedFrom = from;

    if(getUserParam<scalar_t>().debug > 1) {
      std::cout << "[log proba emiss precomputed from " << preComputedFrom << "]\n";
    }
  }

  std::vector<scalar_t> & getLogEmiss(size_t i) {
    if(i < preComputedFrom || i >= preComputedFrom + nbPrecomp) {
      preCompute(i);
    }
    return LEMISS[i - preComputedFrom];
  }
  
};

#endif
