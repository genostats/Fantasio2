#' Wrapper for the package Fantasio
#' 
#' This function is used as a wrapper for the package Fantasio to create segments list and submaps
#' 
#' @param bedmatrix a bed.matrix 
#' @param segment.options a list of arguments to the function that will create the segments list
#' @param n the number of submaps (default is 100)
#' @param min.quality minimal quality (in \%) to include an inbred individual into the analysis (default is 95)
#' @param list.id a list of individuals of interest ('famid:id') (default = no list id)
#' @param allele.freq a vector of allele frequencies (for allele A2), if \code{bedmatrix@p} isn't appropriate
#' @param recap if you want the summary of probabilities by snps or by segments (only by SNPs for the moment)
#' @param q assumed frequency of the mutation involved in the disease for each individual (default is 0.0001)
#' @param epsilon genotype error rate (default is 0.001)
#' @param median define the f and a parameters used to compute pHBD and FLOD
#'	   - if FALSE : f and a estimated on each submap
#'	   - if TRUE : median value of estimations on all submaps of f and a (default) 
#' @param dense.recap define how pHBD and FLOD are combined over the submaps
#'	   - if FALSE : pHBD and FLOD are combined SNP by SNP with a mean on the number of submaps that includes this SNP
#'	   - if TRUE : mean of pHBD and FLOD on all submaps (default)



#' @details This function is a wrapper to make the usage of the package easier. The function calls different functions: 
#' @details The first function, `segments.list.by.hotspots` is used to create a list of segments. 
#' @details The second function, `atlas` is used to create submaps based on recombination hotspots (for the moment).
#' @details The arguments that can be included in `segment.options` are described in `segments.list.by.hotspots`.
#' @details If `recap = 'SNP'`, the quantities such as HBD probabilities, FLOD, HFLOD, 
#'   are recapitulated SNP by SNP (default).


#' @export Fantasio


# pour l'instant, que "by hotspots" avec un summary "by SNPs"
Fantasio <- function(bedmatrix, segment.options, n = 100, min.quality = 95, allele.freq, 
                     recap = c("SNP", "segment"), q = 1e-4, 
                     epsilon = 1e-3, epsilon2 = 1e-3, median = TRUE, dense.recap = TRUE) {

  recap <- match.arg(recap)
  if(recap != "SNP") stop("Not yet implemented")

  if (!missing(allele.freq)) {
    if(length(allele.freq) != ncol(bedmatrix)) {
      stop("allele.freq length should be equal to the number of SNPs in bedmatrix")
    }
    if(any(allele.freq < 0, na.rm = T) | any(allele.freq > 1, na.rm = T)) {
      stop("allele frequencies should be between 0 and 1")
    }
    bedmatrix@p <- allele.freq
  }
  
  #NA [est-ce nécessaire ?]
  if(any(is.na(bedmatrix@p))) {
    bedmatrix <- bedmatrix[, !is.na(bedmatrix@p) ]
  }
  #p=1
  bedmatrix@p <- ifelse(bedmatrix@p == 1, 1 - epsilon2, bedmatrix@p) 
  #p=0
  bedmatrix@p <- ifelse(bedmatrix@p == 0, epsilon2, bedmatrix@p) 

  if (missing(segment.options))
    segment.options <- list()
  segments.list <- do.call(segments.list.by.hotspots, c(bedmatrix = bedmatrix, segment.options))

  # le constructeur atlas() fait à peu près ce que faisait make Atlas suivi de festim
  # les slots "remplis" sont bedmatrix, seeds, epsilon, segments_list, estimations, submap_summary
  x <- atlas(bedmatrix, segments.list, n, min.quality, epsilon)
  
  # détermine les indices des individus sur lesquels on calcule HBD et FLOD 
  # et parmi ceux ci lequels seront à prendre en compte pour le HFLOD
  # (cas consanguins ou autre selon les valeurs de min.quality list.id et phen.code...)
  indexes <- which( x@submap_summary$inbred )
  keep.inds <- seq_len(nrow(bedmatrix)) %in% indexes

  # ceci remplit HBD_recap et FLOD_recap
  if(dense.recap)
    x <- recap.HBD.FLOD.dense(x, keep.inds, q, recap, median)
  else
    x <- recap.HBD.FLOD.sparse(x, keep.inds, q, recap, median)

  x@HBDsegments <- HBD.segments(x, n.consecutive.markers = 5, threshold = 0.5)
  
  x
}


