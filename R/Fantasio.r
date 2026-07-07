#' Wrapper for the package Fantasio
#' 
#' This function is used as a wrapper for the package Fantasio to create segments list and submaps
#' 
#' @param bedmatrix a bed.matrix 
#' @param segment.options a list of arguments to the function that will create the segments list
#' @param n the number of submaps (default is 100)
#' @param min.quality minimal quality (in \%) to include an inbred individual into the analysis (default is 95)
#' @param allele.freq a vector of allele frequencies (for allele A2), if \code{bedmatrix@p} isn't appropriate
#' @param q assumed frequency of the mutation involved in the disease for each individual (default is 0.0001)
#' @param epsilon genotype error rate (default is 0.001)
#' @param epsilon2 shift in border allele frequency p (p = 1 will be changed to 1-epsilon2 ; p = 0 will be changed to epsilon2)
#' @param median define the f and a parameters used to compute pHBD and FLOD
#'	   - if FALSE : f and a estimated on each submap
#'	   - if TRUE : median value of estimations on all submaps of f and a (default) 
#' @param dense.recap define how pHBD and FLOD are combined over the submaps
#'	   - if FALSE : pHBD and FLOD are combined SNP by SNP with a mean on the number of submaps that includes this SNP
#'	   - if TRUE : mean of pHBD and FLOD on all submaps (default)



#' @details This function is a wrapper to make the usage of the package easier. The function calls different functions: 
#' @details The first function, `segments.list.by.hotspots` is used to create a list of segments. 
#' @details The second function, `atlas` is used to create submaps based on recombination hotspots.
#' @details The arguments that can be included in `segment.options` are described in `segments.list.by.hotspots`.


#' @export Fantasio


# pour l'instant, que "by hotspots" avec un summary "by SNPs"
Fantasio <- function(bedmatrix, segment.options, n = 100, min.quality = 95, allele.freq, q = 1e-4, 
                     epsilon = 1e-3, epsilon2 = 1e-3, median = TRUE, dense.recap = TRUE, basename) {

  # check if file exist before anything
  if(!missing(basename)) {
    if(!dense.recap) 
      warning("With dense.recap = FALSE, you can't use memory mapped matrices")
    else
      check_file_exist(basename)
  }

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

  verbose <- Fantasio.parameters("verbose")
  if(verbose) cat("* Calling segments.list.by.hotspots\n")
  segments.list <- do.call(segments.list.by.hotspots, c(bedmatrix = bedmatrix, segment.options))

  # le constructeur atlas() fait à peu près ce que faisait make Atlas suivi de festim
  # les slots "remplis" sont bedmatrix, seeds, epsilon, segments_list, estimations, submap_summary
  if(verbose) cat("\n* Calling atlas\n")
  x <- atlas(bedmatrix, segments.list, n, min.quality, epsilon)
  
  # détermine les indices des individus sur lesquels on calcule HBD et FLOD 
  # et parmi ceux ci lequels seront à prendre en compte pour le HFLOD
  # (cas consanguins ou autre selon les valeurs de min.quality list.id et phen.code...)
  indexes <- which( x@submap_summary$inbred )
  keep.inds <- seq_len(nrow(bedmatrix)) %in% indexes

  # ceci remplit HBD_recap et FLOD_recap

  if(verbose) cat("\n* Computing HBD and FLOD matrices\n")
  if(dense.recap)
    x <- recap.HBD.FLOD.dense(x, keep.inds, q, median, basename)
  else
    x <- recap.HBD.FLOD.sparse(x, keep.inds, q, median)

  if(verbose) cat("\n* Construction of HBD segments (5 consecutive markers with threshold > 0.5)\n")
  x@HBD_segments <- HBD.segments(x, n.consecutive.markers = 5, threshold = 0.5)
  
  x
}


