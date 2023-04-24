#' @export Fantasio
# pour l'instant, que "by hotspots" avec un summary "by SNPs"
Fantasio <- function(bedmatrix, segment.options, n = 100, min.quality = 95, list.id, 
                     recap = c("SNP", "segment"), phen.code = c("plink", "R"), q = 1e-4, epsilon = 1e-3) {

  phen.code <- match.arg(phen.code)

  recap <- match.arg(recap)
  if(recap != "SNP") stop("Not yet implemented")


  if (missing(segment.options))
    segment.options <- list()
  segments.list <- do.call(segments.list.by.hotspots, c(bedmatrix = bedmatrix, segment.options))

  # le constructeur atlas() fait à peu près ce que faisait make Atlas suivi de festim
  # les slots "remplis" sont bedmatrix, seeds, epsilon, segments_list, estimations, submap_summary
  x <- atlas(bedmatrix, segments.list, n, epsilon)
  
  # détermine les indices des individus sur lesquels on calcule HBD et FLOD 
  # et parmi ceux ci lequels seront à prendre en compte pour le HFLOD
  # (cas consanguins ou autre selon les valeurs de min.quality list.id et phen.code...)
  indexes <- which.inbreds(x@submap_summary, min.quality, list.id, phen.code)
  keep.inds <- seq_len(nrow(bedmatrix)) %in% indexes$HBD

  # ceci remplit HBD_recap et FLOD_recap
  x <- recap.HBD.FLOD(x, keep.inds, q, recap, median)

  # ceci remplit HFLOD
  x <- set.HFLOD(x, indexes$HFLOD)

  x@HBDsegments <- HBD.segments(x, n.consecutive.markers = 5, threshold = 0.5)
  
  x
}


