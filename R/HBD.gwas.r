#' Computation of HFLOD scores for HBD-GWAS
#' 
#' This function is used to compute HFLOD scores on individuals in a sample for the HBD-GWAS method
#' 
#' @param atlas an atlas object
#' @param phen the phenotype (default is the bed.matrix phenotype)
#' @param phen.code phenotype coding :
#'        - 'R' : 0:control ; 1:case ; NA:unknown (default)
#'        - 'plink' : 1:control ; 2:case ; 0/-9/NA:unknown
#' if 'plink' the function automatically convert it to 'R' to run logistic regression description
#' 
#' @return the atlas object with its slot HFLOD completed
#' @export

HBD.gwas <- function(atlas, phen, phen.code = c("R", "plink")) {
  # phenotype
  phen.code <- match.arg(phen.code)
  
  if(!missing(phen)){
    atlas@submap_summary$pheno <- phen
  }
  
  w.id <- which.inbreds(atlas@submap_summary, phen.code = phen.code)$HFLOD # which.inbreds gives the id of inbred cases
  
  HFLOD <- get.positions(atlas)
  
  HFLOD_value <- numeric(nrow(HFLOD))
  ALPHA_value <- numeric(nrow(HFLOD))
 
  X <- apply(atlas@FLOD_recap[w.id,], 2, hflod, eps =  0.01) 

  HFLOD$ALPHA <- X[1,]
  HFLOD$HFLOD <- X[2,]
  
  HFLOD
}
