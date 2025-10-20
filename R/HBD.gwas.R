#' Computation of HFLOD scores for HBD-GWAS
#' 
#' This function is used to compute HFLOD scores on individuals in a sample for the HBD-GWAS method
#' 
#' @param atlas a atlas object
#' @param phen the phenotype (default is the bed.matrix phenotype)
#' @param phen.code phenotype coding :
#'        - 'R' : 0:control ; 1:case ; NA:unknown (default)
#'        - 'plink' : 1:control ; 2:case ; 0/-9/NA:unknown
#' if 'plink' the function automatically convert it to 'R' to run logistic regression description
#' 
#' @return the atlas object with its slot HFLOD completed
#' @export

HBD.gwas <- function(atlas, phen, phen.code = c("plink", "R"))
{
  #phenotype
  phen.code <- match.arg(phen.code)
  
  #id <- sub(".*:", "" , row.names(expl.var))
  #id.index <- match ( id, x@bedmatrix@ped$id )
  if(!missing(phen)){
    atlas@submap_summary$pheno <- phen
  }
  
  w.id <- which.inbreds(atlas@submap_summary, phen.code = phen.code)$HFLOD
  
  HFLOD <- get.positions(atlas)
  
  HFLOD_value <- numeric(nrow(HFLOD))
  ALPHA_value <- numeric(nrow(HFLOD))
  
  h <- function(alpha, flods) {
    sum(log10(alpha * 10**flods + (1-alpha) ), na.rm = TRUE)
  }
  for (j in seq_len(nrow(HFLOD))) {
    # optimisation of h(alpha) ; 
    res <- optimize( h, c(0,1), flods = atlas@FLOD_recap[w.id, j], maximum = TRUE, tol = 0.001 )
    
    HFLOD_value[j] <- res$objective # HFLOD = h(alpha max)
    ALPHA_value[j] <- res$maximum   # alpha max 
  }
  HFLOD$HFLOD <- HFLOD_value
  HFLOD$ALPHA <- ALPHA_value
  
  HFLOD
}