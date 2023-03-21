#' Set phenotype
#' @param x a bedmatrix or an atlas
#' @param pheno a vector of phenotypes
#'
#' @details Update the phenotype in \code{x} and returns \code{x}
#' @return A bedmatrix or an atlas
#' @export
set.pheno <- function(x, pheno) {
  if(class(x) == "bedmatrix") {
    x@ped$pheno <- pheno
  }
  if(class(x) == "atlas") {
    x@bedmatrix@ped$pheno <- pheno
    x@submap_summary$pheno <- pheno
  }
  x
}
