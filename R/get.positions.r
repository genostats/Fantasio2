get.positions <- function(atlas) {
  if(atlas@recap != "SNP") stop("Not yet implemented")
  
  names    <- colnames(atlas@HBD_recap)  # les ids des marqueurs où on a le FLOD/ la proba HBD
  index    <- match(names, atlas@bedmatrix@snps$id)
  if(any(is.na(index))) stop("SNP ids don't match!") 
  atlas@bedmatrix@snps[ index, c("id", "chr", "pos", "dist") ]
}
