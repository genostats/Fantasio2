hashProbasToMatrix <- function(h) {

  SNP <- names(h)
  for(snp in SNP) {
    h[[snp]][['HBD']] <- h[[snp]][['HBD']] / h[[snp]][['n']]
    h[[snp]][['FLOD']] <- h[[snp]][['FLOD']] / h[[snp]][['n']]
  }

  SNPi <- as.integer(SNP)
  o <- order(SNPi)
  SNP <- SNP[o]
  HBD <- sapply(SNP, function(snp) h[[snp]][["HBD"]])
  FLOD <- sapply(SNP, function(snp) h[[snp]][["FLOD"]])
  list(HBD = HBD, FLOD = FLOD, snp = SNPi[o] )
}
