hashProbasToMatrix <- function(h) {

  SNP <- names(h)
  for(snp in SNP) {
    h[[snp]][['HBD']] <- h[[snp]][['HBD']] / h[[snp]][['n']]
    h[[snp]][['FLOD']] <- h[[snp]][['FLOD']] / h[[snp]][['n']]
  }

  SNP <- SNP[ order(as.integer(SNP)) ]
  HBD <- sapply(SNP, function(snp) h[[snp]][["HBD"]])
  FLOD <- sapply(SNP, function(snp) h[[snp]][["FLOD"]])
  list(HBD = HBD, FLOD = FLOD)
}
