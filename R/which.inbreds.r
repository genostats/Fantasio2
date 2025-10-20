which.inbreds <- function(summary, phen.code) {
  if (phen.code == 'plink') {
    test <- any( summary$pheno == 2, na.rm = TRUE ) 
  } else {
    test <- any( summary$pheno == 1, na.rm = TRUE ) 
  }
 
  if(test) { # il y a des atteints
    # on calcule les probas HBD et les FLOD sur les individus consanguins avec qualité suffisante
    # le HFLOD sur les atteints parmi ceux là
    w.HBD   <- which( summary$inbred )
    if (phen.code == 'plink') {
      w.HFLOD <- match( which(summary$inbred & summary$pheno == 2), w.HBD )
    } else {
      w.HFLOD <- match( which( summary$inbred & summary$pheno == 1), w.HBD )
    }
  } else {
    # on calcule les probas HBD, les FLOD et les HFLOD sur tous les consanguins 
    # avec qualité
    w.HBD   <- which( summary$inbred )
    w.HFLOD <- seq_along(w.HBD)
    if (phen.code == 'plink') {
      warning("No individual with pheno = 2.\nUsing all inbred individuals with good estimation quality.")
    } else {
      warning("No individual with pheno = 1.\nUsing all inbred individuals with good estimation quality.")
    }
  }
  list(HBD = w.HBD, HFLOD = w.HFLOD)   
}

