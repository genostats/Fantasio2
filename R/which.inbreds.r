which.inbreds <- function(summary, list.id, phen.code) {
  if (phen.code == 'plink') {
    test <- any( summary$pheno == 2 ) 
  } else {
    test <- any( summary$pheno == 1 ) 
  }
  if(missing(list.id)) { # pas de list.id : défaut 
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
  } else { # on calcule sur les individus donnés !
    if(list.id == "all") {
      w.HBD <- seq_len(nrow(summary))
      w.HFLOD <- seq_along(w.HBD)
    } else {
      w.HBD <- match( list.id, uniqueIds(summary$famid, summary$id) )
      w.HFLOD <- seq_along(w.HBD)
    }
  }
  list(HBD = w.HBD, HFLOD = w.HFLOD)   
}

