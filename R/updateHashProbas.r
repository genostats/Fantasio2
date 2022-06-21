# keep = les individus avec une bonne qualité (a < 1)
updateHashProbas <- function(h, submap, keep, HBD, FLOD) {
  submap <- as.character(submap)
  for(i in seq_along(submap)) {
    snp <- as.character(submap[i])
    if(exists(snp, envir = h)) { # updates
      x <- h[[snp]]
      x[["HBD"]] <- x[["HBD"]] + ifelse(keep, HBD[i,], 0)
      x[["FLOD"]] <- x[["FLOD"]] + ifelse(keep, FLOD[i,], 0)
      x[["n"]] <- x[["n"]] + keep
      h[[snp]] <- x
    } else { # creates
      h[[snp]] <- list( HBD = ifelse(keep, HBD[i,], 0), 
                        FLOD = ifelse(keep, FLOD[i,], 0), 
                        n = as.integer(keep) )
    }
  }
  h
}
