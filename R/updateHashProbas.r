updateHashProbas <- function(h, submap, HBD, FLOD) {
  submap <- as.character(submap)
  for(i in seq_along(submap)) {
    snp <- as.character(submap[i])
    if(exists(snp, envir = h)) { # updates
      x <- h[[snp]]
      x[["HBD"]] <- x[["HBD"]] + HBD[i,]
      x[["FLOD"]] <- x[["FLOD"]] + FLOD[i,]
      x[["n"]] <- x[["n"]] + 1
      h[[snp]] <- x
    } else { # creates
      h[[snp]] <- list( HBD = HBD[i,], FLOD = FLOD[i,], n = 1)  
    }
  }
  h
}
