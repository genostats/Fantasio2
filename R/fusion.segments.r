fusion.segments <- function(x, minMarkers) {
  n <- length(x$beg)
  beg <- rep(NA_integer_, n)
  end <- rep(NA_integer_, n)
  beg[1] <- x$beg[1]
  i <- 1L
  k <- 1L
  while(TRUE) {
    for(j in k:n) {
      if( (x$end[j] - beg[i] + 1L) >= minMarkers ) {
        end[i] <- x$end[j]
        if(j < n) { 
          beg[i+1L] <- x$beg[j+1L]
          i <- i + 1L
          k <- j + 1L
          next
        } else { # fin 
          return( list(beg = beg[1:i], end = end[1:i]) )
        }
      }
    }
    # si on est ici le dernier segment est trop court et il faut recommencer la fusion dans l'autre sens ...
    end[i] <- x$end[j]
    for(j in i:1) {
      if( (end[i] - beg[j] + 1L) >= minMarkers ) {
        beg[i] <- beg[j]
        break
      }
    }
    # au pire on prend tout 
    if(j == 1)
      beg[i] <- beg[1]
    # il faut se débarasser du des segments qu'on a rattachés au dernier
    beg <- beg[1:i]
    end <- end[1:i]
    w <- which(beg >= beg[i])
    if(length(w) > 1) {
      w <- w[-length(w)]
      beg <- beg[-w]
      end <- end[-w]
    }
    # et on renvoie
    return(list(beg = beg, end = end))
  }
}
