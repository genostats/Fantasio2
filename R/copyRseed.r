copyRseed <- function() {
  s <- .Random.seed
  if(length(s) != 626L | s[1L] != 10403L) stop("Invalid seed")
  s <- c( s[-(1:2)], s[2L] )
  setSeed(s)
}

