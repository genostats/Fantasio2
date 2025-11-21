#' @export
likelihoods <- function(atlas, f, a) { # keep.inds, q, recap, median) {
  # shortcuts for atlas slots
  n <- ncol(atlas@seeds)
  bedmatrix <- atlas@bedmatrix
  seeds <- atlas@seeds
  segments.list <- atlas@segments_list
  epsilon <- atlas@epsilon

  nb_ind <- nrow(bedmatrix)
  f <- rep_len(f, nb_ind)
  a <- rep_len(a, nb_ind)

  LIK <- matrix(0, nrow = nb_ind, ncol = n)

  verbose <- Fantasio.parameters("verbose")
  for(i in 1:n) { # boucle sur les cartes
    if(verbose) cat("Computing likelihoods on submap", i, "\r")
    # on re génère les cartes
    setSeed(seeds[,i])
    submap <- rsubmap(segments.list)
    d.dist <- delta.dist(bedmatrix, submap)
    LIK[,i] <- likelihoods_(bedmatrix@bed, bedmatrix@p, submap, d.dist, epsilon, a, f)
  }
  if(verbose) cat("\n")
  LIK 
}
