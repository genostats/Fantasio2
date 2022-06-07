

#' @export 
delta.dist <- function(sx) {
  delta.dist.0(sx, sx@submap)
}

delta.dist.0 <- function(bedmatrix, submap) {
  dist <- bedmatrix@snps$dist[ submap ]
  chr  <- bedmatrix@snps$chr[ submap ]
  delta <- diff(dist)
  # treat the change of chromosome 
  I <- cumsum(rle(chr)$length)
  I <- I[-length(I)]
  delta[I] <- -1
  delta
}

