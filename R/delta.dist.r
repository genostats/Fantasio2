

#' @export 
delta.dist <- function(sx) {
  dist <- sx@snps$dist[ sx@submap ]
  chr  <- sx@snps$chr[ sx@submap ]
  delta <- diff(dist)
  # treat the change of chromosome 
  I <- cumsum(rle(chr)$length)
  I <- I[-length(I)]
  delta[I] <- -1
  delta
}

