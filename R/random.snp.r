random.snp <- function(x) {
  beg <- x$beg
  end <- x$end
  submap <- integer(length(beg))
  for(i in seq_along(beg)) 
    submap[i] <- beg[i] + sample.int( end[i] - beg[i] + 1L , 1L) - 1L
  
  return(submap)
}
