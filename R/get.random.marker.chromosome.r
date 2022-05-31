get.random.marker.chromosome <- function(chrSegmentsList)
{
  #find the markers index in the interval
  submap <- integer(length(chrSegmentsList))
  
  #choose a random mkr
  for(i in seq_along(chrSegmentsList))  # throughout the segment
  {
    seg <- chrSegmentsList[[i]]
    if(length(seg) == 0L) next  # empty segment
    if(length(seg) == 1L) 
    {
      s <- seg
    } else { 
      # s <- sample(chrSegmentsList[[i]][1]:chrSegmentsList[[i]][2], 1)
      s <- seg[1L] + sample.int(diff(seg) + 1L, 1L) - 1L
    }
    submap[i] <- s
  }
  
  return(submap)
}
