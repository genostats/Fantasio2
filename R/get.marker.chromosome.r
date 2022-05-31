get.marker.chromosome <- function(chrSegmentsList)
{
  #find the markers index in the interval
  submap <- numeric(length(chrSegmentsList))
  
  #Step 4 : choose a random mkr
  for( i in seq_along(chrSegmentsList))  # throughout the segment
  {
    if(length(chrSegmentsList[[i]]) == 0) next  # empty segment
    if(length(chrSegmentsList[[i]]) == 1) 
    {
      s <- chrSegmentsList[[i]]
    } else { 
      s <- sample(chrSegmentsList[[i]][1]:chrSegmentsList[[i]][2], 1)
    }
    submap[i] <- s
  }
  
  return(submap)
}