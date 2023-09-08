rsubmap <- function(segmentsList){
  
  segmentSummary <- segments.list.summary(segmentsList)
  shift <- cumsum(segmentSummary$number_of_segments)
  shift <- append(0L, shift)
  max <- shift[length(shift)]
  
  submap <- integer(max)
  
  for(chr in seq_along(segmentsList)) {
    chrMarker <- segmentsList[[chr]]
    randomMarkerVector <- random.snp(chrMarker)
    submap[(shift[chr]+1):shift[chr+1]] <- randomMarkerVector
  }
  
  submap
}

