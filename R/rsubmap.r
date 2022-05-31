#' @export

rsubmap <- function(segmentsList){
  
  segmentSummary <- segments.list.summary(segmentsList)
  shift <- cumsum(segmentSummary$number_of_segments)
  shift <- append(0, shift)
  max <- shift[length(shift)]
  
  submap <- numeric(max)
  
  for(chr in seq_along(segmentsList))
  {
    chrMarker <- segmentsList[[chr]]
    randomMarkerVector <- get.marker.chromosome(chrMarker)
    submap[(shift[chr]+1):shift[chr+1]] <- randomMarkerVector
  }
  
  submap
}

