
#' @export
set.submap <- function(sx, segmentsList, ...) {

  restore.seed <- getArg("restore.seed", TRUE, ...)

  # storing current seed
  old.seed <- get(".Random.seed", envir = .GlobalEnv)

  # Note. Le CRAN n'acceptera pas ça. A régler plus tard 
  # (on doit pouvoir les bluffer en passant par Rcpp)
  # the object seed
  assign(".Random.seed", sx@random.seed, envir = .GlobalEnv)

  # cette partie est à remplacer par quelque chose d'intelligent
  
  segmentSummary <- segments.list.summary(segmentsList)
  shift <- cumsum(segmentSummary$number_of_segments)
  shift <- append(0, shift)
  max <- shift[length(shift)]
  
  submap <- numeric(max)
  
  for(chr in seq_along(segmentsList))
  {
    chrMarker <- segmentsList[[chr]]
    randomMarkerVector <- get.marker.chromosom(chrMarker)
    submap[(shift[chr]+1):shift[chr+1]] <- randomMarkerVector
  }
  
  #n <- getArg("n", 100, ...)
  #submap <- sort(sample.int(ncol(sx), n))
  sx@submap <- submap
   
  # restauring seed (meme commentaire)
  if(restore.seed)
    assign(".Random.seed", old.seed, envir = .GlobalEnv)

  sx
}
