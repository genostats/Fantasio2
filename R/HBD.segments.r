#'Find HBD segments
#' 
#'This function creates HBD segments used for creating plots                      
#'
#' @param submaps a atlas object
#' @param n.consecutive.markers the number of consecutive markers with a probabilitie equal or greater to the value of the threshold, to be use to fing HBD segments (default is 5)
#' @param threshold the minimum value of HBD probabilities given for a marker (default is 0.5)
#' 
#' @details The threshold is the minimum value from which we consider the marker is HBD. From this marker we want a minumum number of consecutive markers to create a HBD segment.
#' (argument : n.consecutive.markers)
#' 
#' @return This function returns a list of dataframe with 11 columns : 
#' @return - id
#' @return - famid
#' @return - pheno
#' @return - start
#' @return - end
#' @return - size
#' @return - chromosome
#' @return - start_pos
#' @return - end_pos
#' @return - start_dist
#' @return - end_dist
#' 
#' @keywords internal
HBD.segments <- function(submaps, n.consecutive.markers = 5, threshold = 0.5)
{
  
  if(submaps@recap == 'segment')
  {
    HBD.segments.by.segments(submaps, submaps@HBD_recap, n.consecutive.markers, threshold)
  }
  else
  {
    HBD.segments.by.snps(submaps, submaps@HBD_recap, n.consecutive.markers, threshold)
  }
}