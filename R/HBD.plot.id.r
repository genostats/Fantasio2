#' plot of HBD segment 
#' 
#' This function plots HBDsegments or ROHs for all the chromosoms of a given individual
#' 
#' @param Submaps a list.submap object
#' @param ROH (optional) a data frame from which the segments will be plotted
#' @param unit the unit used to plot, "Bases" or "cM" (default is "CM")
#' @param id the individual id of the individual wanted
#' @param famid the family id of the individual wanted
#' @param regions a specific region to be highlighted in the plot (optional)
#' @param outfile a name for the plot (optional)
#' @param build the value of the genome build to use to plot chromosome in the plot (35, 36, 37,or 38, default is 37)
#' 
#' @details If `ROH` is provided, it must be a data frame containing columns `FID`, `IID`, `CHR`, `POS1`, `POS2`, `PHE`. 
#' Such a data frame can be obtained by reading the report generated by `plink` with the option `--homozyg`.
#' If `ROH` is not provided, the function will plot HBD segments defined in the Submaps object.
#' @details The `regions` argument is optional. If provided, it must be a matrix containing one line per region to be highlighted with in each line : 
#' \itemize{ 
#'   \item{the chromosome number}
#'   \item{start}
#'   \item{end}
#' }
#' 
#' @return A plot of the individual's HBDsegments.
#' 
#' @examples  
#' #Please refer to vignette 
#' 
#' @export
HBD.plot.id <- function(Submaps, ROH, unit= "cM", id, famid, regions, outfile, build = 37)
{
  if(class(Submaps@bedmatrix)[1] != "bed.matrix")
    stop("Need a bed.matrix.")
  
  if(is.null(Submaps@HBD_recap))
    stop("HBD_recap is empty cannot plot, make sure to have atleast one individual considered inbred.")
  
  if(!missing(Submaps) & !missing(ROH))
  {
    plot.ROH.segments.id(Submaps=Submaps, ROH, unit, id=id, famid=famid, regions, outfile=outfile, build=build)
  }else{
    if(!missing(Submaps))
      plot.HBD.segments.id(Submaps = Submaps, id=id, famid=famid, unit=unit, regions, outfile=outfile, build=build)
  }
}