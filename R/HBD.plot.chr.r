#' Plot the HBD segments on a given chromosome
#' 
#' This function plots all the HBD segments for a given chromosome and all the individuals
#' 
#' @param atlas an atlas object
#' @param ROH a data frame from which the segments will be plotted (optional)
#' @param unit the unit used to plot, two options are allowed "Bases", "cM" (default is "cM")
#' @param chr the chromosome number from which to plot HBD segment
#' @param list.ids a vector containing a list of individuals from which only the HBDsegments for this chromosome will be ploted (optional)
#' @param regions a specific region to be enlighted in the plot (optional)
#' @param outfile a name for the plot (optional)
#' @param inbred whether you want to plot only inbred individuals or not (default is FALSE)
#' @param build the value of the build to use to plot chromosome in the plot value accepted are 35, 36, 37, 38 (default is 37)
#' 
#' @details If you use the regions options make sure to pass a matrix containing one line per region to be highlighted with in each line : 
#' @details -the chromosome number 
#' @details -start 
#' @details -end
#' 
#' @seealso Fantasio
#' @seealso set.HFLOD
#' 
#' @return return a plot of the chromosome HBD segments for all individuals
#' 
#' @examples  
#' #Please refer to vignette 
#'
#' 
#' @export
HBD.plot.chr <- function(atlas, ROH, unit="cM", chr, list.ids, regions, outfile, inbred = FALSE, build=37)
{
  if(inbred)
  {
    list.ids <- which(atlas@submap_summary$inbred)
    ind      <- as.vector(atlas@submap_summary$id[list.ids])
    fam      <- as.vector(atlas@submap_summary$famid[list.ids])
    # list.ids <- paste(ind, fam, sep="_")
    list.ids <- uniqueIds(fam, ind)
  }
  
  if(class(atlas@bedmatrix)[1] != "bed.matrix")
    stop("Need a bed.matrix")
  
  if(is.null(atlas@HBD_recap))
    stop("HBD_recap is empty cannot plot, make sure to have atleast one individual considered inbred.")
  
  if(!missing(atlas) & !missing(ROH))
  {
    plot.ROH.segments.chr(ROH = ROH, submaps = atlas, unit = unit, chr = chr, outfile=outfile, listid=list.ids, regions=regions, build=build)
  }else{
    if(!missing(atlas))
      plot.HBD.segments.chr(atlas=atlas, unit=unit, chr=chr, list.ids=list.ids, regions=regions, outfile=outfile, build=build)
  }
}
