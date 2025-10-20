#' Plot of the HFLOD 
#' 
#' This fonction plot the HFLOD score for a chromosome
#' 
#' @param submaps a atlas object
#' @param unit the unit used to plot, two options are allowed "bases", "cM" (default is "cM") 
#' @param chr the chromosome number from which to plot HFLOD score
#' @param regions a matrix containing the value to ve highlighted in the plot
#' @param color2 the color of the regions highlighted (default is "green4")
#' @param MA a boolean indicating whether a red line has to be drawn for the moving average
#' @param nbSNP_MA number of SNP for the moving average (default is 50)
#' 
#' @details If you use the regions options make sure to pass a matrix containing one line per region to be highlighted with in each line : 
#' @details - the chromosome number 
#' @details - start 
#' @details - end 
#' 
#' @seealso setHFLOD
#' 
#' @return This function returns a manhattan plot of all the HFLOD score over all the chromosome
#' 
#' 
#' @examples  
#' #Please refer to vignette 
#'
#' @export
HFLOD.plot.chr <- function(submaps, unit = c("cM", "bases"), chr, regions, color2="green4", MA = TRUE, nbSNP_MA = 50) 
{
  if(class(submaps@bedmatrix)[1] != "bed.matrix")
    stop("Need a bed.matrix.")
  
  if(is.null(submaps@HFLOD))
    stop("HFLOD slots in the object is empty, cannot plot")
  
  HFLOD <- submaps@HFLOD
  chromosome <- HFLOD$chr
  unit <- match.arg(unit)
  
  if(unit == "cM"){
    myxlab <- "Position (cM)"
    coeff  <- 1
  } else {
    myxlab <- "Position (Mb)"
    coeff  <- 1e6
  }
  
  if(missing(regions)) 
    myreg <- NULL
  else { 
    myreg  <- regions
    color2 <- "green4"
    myreg$start <- regions$start*coeff
    myreg$end   <- regions$end*coeff
  }

   #1)graphs per chromosome
  newout   <- NULL
  axis_mp  <- NULL
  chr_pos  <- 5 
  myreg_mp <- NULL
  
  toplot_HFLOD <- HFLOD$HFLOD[HFLOD$chr == chr]
  toplot_MA    <- HFLOD$HFLOD[HFLOD$chr == chr]

  if(unit == "cM")
    toplot_pos   <- HFLOD$dist[HFLOD$chr == chr]
  else
    toplot_pos   <- HFLOD$pos[HFLOD$chr == chr]
  
  ymax <- max(3.3,max(toplot_HFLOD))
  
  plot(toplot_pos,
       toplot_HFLOD,
       pch  = 16,
       ylim = c(0,ymax),
       xlab = myxlab,
       ylab = "HFLOD",
       cex.lab  = 1.4,
       cex.axis = 1.5,
       main     = paste("HFLOD (chromosome ",chr,")",sep=""),
       cex.main = 1.5)
  
  if(!(missing(regions))){
    myreg_chr <-  myreg[which(myreg$chr == chr),]
    if(nrow(myreg_chr) > 0){
      for(i in seq_len(nrow(myreg_chr))){
        polygon(x = myreg_chr[i,c("start","end","end","start")], 
                y = c(rep(-1,2),rep(ymax+1,2)), 
                col    = color2,
                border = color2,
                lwd    = 2)
        
        myreg_mp = rbind(myreg_mp,max(c(0,axis_mp))+10+myreg_chr$start[i]+myreg_chr$end[i])
      }
      points(toplot_pos, toplot_HFLOD, pch=16)
    }
  }
  for (i in seq_len(3)) 
    abline(h=i,col="grey",lwd=1,lty=2)
  
  abline(h=3.3,col="grey",lwd=2)
  
  if(MA)
    lines(toplot_pos, zoo::rollmean(toplot_HFLOD, as.numeric(nbSNP_MA), fill = "extend"), col="red", lwd=2)
  
  axis_mp <- c(axis_mp, max(c(0,axis_mp))+10+toplot_pos)
  chr_pos <- c(chr_pos, max(c(0,axis_mp))+5)
}
