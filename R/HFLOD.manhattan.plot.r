#' Creation of an manhattan plot of the HFLOD
#'
#' This fonction to plot a manhanttan plot of the HFLOD score
#'
#' @param submaps a atlas object
#' @param regions a matrix containing the value to be highlighted in the plot
#' @param unit the unit used to plot, two options are allowed "Bases", "cM" (default is "CM")
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
HFLOD.manhattan.plot <- function(submaps, regions, unit = "cM", MA = TRUE, nbSNP_MA = 50)
{
  if (class(submaps@bedmatrix)[1] != "bed.matrix")
    stop("Need a bed.matrix.")
  
  if (is.null(submaps@HFLOD))
    stop("HFLOD slots in the object is empty, cannot plot")
  
  HFLOD <- submaps@HFLOD
  #to get mean position when working by segments
  if (unit == "cM")
    pos <- HFLOD$dist
  else
    pos <- HFLOD$pos
  
  chromosome <- HFLOD$chr
  
  
  if (missing(regions))
    myreg <- NULL
  else {
    myreg  <- regions
    color2 <- "green4"
    myreg$start = regions$start / 1e6
    myreg$end   = regions$end / 1e6
  }
  
  if (unit == "cM") {
    myxlab <- "Position (cM)"
    coeff  <- 1
  } else{
    myxlab <- "Position (Mb)"
    coeff  <- 1e6
    pos    <- pos / 1e6
  }
  
  
  
  newout   <- NULL
  axis_mp  <- NULL
  chr_pos  <- numeric(length(unique(chromosome)+1)); chr_pos[1] <- 5
  myreg_mp <- NULL
  
  if (!missing(regions)) {
    myreg_chr <-  myreg[which(myreg$chr == c), ]
    if (nrow(myreg_chr) > 0) {
      for (i in seq_len(nrow(myreg_chr))) {
        polygon(
          x = myreg_chr[i, c("start", "end", "end", "start")],
          y = c(rep(-1, 2), rep(ymax + 1, 2)),
          col    = color2,
          border = color2,
          lwd    = 2
        )
        
        myreg_mp = rbind(myreg_mp,
                         max(c(0, axis_mp)) + 10 + myreg_chr$start[i] + myreg_chr$end[i])
      }
    }
  }
  
  #2)Manhattan plot
  
  
  ymax <- max(3.3, max(HFLOD$HFLOD))
  mycol <- rep(c("cadetblue2", 8), 11)
  
  
  for (i in unique(chromosome))
  {
    pos_chr <- pos[chromosome == i]
    chr_pos[i+1] <- pos_chr[length(pos_chr)]
    axis_mp <- c(axis_mp, max(c(0, axis_mp), na.rm = TRUE) + 10 + pos_chr)
  }
  chr_pos  <- cumsum(chr_pos + 10)
  chr_axis <- sapply(seq_along(chr_pos)-1, function(i) mean(c(chr_pos[i], chr_pos[i+1])))
  chr_axis <- chr_axis[-1]
  
  plot (
    axis_mp,
    HFLOD$HFLOD,
    pch = 16,
    ylim = c(0, ymax),
    xlab = "",
    ylab = "HFLOD",
    cex.lab = 1.4,
    cex.axis = 1.5,
    col = mycol[chromosome],
    xaxt = "n",
    cex = 0.75
  )
  
  if (!missing(regions)) {
    for (i in seq_len(nrow(myreg_mp))) {
      polygon(
        x = myreg_mp[i, c("start", "end", "end", "start")] / coeff,
        y = c(rep(-10, 2), rep(max(HFLOD$HFLOD) + 10, 2)),
        col = color2,
        border = color2,
        lwd = 2
      )
    }
    points(axis_mp,
           HFLOD[, 1],
           pch = 16,
           col = mycol[chromosome],
           cex = 0.75)
  }
  
  if(MA)
    lines(axis_mp, zoo::rollmean(HFLOD$HFLOD, as.numeric(nbSNP_MA), fill = "extend"), col="red", lwd=2)
  
  for (i in seq_along(unique(chromosome))) {
    abline(v = chr_pos[i], col = "grey", lwd = 2)
    axis(1, at = chr_axis[i], i, col.ticks = 0, cex.axis = 1.5)
  }
  
  abline(v = chr_pos[ length(chr_pos) ], col = "grey", lwd = 2)
  
  for (i in seq_len(3))
    abline(
      h = i,
      col = "grey",
      lwd = 1,
      lty = 2
    )
  
  
  abline(h = 3.3, col = "grey", lwd = 2)
  
  
  
}
