#' Manhattan Plot for glm on HBD prob or FLOD
#' 
#' @param x a data frame such as sent by HBD.glm
#' @param test which test to plot
#' @param chrom.col chromosome colors
#' @param ... extra arguments for the plot
#' 
#' @seealso \code{\link{HBD.glm}}
#' 
#' @export

glm.HBD.plot = function (x, test = c("bilateral", "right", "left"), chrom.col = c("darksalmon", "darkturquoise"), ... ) {
  
  test <- match.arg(test)
  p.name <- paste0("p.", test)
  x$p <- x[[p.name]]

  # Manhattan Plot
  # Change colnames to fit with gaston
  # colnames(x)[colnames(x) == 'pos_Bp'] <- 'pos'

  gaston::manhattan(x, chrom.col = chrom.col , ...)
}
