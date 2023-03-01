
#' QQ-Plot & Manhattan Plots for glm on HBD prob or FLOD
#' 
#' @param x an atlas object
#' @param expl_var the explanatory variable 'FLOD' or 'pHBD'
#' @param n.cores number of cores for parallelization calculation (default = 1)
#' @param plot choose to plot 'all' = adj + unadj, 'adj' only, 'unadj' only (default = all)
#' if 'all' but only 'adj' or 'unadj' data are available the function will only plot available data
#' @param qq choose to plot qqplot (default = FALSE)
#' @param save choose if plot are saved or not (.png format) (default = FALSE)
#' 
#' @seealso \code{\link{glmHBD}}
#' 
#' @export

glm.HBD.plot = function ( x, expl_var, plot = c('all', 'unadj', 'adj'), qq = FALSE, save = FALSE) {
  
  plot <- match.arg(plot)
  
  if (plot == 'all') {
    if ('unadj' %in% names(x@logisticRegression)){
      unadj 	<- x@logisticRegression$unadj
      unadj <- unadj [which(unadj$p_value != 0), ]
    } else {
      stop('the operator $unadj is missing in atlas@logisticRegression - see documentation of `glmHBD` function to generate unadjusted data')}
    
    if ('adj' %in% names(x@logisticRegression)){
      adj <- x@logisticRegression$adj
      adj <- adj [which(adj$p_value != 0), ]
    } else {
      stop('the operator $adj is missing in atlas@logisticRegression - see documentation of `glmHBD` function to generate adjusted data')}
    
    # Manhattan Plot
    # Change colnames to fit with gaston
    colnames(unadj)[colnames(unadj) == 'pos_Bp'] <- 'pos'
    colnames(unadj)[colnames(unadj) == 'p_value'] <- 'p'
    
    colnames(adj)[colnames(adj) == 'pos_Bp'] <- 'pos'
    colnames(adj)[colnames(adj) == 'p_value'] <- 'p'
    
    
    treshold = -log10(0.05/sum(segments.list.summary(x@segments_list)$number_of_segments))
    lim <- round(max(-log10(adj$p), -log10(unadj$p), treshold))+1
    
    if (save == FALSE) { # Default, just print plots on different windows
      
      # QQ Plot
      if(qq) {
        gaston::qqplot.pvalues(unadj$p, main = paste("QQ-plot GLM with", expl_var, "- unadjusted data"), ylim = c(0,lim) )
        gaston::qqplot.pvalues(adj$p, main = paste("QQ-plot GLM with", expl_var, "- adjusted data"), ylim = c(0,lim))
      }
      
      #Manhattan plot
      gaston::manhattan(unadj, main = paste("Manhattan Plot \n GLM with ", expl_var, "- unadjusted data"), chrom.col = c("darksalmon", "darkturquoise"), ylim = c(0,lim))
      abline(h=treshold, col = 'red')
      
      gaston::manhattan(adj, main = paste("Manhattan Plot \n GLM with ", expl_var, "- adjusted data"), chrom.col = c("darksalmon", "darkturquoise"), ylim = c(0,lim))
      abline(h=treshold, col = 'red')
      
      
    } else { # save plots in png files 
      
      png(paste('qqplot.GLM.', expl_var,'unadj.png')) 		# U = unadjusted
      gaston::qqplot.pvalues(unadj$p, main = paste("QQ-plot GLM with", expl_var, "- unadjusted data"), ylim = c(0,lim) )
      dev.off()
      
      png(paste('qqplot.GLM.', expl_var, 'adj.png'))		# A = adjusted
      gaston::qqplot.pvalues(adj$p, main = paste("QQ-plot GLM with", expl_var, "- adjusted data"), ylim = c(0,lim))
      dev.off()
      
      png(paste( 'manhattanplot.GLM.', expl_var,'unadj.png'))
      gaston::manhattan(unadj, main = paste("Manhattan Plot \n GLM with ", expl_var, "- unadjusted data"), chrom.col = c("darksalmon", "darkturquoise"), ylim = c(0,lim))
      abline(h=treshold, col = 'red')
      dev.off()	
      
      png(paste( 'manhattanplot.GLM.', expl_var,'adj.png'))
      gaston::manhattan(adj, main = paste("Manhattan Plot \n GLM with ", expl_var, "- adjusted data"), chrom.col = c("darksalmon", "darkturquoise"), ylim = c(0,lim))
      abline(h=treshold, col = 'red')
      dev.off()	
      
    }
  }
  
  else if (plot == 'unadj') {
    if ('unadj' %in% names(x@logisticRegression)){
      unadj 	<- x@logisticRegression$unadj
      unadj <- unadj [which(unadj$p_value != 0), ]
    } else {
      stop('the operator $unadj is missing in atlas@logisticRegression - see documentation of `glmHBD` function to generate unadjusted data')}
    
    # Manhattan Plot
    # Change colnames to fit with gaston
    colnames(unadj)[colnames(unadj) == 'pos_Bp'] <- 'pos'
    colnames(unadj)[colnames(unadj) == 'p_value'] <- 'p'
    
    treshold = -log10(0.05/sum(segments.list.summary(x@segments_list)$number_of_segments))
    lim <- round(max( -log10(unadj$p), treshold))+1
    
    if (save == FALSE) { # Default, just print plots on different windows
      
      # QQ Plot
      if(qq) {
        gaston::qqplot.pvalues(unadj$p, main = paste("QQ-plot GLM with", expl_var, "- unadjusted data"), ylim = c(0,lim) )
      }
      
      #Manhattan plot
      gaston::manhattan(unadj, main = paste("Manhattan Plot \n GLM with ", expl_var, "- unadjusted data"), chrom.col = c("darksalmon", "darkturquoise"), ylim = c(0,lim))
      abline(h=treshold, col = 'red')
      
      
    } else { # save plots in png files 
      
      png(paste('qqplot.GLM.', expl_var,'unadj.png')) 		# U = unadjusted
      gaston::qqplot.pvalues(unadj$p, main = paste("QQ-plot GLM with", expl_var, "- unadjusted data"), ylim = c(0,lim) )
      dev.off()
      
      png(paste( 'manhattanplot.GLM.', expl_var,'unadj.png'))
      gaston::manhattan(unadj, main = paste("Manhattan Plot \n GLM with ", expl_var, "- unadjusted data"), chrom.col = c("darksalmon", "darkturquoise"), ylim = c(0,lim))
      abline(h=treshold, col = 'red')
      dev.off()	
      
    }
  }
  
  else if (plot == 'adj') {
    
    if ('adj' %in% names(x@logisticRegression)){
      adj <- x@logisticRegression$adj
      adj <- adj [which(adj$p_value != 0), ]
    } else {
      stop('the operator $adj is missing in atlas@logisticRegression - see documentation of `glmHBD` function to generate adjusted data')}
    
    # Manhattan Plot
    # Change colnames to fit with gaston
    colnames(adj)[colnames(adj) == 'pos_Bp'] <- 'pos'
    colnames(adj)[colnames(adj) == 'p_value'] <- 'p'
    
    treshold = -log10(0.05/sum(segments.list.summary(x@segments_list)$number_of_segments))
    lim <- round(max(-log10(adj$p), treshold))+1
    
    if (save == FALSE) { # Default, just print plots on different windows
      
      # QQ Plot
      if(qq) {
        gaston::qqplot.pvalues(adj$p, main = paste("QQ-plot GLM with", expl_var, "- adjusted data"), ylim = c(0,lim))
      }
      
      #Manhattan plot
      gaston::manhattan(adj, main = paste("Manhattan Plot \n GLM with ", expl_var, "- adjusted data"), chrom.col = c("darksalmon", "darkturquoise"), ylim = c(0,lim))
      abline(h=treshold, col = 'red')
      
      
    } else { # save plots in png files 
      
      png(paste('qqplot.GLM.', expl_var, 'adj.png'))		# A = adjusted
      gaston::qqplot.pvalues(adj$p, main = paste("QQ-plot GLM with", expl_var, "- adjusted data"), ylim = c(0,lim))
      dev.off()
      
      png(paste( 'manhattanplot.GLM.', expl_var,'adj.png'))
      gaston::manhattan(adj, main = paste("Manhattan Plot \n GLM with ", expl_var, "- adjusted data"), chrom.col = c("darksalmon", "darkturquoise"), ylim = c(0,lim))
      abline(h=treshold, col = 'red')
      dev.off()	
      
    }
  }
}
