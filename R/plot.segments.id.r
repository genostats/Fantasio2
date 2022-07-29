plot.segments.id <- function (byROHfile = FALSE, fileOrSubmaps, unit = "cM",
                            regions, color2 = "green4", main, build=37) {
  # chose which unit to use for the plot
  l <- unit.plot.id(file = fileOrSubmaps, unit = unit, byROHfile)
  ecart <- l$ecart
  larg  <- l$larg
  pos1  <- l$pos1
  pos2  <- l$pos2
  
  #empty plot
  if (missing(main))
    main = NULL
  plot( x = c(larg * 0.5, larg * 11.5),
        y = c( 0,
               length.chromosome(1, unit, build) + lengthChromosome(12, unit, build) + 1.25 * ecart),
        type = "n", xaxt = "n", yaxt = "n", xlab = "", ylab = "", main = main, cex.main = 1.5)
  
  
  #chromosome loop
  for (i in seq_len(22)) {
    if (i < 12) {
      i2 <- i
      offset_y <- length.chromosome(12, unit, build) + ecart
    } else {
      i2 <- i - 11
      offset_y <- 0
    }
    
    #dessiner le chromosome
    paint.cytobands( i, units = unit,
                    pos = c(larg * i2, offset_y + length.chromosome(i, unit, build) + ecart / 4),
                    orientation = "v", legend = FALSE, build = build)
    
    text(larg * i2, offset_y - 0.25 * ecart, i)
    
    #traitement de l'options regions
    if (!is.null(regions)) {
      regions_chr <- regions[regions$chr == i,]
      if (nrow(regions_chr) > 0) {
        for (k in seq_len(nrow(regions_chr))) {
          xx <- c(rep(larg * i2, 2), rep(larg * (i2 + 0.5), 2))
          
          yy <- c( 0.25 * ecart + abs( length.chromosome(i, unit, build) - c(regions_chr$start[k], regions_chr$end[k])),
                   0.25 * ecart + abs( length.chromosome(i, unit, build) - c(regions_chr$end[k], regions_chr$start[k])))
          
          polygon( x = xx, y = yy + offset_y, col = color2, border = color2, lwd = 1.5)
        }
      }
    }
    
    #recuperer les lignes pour un individus et un chromosome specifique
    if (byROHfile) {
      seg_chr <- fileOrSubmaps[fileOrSubmaps$CHR == i, ]
    } else {
      seg_chr <- fileOrSubmaps[fileOrSubmaps$chromosome == i, ]
    }
    
    if (nrow(seg_chr) > 0) {
      for (k in seq_len(nrow(seg_chr))) {
        xx <- c(rep(larg * i2, 2), rep(larg * (i2 + 0.5), 2))
        
        yy <- c( 0.25 * ecart + abs( length.chromosome(i, unit, build) - c(seg_chr[k, pos1], seg_chr[k, pos2])),
                 0.25 * ecart + abs( length.chromosome(i, unit, build) - c(seg_chr[k, pos2], seg_chr[k, pos1])))
        
        polygon( x = xx, y = yy + offset_y, col = ifelse(byROHfile, color(seg_chr$PHE[k]), color(seg_chr$pheno[k])),
                 border = ifelse( byROHfile, color(seg_chr$PHE[k]), color(seg_chr$pheno[k])),
                 lwd    = 1)
      }
    }
  }
}