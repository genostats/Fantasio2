paint.cytobands <- function (chr, pos = c(0, 0), units = c("bases", "cM"), build = 37,
                            width = 0.4, length.out, orientation = c("h", "v"), legend = TRUE, 
                            cex.leg = 0.7, ...) {
  
  units <- match.arg(units)
  if(!(build %in% 35:38))
    stop("Argument 'build' should be 35, 36, 37 or 38")
  
  orientation <- match.arg(orientation)
  
  x <- if(build == 35) Fantasio2::cytobands.b35
  else if(build == 36) Fantasio2::cytobands.b36
  else if(build == 37) Fantasio2::cytobands.b37
  else if(build == 38) Fantasio2::cytobands.b38
  
  x <- x[x$chr == chr, ]
  n <- nrow(x)
  if(n == 0)
    stop("No data for chr = ", chr)
  
  x$end <- if(units == "bases") x$end.bp else x$end.cM
  x$beg <- if(units == "bases") x$beg.bp else x$beg.cM
  if(!missing(length.out)) {
    M <- max(x$end)
    x$beg <- x$beg / M * length.out
    x$end <- x$end / M * length.out
  }
  
  band <- as.integer(x$stain)
  band.col <- gray( c(0.5, 1, 0.2, 0.6, 0.75) )[band]
  band.dens <- c(30, -1, -1, -1, 10)[band]
  band.bord <- gray(c(0, 0, 0, 0, 1))[band]
  
  ctr <- max(which( x$arm == "p"))
  idx <- c(2:(ctr - 1), (ctr + 2):(n - 1))
  
  t <- seq(-pi/2, pi/2, length = 25)
  if(orientation == "h") {
    rect(pos[1] + x$beg[idx], pos[2], pos[1] + x$end[idx], pos[2] - width, col = band.col[idx], density = band.dens[idx], border = band.bord[idx])
    # left end
    ra <- x$end[1] - x$beg[1]
    polygon( pos[1] + x$beg[2] - ra*cos(t), pos[2] - width/2 + width/2*sin(t), col = band.col[1], density = band.dens[1], border = band.bord[1])
    # right end
    ra <- x$end[n] - x$beg[n]
    polygon( pos[1] + x$end[n-1] + ra*cos(t), pos[2] - width/2 + width/2*sin(t), col = band.col[n], density = band.dens[n], border = band.bord[n])
    # fin du bras p
    ra <- x$end[ctr] - x$beg[ctr]
    polygon( pos[1] + x$end[ctr-1] + ra*cos(t), pos[2] - width/2 + width/2*sin(t), col = band.col[ctr], density = band.dens[ctr], border = band.bord[ctr])
    # début du bras q
    ra <- x$end[ctr+1] - x$beg[ctr+1]
    polygon( pos[1] + x$beg[ctr+2] - ra*cos(t), pos[2] - width/2 + width/2*sin(t), col = band.col[ctr+1], density = band.dens[ctr+1], border = band.bord[ctr+1])
    if(legend) 
      text( pos[1] + (x$beg + x$end)/2, pos[2] + width/2, paste0(x$arm, x$band), adj = c(0, .5), srt = 90, cex = cex.leg)
  } else {
    rect(pos[1], pos[2] - x$beg[idx], pos[1] - width, pos[2] - x$end[idx], col = band.col[idx], density = band.dens[idx], border = band.bord[idx])
    # left end (au-dessus sur le dessin)
    ra <- x$end[1] - x$beg[1]
    polygon( pos[1] - width/2 + width/2*sin(t), pos[2] - x$beg[2] + ra*cos(t), col = band.col[1], density = band.dens[1], border = band.bord[1])
    # right end
    ra <- x$end[n] - x$beg[n]
    polygon( pos[1] - width/2 + width/2*sin(t), pos[2] - x$end[n-1] - ra*cos(t), col = band.col[n], density = band.dens[n], border = band.bord[n])
    # fin du bras p
    ra <- x$end[ctr] - x$beg[ctr]
    polygon( pos[1] - width/2 + width/2*sin(t), pos[2] - x$end[ctr-1] - ra*cos(t), col = band.col[ctr], density = band.dens[ctr], border = band.bord[ctr])
    # début du bras q
    ra <- x$end[ctr+1] - x$beg[ctr+1]
    polygon( pos[1] - width/2 + width/2*sin(t), pos[2] - x$beg[ctr+2] + ra*cos(t), col = band.col[ctr+1], density = band.dens[ctr+1], border = band.bord[ctr+1])
    if(legend) 
      text( pos[1] + width/2, pos[2] - (x$beg + x$end)/2, paste0(x$arm, x$band), adj = c(0, .5), srt = 0, cex = cex.leg)
  }
}