length.chromosome <- function (chr, units = c("bases", "cM"), build = 37) {
  units <- match.arg(units)
  if(!(build %in% 35:38))
    stop("Argument 'build' should be 35, 36, 37 or 38")
  
  x <- if(build == 35) Fantasio2::cytobands.b35
  else if(build == 36) Fantasio2::cytobands.b36
  else if(build == 37) Fantasio2::cytobands.b37
  else if(build == 38) Fantasio2::cytobands.b38
  
  x$end <- if(units == "bases") x$end.bp else x$end.cM
  
  tapply(x$end, x$chr, function(a) tail(a,1) )[ as.character(chr) ]
}