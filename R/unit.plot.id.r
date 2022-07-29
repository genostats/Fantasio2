unit.plot.id <- function(file, unit, byROHfile)
{
  if(byROHfile) {
    if (unit=="cM"){
      ecart <- 25
      larg  <- 2
      pos1  <- which(colnames(file)=="POS1_cM")
      pos2  <- which(colnames(file)=="POS2_cM")
      if (length(pos1)==0) {
        stop("no genetic distances in ROH file")
      }
    } else if (unit=="bases") {
      ecart <- 25e6
      larg  <- 2
      pos1  <- which(colnames(file)=="POS1")
      pos2  <- which(colnames(file)=="POS2")
    } else
      stop("units parameter accepts cM and bases only")
  } else { 
    if(unit=="cM"){
      ecart <- 25
      larg  <- 2
      pos1  <- which(colnames(file)=="start_dist")
      pos2  <- which(colnames(file)=="end_dist")
      if (length(pos1)==0) {
        stop("no genetic distances for this individual")
      }
    } else if(unit=="bases"){
      ecart <- 2.5e7
      larg  <- 2
      pos1  <- which(colnames(file)=="start_pos")
      pos2  <- which(colnames(file)=="end_pos")
    } else 
      stop("length option accepts cM and bases only")
  }
  
  list("ecart" = ecart, "larg"  = larg, "pos1"  = pos1, "pos2"  = pos2)
}