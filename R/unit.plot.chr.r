unit.plot.chr <- function(file, unit, byROHfile)
{
  if(byROHfile)
  {
    if (unit=="cM"){
      pos1 <- which(colnames(file)=="POS1_cM")
      pos2 <- which(colnames(file)=="POS2_cM")
      myxlab <- "Position (cM)"
      coeff=1
      if (length(pos1)==0)
        stop("no genetic distances in ROH file")
    }
    else if(unit=="bases"){
      pos1 <- which(colnames(file)=="POS1")
      pos2 <- which(colnames(file)=="POS2")
      myxlab <- "Position (Mb)"
      coeff <- 1e6
    }
    else {
      stop("units parameter accepts 'cM' and 'bases' only")
    }
  }else{
    #choose which unit will be use for plot
    
    if (unit=="cM") {
      pos1 <- which(colnames(file)=="start_dist"); 
      pos2 <- which(colnames(file)=="end_dist");
      myxlab="Position (cM)"; 
      coeff <- 1; 
      if (length(pos1)==0) 
        stop("no genetic distances in HBD file")
    }
    else if (unit=="bases") {
      pos1 <- which(colnames(file)=="start_pos"); 
      pos2 <- which(colnames(file)=="end_pos"); 
      myxlab <- "Position (Mb)"; 
      coeff <- 1e6 # pour des unites en Mb
    }
    else {
      stop("units parameter accepts 'cM' and 'bases' only")
    }
  }
  
  returnList <- list("pos1"   = pos1,
                     "pos2"   = pos2, 
                     "myxlab" = myxlab, 
                     "coeff"  = coeff)
  
  return(returnList)
}