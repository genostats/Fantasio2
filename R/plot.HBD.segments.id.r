##################################################################################
#This function plots the HBD probabilities by segments                           #
#                                                                                #
#!!! Submaps : the list of object                                                #                                       
#!!! unit : cM or Bases                                                          #
#!!! regions : a region to be emphasize in the plot                              #
#!!! outfile: (optional) a name for the plot                                     #
#!!! famid : the family id                                                   #
#!!! id = the individual id                                           #
#                                                                                #
#*** return a plot                                                               #
##################################################################################

plot.HBD.segments.id <- function(Submaps, unit= "cM", id, famid, regions, quality = 95, outfile, build)
{
  if(!is.character(id))
    return("Need individual id as character")
  if(!is.character(famid))
    return("Need family id as character")
  
  HBD.recap <- Submaps@HBD_recap
  HBD_segments <- Submaps@HBD_segments
  
  HBD_segments_rbind <- do.call(rbind, HBD_segments) #binding lines 
  
  HBD <- HBD_segments_rbind[which(HBD_segments_rbind$id==id & HBD_segments_rbind$famid==famid),]
  
  if(nrow(HBD) == 0)
    if ((id %in% Submaps@submap_summary$id & famid %in% Submaps@submap_summary$famid ) == FALSE)
      stop("No individual found, check spelling of id and famid")
    else if (Submaps@submap_summary$quality[which(Submaps@submap_summary$id == id & Submaps@submap_summary$famid == famid)]<= quality)
      stop("No HBD segment evaluation for this individual because of low QUALITY (<=",quality,")")
    else if(Submaps@submap_summary$inbred[which(Submaps@submap_summary$id == id & Submaps@submap_summary$famid == famid)] == FALSE)
      stop("Individual with good QUALITY (>=",quality,") but No HBD segment evaluation because he is not inbred")
  
  #regions options
  if (missing(regions)) 
    myreg <- NULL
  else
    myreg <- regions
  
  #name the file
  if (missing(outfile)) 
    outfile <- paste("HBD_", id,"_",unit,".png",sep="")
  else {
    outfile <- paste(outfile,".png",sep="") 
  }
  
  plot.segments.id(fileOrSubmaps=HBD, unit = unit, regions = myreg, main=paste("HBD segments of", uniqueIds(famid, id)), build=build)
}
