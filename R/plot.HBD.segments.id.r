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

plot.HBD.segments.id <- function(Submaps, unit= "cM", id, famid, regions, outfile, build)
{
  if(!is.character(id))
    return("Need individual id as character")
  if(!is.character(famid))
    return("Need family id as character")
  
  HBD.recap <- Submaps@HBD_recap
  HBDsegments <- Submaps@HBDsegments
  
  HBDsegments_rbind <- do.call(rbind, HBDsegments) #binding lines 
  
  HBD <- HBDsegments_rbind[which(HBDsegments_rbind$id==id & HBDsegments_rbind$famid==famid),]
  
  if(nrow(HBD) == 0)
    stop("No individual found")
  
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
  
  plot.segments.id(fileOrSubmaps=HBD, unit = unit, regions = myreg, main=paste("HBDsegments of", uniqueIds(famid, id)), build=build)
}
