##################################################################################
#This function is uses to create a plot of the HBD for a chromosome and for all  #
#the individual on the submap                                                    #
#                                                                                #
#!!! Submaps : the list of object                                                #                                       
#!!! unit : bases or cM                                                          #
#!!! chr : the number of the chromosome wanted                                   #
#!!! list.ids : (optional) a list of indivial to plot with                       #
#!!! regions  : a regions to be emphasize                                        #
#!!! outfile : the name of the plot                                              #
#                                                                                #
#*** return a new submap object                                                  #
##################################################################################

plot.HBD.segments.chr <- function(Submaps, unit, chr, list.ids, regions, outfile, build)
{
  HBD_segments <- Submaps@HBD_segments
  
  if(missing(regions)) 
    myreg <- NULL
  else { 
    myreg <- regions[regions$chr == chr,]
  }
  
  #Get the lines for the wanted chromosome in the HBD_segments dataframe
  HBD_segments_rbind <- do.call(rbind, HBD_segments) #binding lines
  
  HBD <- subset(HBD_segments_rbind, HBD_segments_rbind$chromosome==chr)#only the wanted lines
  
  HBD$id <- as.character(HBD$id) #otherwise factors level in the vector
  HBD$famid <- as.character(HBD$famid)
  
  if (missing(list.ids)) {
    list.ids <- unique(paste(HBD$famid, HBD$id, sep = ":"))  
  }
  
  #name the file
  if (missing(outfile) )
    outfile <- paste("HBD_chr_",chr,"_",unit,".png",sep="") 
  else 
    outfile <- paste(outfile,".png",sep="") 
  
  plot.segments.chr(fileOrSubmaps=HBD,unit= unit,chr= chr,list_id = list.ids,regions = myreg, build=build)
}
