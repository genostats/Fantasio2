##################################################################################
#This function creates a dataframe of HBD probabilities mean by snps             #
#                                                                                #
#!!! submaps : the list of object                                                #                                       
#!!! HBD_recap : the HBD_recap dataframe                                         #
#!!! n.consecutive.markers : the number of consecutives markers with an HBD >      #
#    threshold                                                                   #
#!!! threshold : a threshold for the number of consecutive markers HBD            #
#                                                                                #
#*** return a list of dataframe with HBD segment by snps                         #
##################################################################################

HBD.segments.by.snps <- function(submaps, HBD_recap, n.consecutive.markers, threshold)
{ 
  l <- list()
  
  un.ids <- rownames(HBD_recap) # famid:id
  status <- submaps@bedmatrix@ped$pheno[ match( un.ids, uniqueIds(submaps@bedmatrix@ped$famid, submaps@bedmatrix@ped$id) ) ]
  individuals_name <- get.id(un.ids)
  family_id <- get.famid(un.ids)
  
  # ### fin modifs ####
  
  marker <- colnames(HBD_recap)#save the marker
  
  correspondance <- match(colnames(HBD_recap), submaps@bedmatrix@snps$id)#match between marker's name in HBD_recap and the bedmatrix
  chr <- submaps@bedmatrix@snps$chr[correspondance]                      #chromosome on which we have the marker
  min_segment_size <- n.consecutive.markers                               #minimum size of the marker
  
  for(i in seq_len(nrow(HBD_recap)))
  {
    data<-c(HBD_recap[i,])#save the line 
    test<- (data >= threshold)#test
    
    
    # get the segments
    
    all_segments<-rle(test)
    
    good_segments<- as.numeric(which( all_segments$length >= min_segment_size & all_segments$value )-1) #first marker
    
    if(length(good_segments) == 0)
      next()
    
    if(good_segments[1] == 0)
      good_segments[1] <- 1
    
    good_segments_length <-as.numeric(all_segments$length[ good_segments+1 ]) 
    
    good_segments_start<- as.numeric(cumsum(all_segments$lengths)[ good_segments ]+1)#segment start
    
    good_segments_end<-as.numeric(good_segments_start+good_segments_length-1)
    
    
    
    
    
    #finding distance and position 
    
    start_pos <- submaps@bedmatrix@snps$pos[correspondance[as.numeric(good_segments_start)]]
    end_pos <- submaps@bedmatrix@snps$pos[correspondance[as.numeric(good_segments_end)]]
    
    
    start_dist <- submaps@bedmatrix@snps$dist[correspondance[as.numeric(good_segments_start)]]
    end_dist <- submaps@bedmatrix@snps$dist[correspondance[as.numeric(good_segments_end)]]
    
    
    #dataframe
    
    segment_dataframe<-data.frame(id         = rep(individuals_name[i], length(start_pos)),
                                  famid      = rep(family_id[i], length(start_pos)),
                                  pheno      = rep(status[i], length(start_pos)),
                                  start      = good_segments_start, 
                                  end        = good_segments_end ,
                                  size       = good_segments_length,
                                  chromosome = chr[as.numeric(good_segments_start)],
                                  start_pos  = start_pos,
                                  end_pos    = end_pos,
                                  start_dist = start_dist,
                                  end_dist   = end_dist)
    
    #treating the case when segments overlaps two different chromosomes
    overlap <- which(segment_dataframe$start_dist > segment_dataframe$end_dist)
    
    
    if(length(overlap) != 0)
    {
      segment_dataframe <- segment_dataframe[-overlap,,drop=F]
    }
    
    l[[i]] <- segment_dataframe
  }
  return(l)
}