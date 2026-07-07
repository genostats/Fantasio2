#'Find HBD segments
#' 
#'This function creates HBD segments used for creating plots                      
#'
#' @param x a atlas object
#' @param n.consecutive.markers the number of consecutive markers with a probabilitie equal or greater to the value of the threshold, to be use to fing HBD segments (default is 5)
#' @param threshold the minimum value of HBD probabilities given for a marker (default is 0.5)
#' 
#' @details The threshold is the minimum value from which we consider the marker is HBD. From this marker we want a minumum number of consecutive markers to create a HBD segment.
#' (argument : n.consecutive.markers)
#' 
#' @return This function returns a list of dataframe with 11 columns : 
#' @return - id
#' @return - famid
#' @return - pheno
#' @return - start
#' @return - end
#' @return - size
#' @return - chromosome
#' @return - start_pos
#' @return - end_pos
#' @return - start_dist
#' @return - end_dist
#' 
#' @keywords internal
HBD.segments <- function(x, n.consecutive.markers = 5, threshold = 0.5) {
  HBD_recap <- x@HBD_recap
  L <- list()
  
  # tester s'il y a des inds consanguins
  # si non, return NULL
  if(is.null(HBD_recap)) {
    return(NULL)
  }

  un.ids <- rownames(HBD_recap) # famid:id
  status <- x@bedmatrix@ped$pheno[ match( un.ids, uniqueIds(x@bedmatrix@ped$famid, x@bedmatrix@ped$id) ) ]
  individuals_name <- get.id(un.ids)
  family_id <- get.famid(un.ids)

  # ### fin modifs ####

  marker <- colnames(HBD_recap)#save the marker

  correspondance <- match(colnames(HBD_recap), x@bedmatrix@snps$id)#match between marker's name in HBD_recap and the bedmatrix
  chr <- x@bedmatrix@snps$chr[correspondance]                      #chromosome on which we have the marker
  min_segment_size <- n.consecutive.markers                               #minimum size of the marker

  for(i in seq_len(nrow(HBD_recap))) {
    data <- as.vector(HBD_recap[i,])       #save the line 
    test <- (data >= threshold) #test
  
    # get the segments
    all_segments <- rle(test)
    good_segments <- as.numeric(which( all_segments$length >= min_segment_size & all_segments$value )-1) #first marker
  
    if(length(good_segments) == 0) next
  
    if(good_segments[1] == 0) good_segments[1] <- 1
    
    good_segments_length <- as.numeric(all_segments$length[ good_segments+1 ]) 
    good_segments_start <- as.numeric(cumsum(all_segments$lengths)[ good_segments ]+1)#segment start
    good_segments_end <- as.numeric(good_segments_start+good_segments_length-1)
  
    # finding distance and position 
  
    start_pos <- x@bedmatrix@snps$pos[correspondance[as.numeric(good_segments_start)]]
    end_pos <- x@bedmatrix@snps$pos[correspondance[as.numeric(good_segments_end)]]
  
    start_dist <- x@bedmatrix@snps$dist[correspondance[as.numeric(good_segments_start)]]
    end_dist <- x@bedmatrix@snps$dist[correspondance[as.numeric(good_segments_end)]]
  
    # dataframe
  
    segment_dataframe <- data.frame(id         = rep(individuals_name[i], length(start_pos)),
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
  
    # treating the case when segments overlaps two different chromosomes
    overlap <- which(segment_dataframe$start_dist > segment_dataframe$end_dist)
  
    if(length(overlap) != 0) {
      segment_dataframe <- segment_dataframe[-overlap,,drop=F]
    }
  
    L[[i]] <- segment_dataframe
  }
  L
}
