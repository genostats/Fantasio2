segments.list.summary <- function(segmentList)
{
  if(class(segmentList)[1] == "snpsSegments")
    segmentList <- segmentList@snpsSegments
  else if(class(segmentList)[1] != "HostspotsSegments" & class(segmentList)[1] != "list")
    stop("Argument must be of class 'snpsSegments', 'HostspotsSegments' or 'list'")
  
  #number of segments
  n_seg <- numeric(length(segmentList))
  
  for(i in seq_along(segmentList))
    n_seg[i] <- length(segmentList[[i]])
  
  #number of markers 
  n_mark <- numeric(length(segmentList))
  for(i in seq_along(segmentList))
  {
    res <- numeric(length(segmentList[[i]]))
    for(j in seq_along(segmentList[[i]]))
    {
      if (is.list(segmentList[[i]][[j]])) {
        l <- sapply(segmentList[[i]][[j]], function(k) length(k[1]: k[2]) )
        res[j] <- sum(l)
      } else {
        if (length(segmentList[[i]][[j]]) == 0) next
        else if (length(segmentList[[i]][[j]]) == 1) { 
          res[j] <- 1 
        } else { res[j] <- length(segmentList[[i]][[j]][1]:segmentList[[i]][[j]][2])
        }
      }
    }
    res <- sum(res)
    n_mark[i] <- res
  }
  
  
  #dataframe
  df <- data.frame(
    chromosome = if(class(segmentList)[1] == "HostspotsSegments") seq_along(segmentList) else getOption("gaston.autosomes"),
    number_of_segments = n_seg, 
    number_of_markers= n_mark
  )
  df
}