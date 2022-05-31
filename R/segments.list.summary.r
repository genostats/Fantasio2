segments.list.summary <- function(segment.list)
{
  if(class(segment.list)[1] == "snpsSegments")
    segment.list <- segment.list@snpsSegments
  else if(class(segment.list)[1] != "HostspotsSegments" & class(segment.list)[1] != "list")
    stop("Argument must be of class 'snpsSegments', 'HostspotsSegments' or 'list'")
  
  #number of segments
  n_seg <- integer(length(segment.list))
  
  for(i in seq_along(segment.list))
    n_seg[i] <- length(segment.list[[i]])
  
  #number of markers 
  n_mark <- integer(length(segment.list))
  for(i in seq_along(segment.list)) {
    res <- integer(length(segment.list[[i]]))
    for(j in seq_along(segment.list[[i]]))
    {
      if (is.list(segment.list[[i]][[j]])) {
        l <- sapply(segment.list[[i]][[j]], function(k) length(k[1]: k[2]) )
        res[j] <- sum(l)
      } else {
        if (length(segment.list[[i]][[j]]) == 0L) next
        else if (length(segment.list[[i]][[j]]) == 1L) { 
          res[j] <- 1L
        } else { res[j] <- length(segment.list[[i]][[j]][1]:segment.list[[i]][[j]][2])
        }
      }
    }
    res <- sum(res)
    n_mark[i] <- res
  }
  
  
  #dataframe
  df <- data.frame(
    chromosome = if(class(segment.list)[1] == "HostspotsSegments") seq_along(segment.list) else getOption("gaston.autosomes"),
    number_of_segments = n_seg, 
    number_of_markers= n_mark
  )
  df
}
