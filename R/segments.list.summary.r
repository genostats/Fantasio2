segments.list.summary <- function(segments.list)
{
  if(class(segments.list)[1] == "snpsSegments")
    segments.list <- segments.list@snpsSegments
  else if(class(segments.list)[1] != "HostspotsSegments" & class(segments.list)[1] != "list")
    stop("Argument must be of class 'snpsSegments', 'HostspotsSegments' or 'list'")
  
  #number of segments
  n_seg <- integer(length(segments.list))
  
  for(i in seq_along(segments.list))
    n_seg[i] <- length(segments.list[[i]])
  
  #number of markers 
  n_mark <- integer(length(segments.list))
  for(i in seq_along(segments.list)) {
    res <- integer(length(segments.list[[i]]))
    for(j in seq_along(segments.list[[i]]))
    {
      if (is.list(segments.list[[i]][[j]])) {
        l <- sapply(segments.list[[i]][[j]], function(k) length(k[1]: k[2]) )
        res[j] <- sum(l)
      } else {
        if (length(segments.list[[i]][[j]]) == 0L) next
        else if (length(segments.list[[i]][[j]]) == 1L) { 
          res[j] <- 1L
        } else { res[j] <- length(segments.list[[i]][[j]][1]:segments.list[[i]][[j]][2])
        }
      }
    }
    res <- sum(res)
    n_mark[i] <- res
  }
  
  
  #dataframe
  df <- data.frame(
    chromosome = if(class(segments.list)[1] == "HostspotsSegments") seq_along(segments.list) else getOption("gaston.autosomes"),
    number_of_segments = n_seg, 
    number_of_markers= n_mark
  )
  df
}
