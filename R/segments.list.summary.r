segments.list.summary <- function(segments.list)
{
  if( class(segments.list)[1] != "snpsSegments" & class(segments.list)[1] != "HostspotsSegments" & class(segments.list)[1] != "list")
    stop("Argument must be of class 'snpsSegments', 'HostspotsSegments' or 'list'")
  
  chr <- names(segments.list)
  n.seg <- sapply(segments.list, function(x) length(x$beg))
  n.mrk <- sapply(segments.list, function(x) sum(x$end - x$beg + 1L))

  data.frame(chromosome = chr, number_of_segments = n.seg, number_of_markers = n.mrk)
}
