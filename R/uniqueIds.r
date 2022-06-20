uniqueIds <- function(famid, id) paste(famid, id, sep = ":") 

get.famid <- function(x) {
  if(is.null(x)) 
    NULL
  else
    sapply( strsplit(x, ":"), function(x) x[1] )
}

get.id <- function(x) {
  if(is.null(x)) 
    NULL
  else
    sapply( strsplit(x, ":"), function(x) x[2] )
}
