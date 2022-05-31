clean.hotspots <- function(l, number_of_marker)
{
  index <- c()
  for(i in seq_along(l)) {
    if(length(l[[i]]) < number_of_marker) {
      if(i == length(l))
        l[[i-1]] <- c(l[[i-1]], l[[i]])
      else
        l[[i+1]] <- c(l[[i]], l[[i+1]])
    } else {
      index <- c(index, i)
    }
  }
  l <- l[index]
  l
}