
#' @export
set.submap <- function(sx, segment.list, restore.seed = TRUE) {

  # storing current seed
  old.seed <- get(".Random.seed", envir = .GlobalEnv)

  # Note. Le CRAN n'acceptera pas ça. A régler plus tard 
  # (on doit pouvoir les bluffer en passant par Rcpp)
  # the object seed
  assign(".Random.seed", sx@random.seed, envir = .GlobalEnv)

  sx@submap <- rsubmap(segment.list)
   
  # restauring seed (meme commentaire)
  if(restore.seed)
    assign(".Random.seed", old.seed, envir = .GlobalEnv)

  sx
}
