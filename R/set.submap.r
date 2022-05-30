
#' @export
set.submap <- function(sx, ...) {

  restore.seed <- getArg("restore.seed", TRUE, ...)

  # storing current seed
  old.seed <- get(".Random.seed", envir = .GlobalEnv)

  # Note. Le CRAN n'acceptera pas ça. A régler plus tard 
  # (on doit pouvoir les bluffer en passant par Rcpp)
  # the object seed
  assign(".Random.seed", sx@random.seed, envir = .GlobalEnv)

  # cette partie est à remplacer par quelque chose d'intelligent
  n <- getArg("n", 100, ...)
  submap <- sort(sample.int(ncol(sx), n))
  sx@submap <- submap
   
  # restauring seed (meme commentaire)
  if(restore.seed)
    assign(".Random.seed", old.seed, envir = .GlobalEnv)

  sx
}
