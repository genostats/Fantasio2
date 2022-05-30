
#' @exportClass sub.bed.matrix
setClass("sub.bed.matrix", contains = "bed.matrix", slots = list(random.seed = "integer", submap = "integer"))

init.sub.bed.matrix <- function(.Object, ...) {
  .Object <- callNextMethod() 
  # step RNG (to avoid using twice the same seed
  runif(1)
  .Object@random.seed <- get(".Random.seed", envir = .GlobalEnv)
  .Object
}

setMethod("initialize", "sub.bed.matrix", init.sub.bed.matrix)
