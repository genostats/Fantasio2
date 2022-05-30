
#' @exportClass sub.bed.matrix
setClass("sub.bed.matrix", contains = "bed.matrix", slots = list(random.seed = "integer", submap = "integer"))

init.sub.bed.matrix <- function(.Object, ...) {
  .Object <- callNextMethod() 
  .Object@random.seed <- get(".Random.seed", envir = .GlobalEnv)
  arg <- list(...) 
  arg$restore.seed <- FALSE
  arg$sx <- .Object
  do.call(set.submap, arg)
}

setMethod("initialize", "sub.bed.matrix", init.sub.bed.matrix)
