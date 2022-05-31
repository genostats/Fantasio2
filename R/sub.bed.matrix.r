
#' @exportClass sub.bed.matrix
setClass("sub.bed.matrix", contains = "bed.matrix", slots = list(random.seed = "integer", submap = "integer"))

#' @export
sub.bed.matrix <- function(bed.matrix, segment.list) {
  sx <- new("sub.bed.matrix", bed.matrix, random.seed = get(".Random.seed", envir = .GlobalEnv))
  set.submap(sx, segment.list, restore.seed = FALSE)
}
