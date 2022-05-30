
#' @export
unset.submap <- function(sub.bed.matrix) {
  sub.bed.matrix@submap <- integer(0)
  sub.bed.matrix
}
