#' @useDynLib Fantasio2, .registration=TRUE
#' @importFrom Rcpp evalCpp
NULL

# setting two parameters of the optimization algorithm.
# The default lower bound for a is 0.01. 
# In many cases when both a and f are near to 0 the likelihood is very flat,
# and the optimization algorithm might fail to converge.
.onLoad <- function(libname, pkgname) {
  Fantasio.parameters(lower = c(0.01, 0), upper = c(Inf, 0.999))
}
