#' Controlling optimization algorithm
#'
#' @param ... 
#'
#' @details If no arguments are given, the list of current parameter 
#' values is returned. A description of the parameters can be seen at
#' \url{https://lbfgspp.statr.me/doc/classLBFGSpp_1_1LBFGSBParam.html}.
#' @details Additionally, the vectors of length 2 `lower` and `upper` 
#' allow to fix the bounds for `a` and `f`.
#' @details When some arguments are given, they are used to update 
#' the current values (cf example).
#' 
#' @return A list or `NULL`.
#' @export
#'
#' @examples
#' optimization.control()
#' optimization.control(epsilon = 1e-6)
#' optimization.control()

optimization.control <- function(...) {
  params <- getUserParam()
  update <- list(...)
  if(length(update) == 0) 
    return(params)

  for(a in names(update)) {
    if( !(a %in% names(params)) ) 
      stop(a, " is not a optimization parameter")
    params[[a]] <- update[[a]]
  }
  do.call(setUserParam, params)
}
