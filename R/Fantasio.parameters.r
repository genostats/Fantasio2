#' Controlling Fantasio behaviour
#'
#' @param ... 
#'
#' @details If no arguments are given, the list of current parameters
#' values is returned. 
#' @details When some arguments are given, they are used to update 
#' the current values (cf example).
#' @details The parameters `m`, `epsilon`, `past`, `delta`, `max_iterations`, `max_submin`,
#' `max_linesearch`, `min_step`, `max_step`, `ftol` and `wolfe` are for L-BFGS-B optimization
#' algorithm control, and are described at 
#' \url{https://lbfgspp.statr.me/doc/classLBFGSpp_1_1LBFGSBParam.html}.
#' @details The parameter `max_retries` gives the number of times the algorithm can be
#' re-run after a runtimer error. The vectors of length 2 `lower` and `upper` give the 
#' lower and upper bounds for \eqn{a} and \eqn{f}.
#' @details The parameter `n_threads` controls the number of threads; and `debug` to a positive value
#' will turn on the displaying of various amounts of information.
#' 
#' @return A list, or `NULL`.
#' @export
#'
#' @examples
#' Fantasio.parameters()
#' Fantasio.parameters(delta = 1e-8)
#' Fantasio.parameters()

Fantasio.parameters <- function(...) {
  params <- getUserParam()
  update <- list(...)
  if(length(update) == 0) 
    return(params)

  for(a in names(update)) {
    if( !(a %in% names(params)) ) 
      stop(a, " is not a valid parameter")
    params[[a]] <- update[[a]]
  }
  if(params$n_threads > 1 & !checkOpenMP()) {
    params$n_treads <- 1
    warning("Fantasio was not compiled with OpenMP, multithreading is impossible")
  }
  do.call(setUserParam, params)
}
