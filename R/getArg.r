
getArg <- function(.name, .default, ...) {
  arg <- list(...)
  if(is.null(arg[[.name]]))
    .default
  else
    arg[[.name]]
}
