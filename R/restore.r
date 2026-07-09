#' Restore HBD and FLOD matrices
#'
#' @rdname restore
#' @param atlas an atlas
#' 
#' @description if atlas was built with memory mapped objects for HBD and FLOD matrices, 
#' and their pointers is broken, this method will attempt to recreate a valid pointer, if the 
#' file still exists.
#'
#' @return an atlas
#'

#' @rdname restore
setMethod("restore", "atlas", 
  function(object) {
    if(is(object@HBD_recap, "mmatrix")) object@HBD_recap <- restore(object@HBD_recap)
    if(is(object@FLOD_recap, "mmatrix")) object@FLOD_recap <- restore(object@FLOD_recap)
    object
  }
)
