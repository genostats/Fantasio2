#' Computation of HFLOD scores
#' 
#' This function is used to compute HFLOD scores on individuals in a sample
#' 
#' @param atlas a atlas object
#' @param w.id vector of indices to be considered in the FLOD_recap matrix
#' 
#' @return the atlas object with its slot HFLOD completed
#' @export
set.HFLOD <- function(atlas, w.id)
{
  HFLOD <- get.positions(atlas)

  HFLOD_value <- numeric(nrow(HFLOD))
  ALPHA_value <- numeric(nrow(HFLOD))
  
  h <- function(alpha, flods) {
     sum(log10(alpha * 10**flods + (1-alpha) ), na.rm = TRUE)
  }
  for (j in seq_len(nrow(HFLOD))) {
    # optimisation of h(alpha) ; 
    res <- optimize( h, c(0,1), flods = atlas@FLOD_recap[w.id, j], maximum = TRUE, tol = 0.001 )
    
    HFLOD_value[j] <- res$objective # HFLOD = h(alpha max)
    ALPHA_value[j] <- res$maximum   # alpha max 
  }
  HFLOD$HFLOD <- HFLOD_value
  HFLOD$ALPHA <- ALPHA_value

  atlas@HFLOD <- HFLOD
  atlas 
}
