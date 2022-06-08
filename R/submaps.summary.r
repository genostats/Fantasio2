#' Submap summary
#' 
#' This function creates a summary on the submaps created. 
#' 
#' @param bedmatrix the original bedmatrix
#' @param a.threshold  the maximum value for a (default is 1)
#' 
#' @details This function gives for each genotyped individual summary statistics about the calculations.
#' @details This function returns a dataframe with 10 columns :
#' @details - famid : family identifier
#' @details - id : individual identifier
#' @details - pheno : phenotype (1 non affected, 2 affected, 0 unknown)
#' @details - submaps : number of valid submaps (i.e. submaps with a < a.threshold)
#' @details - quality: percentage of valid submaps 
#' @details - f_median: median f on valid submaps (recommended to estimate f)
#' @details - a_median: median a on valid submaps (recommended to estimate a)
#' @details - pLRT_median: median p-value of LRT tests on valid submaps
#' @details - inbred: a flag indicating if the individual is inbred (pLRT_median <0.05) or not
#' @details - pLRT_<0.05: number of valid submaps with a LRT having a p-value below 0.05

#' @return this function returns a dataframe.
#' 
#' @seealso setSummary
#' 
#' @export

submaps.summary <- function (bedmatrix, a, f, a.threshold = 1) {
  
  w.a <- (a > a.threshold)
  f[w.a] <- NA
  a[w.a] <- NA
  p.lrt[w.a] <- NA
  l <- (p.lrt < 0.05)
  pLRT_median <- apply(p, 1, median, na.rm=TRUE)
  
  nValidSubmap <- numeric(nrow(l))
  for (i in seq_len(nrow(l)))
  {
    nValidSubmap[i] <- sum(l[i,], na.rm = TRUE)
  }
  
  submaps_used <- rowSums(sapply(as.data.frame(a), function(x) x <= a.threshold), na.rm = TRUE )
  quality <- (submaps_used*100)/dim(f)[2]
  
  df <- data.frame(famid  = bedmatrix@ped$famid,
                   id     = bedmatrix@ped$id,
                   pheno  = bedmatrix@ped$pheno,
                   submaps       = submaps_used,
                   quality       = quality,
                   f_median      = apply(f, 1, median, na.rm=TRUE),
                   a_median      = apply(a, 1, median, na.rm=TRUE),
                   pLRT_median   = pLRT_median,
                   inbred        = pLRT_median < 0.05,
                   pLRT_inf_0.05 = nValidSubmap
  )
  for(i in seq_len(nrow(df))) {
    if(!is.finite(df$f_median[i]))
      df[i,5:10] <- NA
  }
  
  return(df)
}

