#' 
#' This function calculate the threshold of the association analysis with permutation method for the HBD-GWAS
#' 
#' @param   atlas an atlas object
#' @param 

#' @export threshold.permutations.HFLOD



threshold.permutations.HFLOD <- function(atlas, nb.perm = 1000, phen, phen.code = c("R", "plink"), cores){

  # phenotype coding

  phen_code <- match.arg(phen.code)
  

  # recover phenotype

  if(missing(phen)) {
    pheno <- atlas@bedmatrix@ped$pheno
  } else {
    pheno <- phen
  }

  if (phen.code == 'plink') {
    pheno <- ifelse(pheno == 1, 0, ifelse(pheno == 2, 1, NA))# Translate phenotype
  }

  NAs <- is.na(pheno)
  if(any(NAs)) {
    keep <- which(!NAs)

    atlas@bedmatrix@ped <- atlas@bedmatrix@ped[keep, ]
    atlas@submap_summary <- atlas@submap_summary[keep, ]

    # remove them also from FLOD matrix
    id.k <- paste0(atlas@submap_summary$famid, ":", atlas@submap_summary$id)
    HBD.k <- which(rownames(atlas@FLOD_recap) %in% id.k)

    # we just need to perform the extraction on the FLOD matrix

    atlas@FLOD_recap <- atlas@FLOD_recap[HBD.k, ]

    pheno <- pheno[keep] 
  }

  # to ease the subsequent steps we keep only the inbred individuals in the atlas / covar / pheno
  keep <- which(atlas@submap_summary$inbred)
  atlas@bedmatrix@ped <- atlas@bedmatrix@ped[keep, ]
  atlas@submap_summary <- atlas@submap_summary[keep, ]
  pheno <- pheno[keep] 

  
  #run HBD.gwas on real phenotype
  hflod <- HBD.gwas(atlas = atlas, phen = pheno, phen.code = "R")
  
  # prepare a set of parameters with n_threads = 1 (no multithreading in the cluster nodes !)
  fp <- Fantasio.parameters()
  fp$n_threads <- 1
 
  # the function which computes one (permutated) value of z.min and z.max
  get.hflod <- function(iteration) {
    pheno <- sample(pheno)
    hg <- HBD.gwas(atlas = atlas, phen = pheno, phen.code = "R")
    return(hflodmax = max(hg$HFLOD))
  }

  calc.hflod <- function(fin){
    deb <- 1
  
    library(parallel)
    cat("cores =", cores, "\n")
    cl <- makeCluster(cores) # créer le cluster
    on.exit(stopCluster(cl), add=TRUE)
    clusterSetRNGStream(cl) # L ecuyer
    
    # no multithreading in the nodes
    parLapply(cl, 1:cores, function(i) do.call(Fantasio.parameters, fp) )
    
    results.all.hflod <- parLapply(cl, deb:fin, get.hflod) # boucle for inclue dans parLapply
  
    return(results.all.hflod)
  }
  
  all.hflod <- calc.hflod(nb.perm)
  
  HFLOD.max <- unlist(all.hflod)
  HFLOD.max.95 <- quantile(HFLOD.max, 0.95)
  
  true.max <- max(hflod$HFLOD)

  res <- list(HFLOD.max = HFLOD.max, threshold = HFLOD.max.95, signif = (true.max > HFLOD.max.95))

  res

}
