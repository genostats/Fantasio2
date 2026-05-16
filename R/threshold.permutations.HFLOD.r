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

  # now the coding is 0:control / 1:case / NA

  nb.cases <- sum(pheno == 1, na.rm=TRUE)
  nb.controls <- sum(pheno == 0, na.rm=TRUE)

  HFLOD.max <- list()
  
  #run HBD.gwas on real phenotype
  hflod <- HBD.gwas(atlas = atlas, phen = pheno, phen.code = "R")
  
#  if(score) {
#    #first get the variance for all permutations
#    cases <- sample(which(pheno == 0 | pheno == 1), nb.cases, replace = FALSE) #only on non-NA phenotypes
#    controls <- sample(which(pheno == 0 | pheno == 1)[-cases], nb.controls, replace = FALSE) #only on non-NA phenotypes
#    pheno[cases] <- 1
#    pheno[controls] <- 0
    
#    reg <- HBD.glm(x = atlas, expl_var = expl_var, phen = pheno, phen.code = "R", score = TRUE, pval = FALSE)
#    z.max.min[[1]] <- c(zmax = max(reg$z.value), zmin = min(reg$z.value))
#    sigma2 <- reg$variance
    
    #then compute for the rest of the permutations    
#  } else {
#    force(atlas)
#  }
  
  fp <- Fantasio.parameters()
  fp$n_threads <- 1
  
  get.HFLOD.max.unit <- function(iteration) {
    cases <- sample(which(pheno == 0 | pheno == 1), nb.cases, replace = FALSE) #only on non-NA phenotypes
    controls <- sample(which(pheno == 0 | pheno == 1)[-cases], nb.controls, replace = FALSE) #only on non-NA phenotypes
    pheno[cases] <- 1
    pheno[controls] <- 0

    hg <- HBD.gwas(atlas = atlas, phen = pheno, phen.code = "R")

    return(hflodmax = max(hg$HFLOD))
    #z.max[iteration] <- max(reg$z.value)
  }


  calc.HFLOD.max <- function(fin){
    deb <- 1
    
  
    library(parallel)
    cat("cores =", cores, "\n")
    cl <- makeCluster(cores) # créer le cluster
    on.exit(stopCluster(cl), add=TRUE)
    clusterSetRNGStream(cl) # L ecuyer
    
    
    #clusterExport(cl, "fp")
    parLapply(cl, 1:cores, function(i) do.call(Fantasio.parameters, fp) )
    
    results.all.HFLOD.max <- parLapply(cl, deb:fin, get.HFLOD.max.unit) # boucle for inclue dans parLapply
    
  
    return(results.all.HFLOD.max)
    #z.max <- c(z.max, results.all.z.max)  
  }
  
  all.HFLOD.max <- calc.HFLOD.max(nb.perm)

  HFLOD.max <- unlist(all.HFLOD.max)
  HFLOD.max.95 <- quantile(HFLOD.max, 0.95)
  
  true.max <- max(hflod$HFLOD)

  res <- list(HFLOD.max = HFLOD.max, threshold = HFLOD.max.95, signif = (true.max > HFLOD.max.95))

  res

}
