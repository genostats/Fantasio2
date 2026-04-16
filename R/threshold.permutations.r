#' 
#' This function calculate the threshold of the association analysis with permutation method
#' 
#' @param   atlas an atlas object
#' @param 

#' @export threshold.permutations



threshold.permutations <- function(atlas, nb.perm = 1000, expl.var = c("FLOD", "pHBD"), phen, phen.code = c("R", "plink"), score, cores){

  # phenotype coding

  phen_code <- match.arg(phen.code)
  
  # explanatory variable
  
  expl_var <- match.arg(expl.var)

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

  nb.cases <- sum(pheno == 1)
  nb.controls <- sum(pheno == 0)

  z.max <- list()
  
  
  if(score) {
    #first get the variance for all permutations
    cases <- sample(which(pheno == 0 | pheno == 1), nb.cases, replace = FALSE) #only on non-NA phenotypes
    controls <- sample(which(pheno == 0 | pheno == 1)[-cases], nb.controls, replace = FALSE) #only on non-NA phenotypes
    pheno[cases] <- 1
    pheno[controls] <- 0
    
    reg <- HBD.glm(x = atlas, expl_var = expl_var, phen = pheno, phen.code = "R", score = TRUE, pval = FALSE)
    z.max[1] <- max(reg$z.value)
    sigma2 <- reg$variance
    
    #then compute for the rest of the permutations    
  } else {
    force(atlas)
  }
  
  fp <- Fantasio.parameters()
  fp$n_threads <- 1
  
  get.z.max.unit <- function(iteration) {
    cases <- sample(which(pheno == 0 | pheno == 1), nb.cases, replace = FALSE) #only on non-NA phenotypes
    controls <- sample(which(pheno == 0 | pheno == 1)[-cases], nb.controls, replace = FALSE) #only on non-NA phenotypes
    pheno[cases] <- 1
    pheno[controls] <- 0
    
    variance <- if (score) sigma2 else NULL

    reg <- HBD.glm(x = atlas, expl_var = expl_var, phen = pheno, phen.code = "R", score = score, variance = variance, pval = FALSE)

    return(max(reg$z.value))
    #z.max[iteration] <- max(reg$z.value)
  }


  calc.z.max <- function(fin, score){
    deb <- ifelse(score, 2,1)
    
  
    library(parallel)
    cat("cores =", cores, "\n")
    cl <- makeCluster(cores) # créer le cluster
    on.exit(stopCluster(cl), add=TRUE)
    clusterSetRNGStream(cl) # L ecuyer
    
    
    #clusterExport(cl, "fp")
    parLapply(cl, 1:cores, function(i) do.call(Fantasio.parameters, fp) )
    
    results.all.z.max <- parLapply(cl, deb:fin, get.z.max.unit) # boucle for inclue dans parLapply
    
  
    return(results.all.z.max)
    #z.max <- c(z.max, results.all.z.max)  
  }
  
  all.z.max <- calc.z.max(nb.perm, score)

  z.max <- c(z.max, all.z.max)
  z.max <- unlist(z.max)
  z.max.95 <- quantile(z.max, 0.95)

  res <- list(z.max = z.max, threshold = z.max.95, signif = (z.max > z.max.95))

  res

}

