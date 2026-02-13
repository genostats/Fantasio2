#' 
#' This function calculate the threshold of the association analysis with permutation method
#' 
#' @param   atlas an atlas object
#' @param 

#' @export threshold.permutations

threshold.permutations <- function(atlas, nb.perm = 1000, expl.var = c("FLOD", "pHBD"), phen, phen.code = c("R", "plink"), score){

  # phenotype coding

  phen_code <- match.arg(phen.code)

  # recover phenotype

  if(missing(phen)){
    pheno <- atlas@bedmatrix@ped$pheno
  }
  else{
    pheno <- phen
  }

  if (phen.code == 'plink') {
    pheno <- ifelse(pheno == 1, 0, ifelse(pheno == 2, 1, NA))# Translate phenotype
  }

  # now the coding is 0:control / 1:case / NA

  nb.cases <- sum(pheno == 1)
  nb.controls <- sum(pheno == 0)

  z.max <- vector()
  
  if(score) {
    #first get the variance for all permutations
    cases <- sample(which(pheno == 0 | pheno == 1), nb.cases, replace = FALSE) #only on non-NA phenotypes
    controls <- sample(which(pheno == 0 | pheno == 1)[-cases], nb.controls, replace = FALSE) #only on non-NA phenotypes
    pheno[cases] <- 1
    pheno[controls] <- 0
    
    reg <- HBD.glm(x = atlas, expl_var = expl_var, phen = pheno, phen.code = "R", score = TRUE)
    z.max[1] <- max(reg$z.value)
    sigma2 <- reg$variance
    
    #then compute for the rest of the permutations
    for(j in 2:nb.perm){
      cases <- sample(which(pheno == 0 | pheno == 1), nb.cases, replace = FALSE) #only on non-NA phenotypes
      controls <- sample(which(pheno == 0 | pheno == 1)[-cases], nb.controls, replace = FALSE) #only on non-NA phenotypes
      pheno[cases] <- 1
      pheno[controls] <- 0

      reg <- HBD.glm(x = atlas, expl_var = expl_var, phen = pheno, phen.code = "R", score = TRUE, variance = sigma2)
      z.max[j] <- max(reg$z.value)
    }
    
  } else {
    for(j in 1:nb.perm){
      cases <- sample(which(pheno == 0 | pheno == 1), nb.cases, replace = FALSE) #only on non-NA phenotypes
      controls <- sample(which(pheno == 0 | pheno == 1)[-cases], nb.controls, replace = FALSE) #only on non-NA phenotypes
      pheno[cases] <- 1
      pheno[controls] <- 0

      reg <- HBD.glm(x = atlas, expl_var = expl_var, phen = pheno, phen.code = "R", score = FALSE)
      z.max[j] <- max(reg$z.value)
    }
  }

  z.max
  z.max.95 <- quantile(z.max, 0.95)

  res <- list(z.max = z.max, z.max.95 = z.max.95)

  res

}

