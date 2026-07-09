#' 
#' This function calculate the threshold of the association analysis with permutation method
#' 
#' @param   atlas an atlas object
#' @param 

#' @export threshold.permutations

threshold.permutations <- function(atlas, nb.perm = 1000, expl.var = c("FLOD", "pHBD"), phen, phen.code = c("R", "plink"), covar_df = NULL, covar = NULL, score, 
                                   centered = TRUE, cores){

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

  # if there are any NA's, build a new atlas / covar / pheno without those individuals 
  # this will avoid the repeated extraction of the relevant lines in glm.HBD.0
  # and will also simplify the permutations in get.z.max.unit

  NAs <- is.na(pheno)
  if(any(NAs)) {
    keep <- which(!NAs)

    atlas@bedmatrix@ped <- atlas@bedmatrix@ped[keep, ]
    atlas@submap_summary <- atlas@submap_summary[keep, ]

    # remove them also from HBD/FLOD matrix
    id.k <- paste0(atlas@submap_summary$famid, ":", atlas@submap_summary$id)
    HBD.k <- which(rownames(atlas@FLOD_recap) %in% id.k)

    # we just need to perform the extraction on the relevant matrix
    if(expl_var == "FLOD") {
      atlas@FLOD_recap <- atlas@FLOD_recap[HBD.k, ]
    } else {
      atlas@HBD_recap <- atlas@HBD_recap[HBD.k, ]
    }

    if(!is.null(covar)) covar <- covar[keep, ]
    pheno <- pheno[keep] 
  }

  if(centered) {
    if(expl_var == "FLOD") {
      rowMeansVar <- as.vector(rowMeans(atlas@FLOD_recap))
    } else {
      rowMeansVar <- as.vector(rowMeans(atlas@HBD_recap))
    }
  } else {
    rowMeansVar <- numeric(0) 
  }

  # to ease the subsequent steps we keep only the inbred individuals in the atlas / covar / pheno
  keep <- which(atlas@submap_summary$inbred)
  atlas@bedmatrix@ped <- atlas@bedmatrix@ped[keep, ]
  atlas@submap_summary <- atlas@submap_summary[keep, ]
  if(!is.null(covar)) covar <- covar[keep, ]
  pheno <- pheno[keep] 

  # run HBD.glm on real phenotype
  as <- HBD.glm(x = atlas, expl_var = expl_var, phen = pheno, covar_df = covar_df, covar = covar, phen.code = "R", score = score, pval = FALSE, 
                centered = centered, rowMeansVar = rowMeansVar)

  # keep the variance from 'as'
  sigma2 <- as$variance

  # prepare a set of parameters with n_threads = 1 (no multithreading in the cluster nodes !)
  fp <- Fantasio.parameters()
  fp$n_threads <- 1
 
  # the function which computes one (permutated) value of z.min and z.max
  get.z <- function(iteration) {
    pheno <- sample(pheno)
    variance <- if (score) sigma2 else NULL
    reg <- HBD.glm(x = atlas, expl_var = expl_var, phen = pheno, covar_df = covar_df, covar = covar, phen.code = "R", score = score, 
                   variance = variance, pval = FALSE, centered = centered, rowMeansVar = rowMeansVar)
    return(c(zmax = max(reg$z.value), zmin = min(reg$z.value)))
  }

  calc.z <- function(fin, score){
    deb <- 1
  
    library(parallel)
    cat("cores =", cores, "\n")
    cl <- makeCluster(cores) # créer le cluster
    on.exit(stopCluster(cl), add=TRUE)
    clusterSetRNGStream(cl) # L ecuyer
    
    # no multithreading in the nodes
    parLapply(cl, 1:cores, function(i) do.call(Fantasio.parameters, fp) )
    
    results.all.z <- parLapply(cl, deb:fin, get.z) # boucle for inclue dans parLapply
  
    return(results.all.z)
  }
  
  z.max.min <- calc.z(nb.perm, score)

  z.max <- sapply(z.max.min, function(x) x['zmax'])
  z.min <- sapply(z.max.min, function(x) x['zmin'])
  z.max.95 <- quantile(z.max, 0.95)
  z.min.95 <- quantile(z.min, 0.05)
  p.left <- -log10(pnorm(z.min.95))
  p.bil <- -log10(pchisq(max(abs(z.min.95), abs(z.max.95))**2, df = 1, lower.tail = FALSE))
  p.right <- -log10(1-pnorm(z.max.95))
  true.max <- max(as$z.value)

  res <- list(z.max = z.max, threshold.z.max = z.max.95,
              z.min = z.min, threshold.z.min = z.min.95,
              p.left = p.left, p.bilateral = p.bil, p.right = p.right,
              signif = (true.max > z.max.95))

  res
}

