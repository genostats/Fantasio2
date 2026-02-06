#' Sparse summary of pHBD and FLOD
#' 
#' pHBD and FLOD are combined over the submaps SNP by SNP with a mean on the number of submaps that includes this SNP
#'
#' @param atlas an atlas object
#' @param keep.inds the vector of index of consanguineous individuals
#' @param q assumed frequency of the mutation involved in the disease for each individual (default is 0.0001)
#' @param recap if you want the summary of probabilities by snps or by segments (only by SNPs for the moment)
#' @param median define the f and a parameters used to compute pHBD and FLOD
#'	   - if FALSE : f and a estimated on each submap
#'	   - if TRUE : median value of estimations on all submaps of f and a (default)
#'
#' @export
recap.HBD.FLOD.sparse <- function(atlas, keep.inds, q, recap, median) {
  if(recap != "SNP") stop("Not yet implemented")

  # shotcuts for atlas slots
  n <- ncol(atlas@seeds)
  bedmatrix <- atlas@bedmatrix
  seeds <- atlas@seeds
  segments.list <- atlas@segments_list
  A <- atlas@estimations$a
  F <- atlas@estimations$f
  summary <- atlas@submap_summary
  epsilon <- atlas@epsilon

  wi <- which(keep.inds)
  h <- new.env()
  
  #tester s'il y a des inds consanguins
  #si oui, remplir les matrices de FLOD et pHBD par NULL
  if(length(wi)==0){
    atlas@HBD_recap <- NULL     
    atlas@FLOD_recap <- NULL  
    atlas@recap <- recap
    atlas@q <- q
    atlas
  }else{  

    if(median) {
      a <- summary$a_median
      f <- summary$f_median
    }

    verbose <- Fantasio.parameters("verbose")

    for(i in 1:n) { # boucle sur les cartes
      if(verbose) cat("Computing HBD and FLOD on submap", i, "\r")
      # on re génère les cartes
      setSeed(seeds[,i])
      submap <- rsubmap(segments.list)
      d.dist <- delta.dist(bedmatrix, submap)

      if(!median) {
        # les a et f pour la carte en cours
        a <- A[, i]
        f <- F[, i]
      }

      # matrice des pHBD [une colonne par individu, une ligne par SNP]
      HBD <- probaHBD(bedmatrix@bed, bedmatrix@p, submap, d.dist, keep.inds, a = a, f = f, epsilon)
      HBD[!is.finite(HBD)] <- 0

      # les f et a pour les individus de keep inds
      ff <- f[wi]
      aa <- a[wi]
      # matrice des FLOD (une colonne par individu)
      FLOD <- log10(HBD + q * (1 - HBD))
      # chaque colonne doit etre divisée par ( f + q * (1 - f) ) [ soustraction à l'échelle log10 ]
      # (on pourrait utiliser sweep mais niveau gestion mémoire ceci doit être plus efficace)
      for(j in 1:ncol(FLOD))
        FLOD[,j] <- FLOD[,j] - log10( ff[j] + q * (1 - ff[j]) )

      h <- updateHashProbas(h, submap, (aa < 1), HBD, FLOD)
    }

    # calcule les matrices moyennes des HBD / FLOD snp par snp
    # ces matrices ont une ligne par individu / une colonne par SNP
    if(verbose) cat("Merging HBD and FLOD values\n")
    x <- hashProbasToMatrix(h)
    rownames(x$HBD) <- rownames(x$FLOD) <- uniqueIds(summary$famid[wi], summary$id[wi])
    colnames(x$HBD) <- colnames(x$FLOD) <- bedmatrix@snps$id[x$snp]

    # c'est fini.
    atlas@HBD_recap <- x$HBD
    atlas@FLOD_recap <- x$FLOD
    atlas@recap <- recap
    atlas@q <- q
    atlas
  }
}
