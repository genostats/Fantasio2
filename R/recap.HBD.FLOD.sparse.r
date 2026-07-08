#' Sparse summary of pHBD and FLOD
#' 
#' pHBD and FLOD are combined over the submaps SNP by SNP with a mean on the number of submaps that includes this SNP
#'
#' @param atlas an atlas object
#' @param keep.inds the vector of index of consanguineous individuals
#' @param q assumed frequency of the mutation involved in the disease for each individual (default is 0.0001)
#' @param median define the f and a parameters used to compute pHBD and FLOD
#'	   - if FALSE : f and a estimated on each submap
#'	   - if TRUE : median value of estimations on all submaps of f and a (default)
#' @export
recap.HBD.FLOD.sparse <- function(atlas, keep.inds, q, median) {

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
  
  # tester s'il y a des inds consanguins
   #si non, remplir les matrices de FLOD et pHBD par NULL
  if(length(wi) == 0) {
    atlas@HBD_recap <- NULL     
    atlas@FLOD_recap <- NULL  
    atlas@q <- q
    return(atlas)
  } 

  # il y en a...
  if(median) {
    a <- summary$a_median
    f <- summary$f_median
  }
  # les f et a pour les individus de keep inds
  ff <- f[wi]
  aa <- a[wi]

  verbose <- Fantasio.parameters("verbose")

  # nb markers par submap
  nbSNPs <- sum( segments.list.summary(atlas@segments_list)$number_of_segments )
  HBD <- matrix(0, nrow = nbSNPs, ncol = sum(keep.inds))
  FLOD <- matrix(0, nrow = nbSNPs, ncol = sum(keep.inds))

  for(i in 1:n) { # boucle sur les cartes
    if(verbose) cat("Computing HBD and FLOD on submap", i, "\r")
    # on re génère les cartes
    setSeed(seeds[,i])
    submap <- rsubmap(segments.list)  # toujours longueur nbSNPs
    d.dist <- delta.dist(bedmatrix, submap)

    if(!median) {
      # les a et f pour la carte en cours
      a <- A[, i]
      f <- F[, i]
      # aa et ff aussi doivent être mis à jour
      ff <- f[wi]
      aa <- a[wi]
    }

    # on calcule les pHBD avec les snps sur les lignes et les inds sur les colonnes
    probaHBD_matrix(bedmatrix@bed, HBD, bedmatrix@p, submap, d.dist, keep.inds, a = a, f = f, epsilon) 

    # les f et a pour les individus de keep inds
    # matrice des FLOD (une colonne par individu)
    FLOD[] <- 0 # il faut mettre à zéro car la fonction C++ ajoute à la matrice FLOD 
    FLOD_update_matrix(HBD, FLOD, ff, q);

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
  atlas@q <- q
  atlas
}
