#' Dense summary of pHBD and FLOD
#' 
#' pHBD and FLOD are combined over the submaps with a mean of pHBD and FLOD on all submaps
#'
#' @param atlas an atlas object
#' @param keep.inds the vector of index of consanguineous individuals
#' @param q assumed frequency of the mutation involved in the disease for each individual (default is 0.0001)
#' @param recap if you want the summary of probabilities by snps or by segments (only by SNPs for the moment)
#' @param median define the f and a parameters used to compute pHBD and FLOD
#'	   - if FALSE : f and a estimated on each submap
#'	   - if TRUE : median value of estimations on all submaps of f and a (default)
#' @param basename if missing, the HBD and FLOD matrices will be R native matrices, else they will be houba matrices
#'
#' @export
recap.HBD.FLOD.dense <- function(atlas, keep.inds, q, recap, median, basename) {
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
  
  # tester s'il y a des inds consanguins
  # si oui, remplir les matrices de FLOD et pHBD par NULL
  if(length(wi) == 0) {
    atlas@HBD_recap <- NULL     
    atlas@FLOD_recap <- NULL  
    atlas@recap <- recap
    atlas@q <- q
    atlas
  } else {  
    if(median) {
      a <- summary$a_median
      f <- summary$f_median
    }
 
    verbose <- Fantasio.parameters("verbose")
    if(verbose) cat("Merging submaps for dense HBD computation\n")

    # première boucle pour créer la grande sous-carte = union des snps tirés dans les n sous-cartes
    big.submap <- as.integer(vector())
    for(i in 1:n){
      setSeed(seeds[,i])
      submap <- rsubmap(segments.list)
      big.submap <- union(big.submap,submap)
    }
    big.submap <- sort(big.submap) 
  
    # matrice des phbd /des flod, avec 1 colonne par individu consanguin et 1 ligne par snp tiré
    if(missing(basename)) {
      big.HBD <- matrix(0, ncol = length(wi), nrow = length(big.submap)) 
      big.FLOD <- matrix(0, ncol = length(wi), nrow = length(big.submap)) 
      # la matrice pour chacune des sous carte
      HBD <- matrix(0, nrow = length(big.submap), ncol = sum(keep.inds))
    } else {
      stop("# houba bla")
    }
  
    for(i in 1:n) { # boucle sur les cartes
      if(verbose) cat("Computing HBD and FLOD using SNPs from submap", i, "\r")
      # on re génère les cartes
      setSeed(seeds[,i])
      submap <- rsubmap(segments.list)
      d.dist <- delta.dist(bedmatrix, big.submap)

      # les a et f pour la carte en cours
      if(!median){
        a <- A[, i]
        f <- F[, i]
      }
    
      # matrice des pHBD [une colonne par individu, une ligne par SNP]
    
      # créer vecteur freq.submap de NA de longueur ncol(bedmatrix) puis remplacer par les freq aux positions de la carte 
      freq.submap <- rep(NA, times = length(bedmatrix@p))
      freq.submap[submap] <- bedmatrix@p[submap]
      # va calculer les pHBD aux positions de big.submap avec les fréqs à NA sauf aux points de la carte courante
      # (freq à NA : proba d'émission mise à 1, équivalent à "tous les génotypes manquants à cette position")

      # on calcule les pHBD avec les snps sur les lignes et les inds sur les colonnes
      probaHBD_matrix(bedmatrix@bed, HBD, p = freq.submap, submap = big.submap, d.dist, keep.inds, a = a, f = f, epsilon) 
   
      # extraction du f pour les individus conservés
      ff <- f[wi]
      # matrice des FLOD (une colonne par individu)
      FLOD <- log10(HBD + q * (1 - HBD))
      # chaque colonne doit etre divisée par ( f + q * (1 - f) ) [ soustraction à l'échelle log10 ]
      # (on pourrait utiliser sweep mais niveau gestion mémoire ceci doit être plus efficace)
      for(j in 1:ncol(FLOD))
        FLOD[,j] <- FLOD[,j] - log10( ff[j] + q * (1 - ff[j]) )

    
      big.HBD <- big.HBD + HBD 
      big.FLOD <- big.FLOD + FLOD  
    }
    if(verbose) cat("\n")
  
    # calcule les matrices moyennes des HBD / FLOD snp par snp
    # ces matrices doivent avoir une ligne par individu / une colonne par SNP
    HBD <- t(big.HBD)/n #on transpose pour la suite
    FLOD <- t(big.FLOD)/n #on transpose pour la suite
    rownames(HBD) <- rownames(FLOD) <- uniqueIds(summary$famid[wi], summary$id[wi])
    colnames(HBD) <- colnames(FLOD) <- bedmatrix@snps$id[big.submap] #snps vus dans grande sous carte

    # c'est fini.
    atlas@HBD_recap <- HBD     
    atlas@FLOD_recap <- FLOD   
    atlas@recap <- recap
    atlas@q <- q
    atlas
  }
}
