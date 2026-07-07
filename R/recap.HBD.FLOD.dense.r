#' Dense summary of pHBD and FLOD
#' 
#' pHBD and FLOD are combined over the submaps with a mean of pHBD and FLOD on all submaps
#'
#' @param atlas an atlas object
#' @param keep.inds the vector of index of consanguineous individuals
#' @param q assumed frequency of the mutation involved in the disease for each individual (default is 0.0001)
#' @param median define the f and a parameters used to compute pHBD and FLOD
#'	   - if FALSE : f and a estimated on each submap
#'	   - if TRUE : median value of estimations on all submaps of f and a (default)
#' @param basename if missing, the HBD and FLOD matrices will be R native matrices, else they will be houba matrices
#'
#' @export
recap.HBD.FLOD.dense <- function(atlas, keep.inds, q = 0.0001, median = TRUE, basename) {

  use.houba <- !missing(basename)
  if(use.houba) {
    check_file_exist(basename)
    hbd.file  <- paste0(path.expand(basename), ".hbd")
    flod.file <- paste0(path.expand(basename), ".flod")
  }

  # shortcuts for atlas slots
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
  # si non, remplir les matrices de FLOD et pHBD par NULL
  if(length(wi) == 0) {
    atlas@HBD_recap <- NULL     
    atlas@FLOD_recap <- NULL  
    atlas@q <- q
    return(atlas)
  }

  ## il y a des consanguins
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
  if(use.houba) {
    dt <- if(Fantasio.parameters("use_float")) "float" else "double"
    big.HBD <- mmatrix(dt, nrow = length(big.submap), ncol = length(wi))
    big.FLOD <- mmatrix(dt, nrow = length(big.submap), ncol = length(wi))
    HBD <- mmatrix(dt, nrow = length(big.submap), ncol = length(wi))
  } else {
    big.HBD <- matrix(0, ncol = length(wi), nrow = length(big.submap)) 
    big.FLOD <- matrix(0, ncol = length(wi), nrow = length(big.submap)) 
    # la matrice pour chacune des sous cartes
    HBD <- matrix(0, nrow = length(big.submap), ncol = sum(keep.inds))
  } 

  # extraction du f pour les individus conservés
  ff <- f[wi]

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
    if(use.houba) {
      probaHBD_mmatrix(bedmatrix@bed, HBD, p = freq.submap, submap = big.submap, d.dist, keep.inds, a = a, f = f, epsilon) 
      FLOD_update_mmatrix(HBD, big.FLOD, ff, q);
      houba::inplace.sum(big.HBD, HBD)
    } else {
      probaHBD_matrix(bedmatrix@bed, HBD, p = freq.submap, submap = big.submap, d.dist, keep.inds, a = a, f = f, epsilon) 
      FLOD_update_matrix(HBD, big.FLOD, ff, q);
      big.HBD <- big.HBD + HBD 
    }
    # NOTE au lieu de calculer FLOD à chaque sous carte puis de prendre la moyenne des FLOD, ou aurait pu
    # calculer FLOD sur la moyenne des P(HBD) ...
  }
  if(verbose) cat("\n")

  # calcule les matrices moyennes des HBD / FLOD snp par snp
  # ces matrices doivent avoir une ligne par individu / une colonne par SNP
  # -> on transpose
  if(use.houba) {
    houba::inplace.div(big.HBD, n)
    HBD <- transpose(big.HBD, hbd.file)
    houba::inplace.div(big.FLOD, n)
    FLOD <- transpose(big.FLOD, flod.file)
  } else {
    HBD <- t(big.HBD)/n   
    FLOD <- t(big.FLOD)/n 
  }

  rownames(HBD) <- rownames(FLOD) <- uniqueIds(summary$famid[wi], summary$id[wi])
  colnames(HBD) <- colnames(FLOD) <- bedmatrix@snps$id[big.submap] #snps vus dans grande sous carte

  # c'est fini.
  atlas@HBD_recap <- HBD     
  atlas@FLOD_recap <- FLOD   
  atlas@q <- q
  atlas
}
