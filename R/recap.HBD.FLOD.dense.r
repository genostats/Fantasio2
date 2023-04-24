#' @export
recap.HBD.FLOD.dense <- function(atlas, keep.inds, q, recap, median) {
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

  old.seed <- getRandomSeed() # storing current seed
  wi <- which(keep.inds)
  

  x <- list() #créer liste pour les 2 grandes matrices phbd et FLOD

  if(median==TRUE) {
    a <- summary$a_median
    f <- summary$f_median
    }
  
  #première boucle pour créer la grande sous-carte = union des snps tirés dans les n sous-cartes
  big.submap <- as.integer(vector())
  for(i in 1:n){
    setRandomSeed(seeds[,i])
    submap <- rsubmap(segments.list)
    big.submap <- union(big.submap,submap)
  }
  big.submap <- sort(big.submap) 
  
  big.HBD <- matrix(0, ncol = length(wi), nrow = length(big.submap)) #matrice des phbd avec 1 colonne par individu consanguin et 1 ligne par snp tiré
  big.FLOD <- matrix(0, ncol = length(wi), nrow = length(big.submap)) #matrice des flod avec 1 colonne par individu consanguin et 1 ligne par snp tiré
  
  
  for(i in 1:n) { # boucle sur les cartes
    # on re génère les cartes
    setRandomSeed(seeds[,i])
    submap <- rsubmap(segments.list)
    d.dist <- delta.dist(bedmatrix, big.submap)

    # les a et f pour les individus qui nous intéressent, pour la carte en cours
    #a <- A[wi, i]
    #f <- F[wi, i]
    
    # les a et f pour la carte en cours
    if(median==FALSE){
      a <- A[, i]
      f <- F[, i]
    }
    

    # matrice des pHBD [une colonne par individu, une ligne par SNP]
    
    #créer vecteur freq.submap de NA de longueur big.submap puis remplacer par les freq aux positions tirées
    freq.submap <- rep(NA, times = length(bedmatrix@p))
    freq.submap[submap] <- bedmatrix@p[submap]
    HBD <- probaHBD(bedmatrix@bed, freq.submap, big.submap, d.dist, keep.inds, a = a, f = f, epsilon) #renvoie les snps sur les lignes et les inds sur les colonnes
    HBD[!is.finite(HBD)] <- 0
   
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
  setRandomSeed(old.seed) # restoring seed

  
  # calcule les matrices moyennes des HBD / FLOD snp par snp
  # ces matrices doivent avoir une ligne par individu / une colonne par SNP
  x$HBD <- t(big.HBD)/n #on transpose pour la suite
  x$FLOD <- t(big.FLOD)/n #on transpose pour la suite
  rownames(x$HBD) <- rownames(x$FLOD) <- uniqueIds(summary$famid[wi], summary$id[wi])
  colnames(x$HBD) <- colnames(x$FLOD) <- bedmatrix@snps$id[big.submap] #snps vus dans grande sous carte

  # c'est fini.
  atlas@HBD_recap <- x$HBD     
  atlas@FLOD_recap <- x$FLOD   
  atlas@recap <- recap
  atlas@q <- q
  atlas
}
