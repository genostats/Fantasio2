#' @export
recap.HBD.FLOD <- function(atlas, keep.inds, q, recap) {
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
  h <- new.env()
  for(i in 1:n) { # boucle sur les cartes
    # on re génère les cartes
    setRandomSeed(seeds[,i])
    submap <- rsubmap(segments.list)
    d.dist <- delta.dist(bedmatrix, submap)

    # les a et f pour les individus qui nous intéressent, pour la carte en cours
    a <- A[wi, i]
    f <- F[wi, i]

    # matrice des pHBD [une colonne par individu, une ligne par SNP]
    HBD <- probaHBD(bedmatrix@bed, bedmatrix@p, submap, d.dist, keep.inds, a = A[,i], f = F[,i], epsilon)
    HBD[!is.finite(HBD)] <- 0

    # matrice des FLOD (une colonne par individu)
    FLOD <- log10(HBD + q * (1 - HBD))
    # chaque colonne doit etre divisée par ( f + q * (1 - f) ) [ soustraction à l'échelle log10 ]
    # (on pourrait utiliser sweep mais niveau gestion mémoire ceci doit être plus efficace)
    for(j in 1:ncol(FLOD))
      FLOD[,j] <- FLOD[,j] - log10( f[j] + q * (1 - f[j]) )

    h <- updateHashProbas(h, submap, (a < 1), HBD, FLOD)
  }
  setRandomSeed(old.seed) # restoring seed

  # calcule les matrices moyennes des HBD / FLOD snp par snp
  # ces matrices ont une ligne par individu / une colonne par SNP
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
