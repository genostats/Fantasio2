require(Fantasio2)

if( !require("HGDP.CEPH") ) {
  install.packages("HGDP.CEPH", repos="https://genostats.github.io/R/")
  require("HGDP.CEPH")
}

if( !require("HumanGeneticMap") ) {
  install.packages("HumanGeneticMap", repos="https://genostats.github.io/R/")
  require("HumanGeneticMap")
}

if(!exists('x.be')) {
  filepath <- system.file("extdata", "hgdp_ceph.bed", package="HGDP.CEPH")
  x <- read.bed.matrix(filepath)
  x <- set.stats(x)
  x <- set.dist(x, HumanGeneticMap::genetic.map.b36) # vérifier le build ?

  x.be <- select.inds(x, population == "Bedouin")
  # x.be <- x.be[1:47,]
}

if(!exists("segment.list")) {
  segment.list <- segments.list.by.hotspots(x.be)
}


Fantasio.parameters(n_threads = 4)

# pour l'instant, que "by hotspots" avec un summary "by SNPs"
set.seed(1)
Fantasio2 <- function(bedmatrix, n, segment.options, min.quality = 95, list.id, probs = TRUE, phen.code = c("plink", "R"), 
                      q = 1e-4, epsilon = 1e-5) {

  if (missing(segment.options))
    segment.options <- list()
  # segment.list <- do.call(segments.list.by.hotspots, c(bedmatrix = bedmatrix, segment.options))
  n.segments <- sum( sapply(segment.list, length) ) 


  Seeds <- matrix( nrow = length( Fantasio2:::getRandomSeed() ), ncol = n )

  # ceci correspond à peu près à ce que faisait make Atlas suivi de festim
  A <- matrix( nrow = nrow(bedmatrix), ncol = n )
  F <- matrix( nrow = nrow(bedmatrix), ncol = n )

  # LIK.0 <- matrix( nrow = nrow(bedmatrix), ncol = n )
  # LIK.1 <- matrix( nrow = nrow(bedmatrix), ncol = n )
  P.LRT <- matrix( nrow = nrow(bedmatrix), ncol = n )

  for(i in 1:n) {
    Seeds[,i] <- Fantasio2:::getRandomSeed()
    submap <- Fantasio2:::rsubmap(segment.list)
    d.dist <- Fantasio2:::delta.dist.0(bedmatrix, submap)
    S <- Fantasio2:::festim(bedmatrix@bed, bedmatrix@p, submap, d.dist, epsilon)
    A[,i] <- S$a
    F[,i] <- S$f
    # LIK.0[,i] <- S$likelihood.0
    # LIK.1[,i] <- S$likelihood.1
    P.LRT[,i] <- pchisq( 2*(S$likelihood.1 - S$likelihood.0), df = 2, lower.tail = FALSE )
  }

  # A et F correspondent au contenu du slot estimation_summary (qui n'est pas un summary)

  # Ceci correspond au slot submap_summary
  summary <- submaps.summary(bedmatrix, A, F, P.LRT)

  # détermine les indices des individus sur lesquels on calcule HBD et FLOD 
  # (cas consanguins ou autre selon les valeurs list.id, probs, phen.code...)
  indexes <- Fantasio2:::which.inbreds(summary, min.quality, list.id, probs, phen.code = match.arg(phen.code))
  whichInds <- seq_len(nrow(bedmatrix)) %in% indexes$HBD

  # on passe à ce que faisait setHBDProbAndFLODBySnps
  old.seed <- Fantasio2:::getRandomSeed() # storing current seed
  h <- new.env()
  for(i in 1:n) {
    # on génère les mêmes cartes
    Fantasio2:::setRandomSeed(Seeds[,i])
    submap <- Fantasio2:::rsubmap(segment.list)
    d.dist <- Fantasio2:::delta.dist.0(bedmatrix, submap)
    # matrice des pHBD
    HBD <- Fantasio2:::probaHBD(bedmatrix@bed, bedmatrix@p, submap, d.dist, whichInds, a = A[,i], f = F[,i], epsilon)
    # matrice des FLOD (une colonne par individu)
    FLOD <- (HBD + q * (1 - HBD)) 
    # chaque colonne doit etre divisée par ( f + q * (1 - f) )
    f <- F[indexes$HBD, i]
    for(j in 1:ncol(FLOD)) 
      FLOD[,j] <- FLOD[,j] / ( f[j] + q * (1 - f[j]) )
    h <- Fantasio2:::updateHashProbas(h, submap, HBD, FLOD)
  }
  Fantasio2:::setRandomSeed(old.seed) # restauring seed
  # contenu de @HBD_recap et @FLOD_recap
  x <- Fantasio2:::hashProbasToMatrix(h)
  rownames(x$HBD) <- rownames(x$FLOD) <- Fantasio2:::uniqueIds(summary$famid[indexes$HBD], summary$id[indexes$HBD])
  colnames(x$HBD) <- colnames(x$FLOD) <- bedmatrix@snps$id[x$snp]
  x
}

debug(Fantasio2)
Fantasio2(x.be, 5, probs = TRUE, phen.code = "plink")

