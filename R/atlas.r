#' @export
atlas <- function(bedmatrix, segments.list, n, epsilon = 1e-3) {

  # si use_froh est vrai, on commence par utiliser gaston::fROH ave les paramètres
  # définis dans Fantasio.parameters()
  # ces valeurs seront utilisées comme point de départ
  pars <- Fantasio.parameters()
  if(pars$use_froh) {
    froh <- gaston::fROH(bedmatrix, "cM", "het", pars$minNbSNPs, pars$minROHlength, pars$minDistHet, pars$maxGapLength)
    # a estimé par a = 1/((1- f) * longueur moyenne HBD) 
    froh$aROH <- froh$nbSegments / (froh$length * (1 - froh$fROH))
    # NOTE : c'est NaN si pas de ROHs détectés. La fonction festim ne fera pas l'estimation pour a = NaN
  } else {
    # si on n'utilise pas fROH, on prendra les pts de départ f = 0.05 et a = 0.05
    froh = data.frame( fROH = rep(0.05, nrow(bedmatrix)), aROH = rep(0.05, nrow(bedmatrix)) )
  }

  # énumération des sous cartes
  # on copie la graine de R
  copyRseed();

  seeds <- matrix( nrow = length(getSeed()), ncol = n )

  # ceci correspond à peu près à ce que faisait make Atlas suivi de festim
  A <- matrix( nrow = nrow(bedmatrix), ncol = n )
  F <- matrix( nrow = nrow(bedmatrix), ncol = n )

  # LIK.0 <- matrix( nrow = nrow(bedmatrix), ncol = n )
  # LIK.1 <- matrix( nrow = nrow(bedmatrix), ncol = n )
  P.LRT <- matrix( nrow = nrow(bedmatrix), ncol = n )

  for(i in 1:n) {
    seeds[,i] <- getSeed()
    submap <- rsubmap(segments.list)
    # pour garder la graine de R synchrone avec la nôtre
    runif( length(submap) )
   
    d.dist <- delta.dist(bedmatrix, submap)
    S <- festim(bedmatrix@bed, bedmatrix@p, submap, d.dist, epsilon, froh$fROH, froh$aROH)
    A[,i] <- S$a
    F[,i] <- S$f
    # LIK.0[,i] <- S$likelihood.0
    # LIK.1[,i] <- S$likelihood.1
    P.LRT[,i] <- pchisq( 2*(S$likelihood.1 - S$likelihood.0), df = 2, lower.tail = FALSE )
  }

  # A et F correspondent au contenu du slot estimation_summary (qui n'est pas un summary)

  # Ceci correspond au slot submap_summary
  summary <- submaps.summary(bedmatrix, A, F, P.LRT)

  if(pars$use_froh) {
    summary <- cbind( froh[, c("fROH", "aROH") ], summary) 
  } 

  new("atlas", bedmatrix, seeds, epsilon, segments.list, list(a = A, f = F, p = P.LRT), summary)
}

