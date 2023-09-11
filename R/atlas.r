#' @export
atlas <- function(bedmatrix, segments.list, n, epsilon = 1e-3) {

  # on copie la graine de R
  copyRseed();

  seeds <- matrix( nrow = 625, ncol = n )

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
    S <- festim(bedmatrix@bed, bedmatrix@p, submap, d.dist, epsilon)
    A[,i] <- S$a
    F[,i] <- S$f
    # LIK.0[,i] <- S$likelihood.0
    # LIK.1[,i] <- S$likelihood.1
    P.LRT[,i] <- pchisq( 2*(S$likelihood.1 - S$likelihood.0), df = 2, lower.tail = FALSE )
  }

  # A et F correspondent au contenu du slot estimation_summary (qui n'est pas un summary)

  # Ceci correspond au slot submap_summary
  summary <- submaps.summary(bedmatrix, A, F, P.LRT)
  new("atlas", bedmatrix, seeds, epsilon, segments.list, list(a = A, f = F, p = P.LRT), summary)
}

