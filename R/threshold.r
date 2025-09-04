

#acfd function
acfd <- function(d, z, Chr, Dist) {
    .Call(`_Fantasio2_acfd`, d, z, Chr, Dist)
}

#max Gaussian function
maxGaussianGenome <- function(nSimus, lambda, sigma, dist, len, seed, nThreads) {
    .Call(`_Fantasio2_maxGaussianGenome`, nSimus, lambda, sigma, dist, len, seed, nThreads)
}


# as = un data frame avec chr, dist, z.value
threshold <- function(as, n.sims = 1e4, n.cores = 10) {
  # cette version n'utilise plus du tout la vraisemblance
  # le lamba est calculé à partir de la fonction d'auto-corrélation
  qq <- quantile(as$z.value, c(0.05, 0.95))
  mu <- mean(qq)
  s2 <- as.numeric( (diff(qq)/3.2897)^2 )
  d <- seq(5, 100, by = 5)
  cc <- sapply(d, acfd, z = ((as$z.value - mu)/sqrt(s2)), Chr = as$chr, Dist = as$dist)
  cc[ cc < 0 ] <- 0.01
  lambda <- -as.numeric(lm(log(cc) ~ d)$coeff[2])
  # peut arriver si beaucoup de bruit
  if(lambda < 0) lambda <- -as.numeric(lm(log(cc[1:10]) ~ d[1:10])$coeff[2])
  if(lambda < 0) lambda <- -as.numeric(lm(log(cc[1:5]) ~ d[1:5])$coeff[2])
  if(lambda < 0) { 
    warning("Failed to estimate lambda, taking lambda = 1")
    lambda <- 1 # loci très peu dépendant, très conservateur a priori
  }

  # on fait des simus pour chaque chromosome en gardant le max
  # longueur des chromosomes
  chrs <- unique(as$chr)
  LE <- NULL
  for(c in chrs) {
    le <- diff( range( as$dist[ as$chr == c ] ) )
    LE <- c(LE, le)
  }
  
  # fonction de simus en cpp
  # dist = 0.1 produit une valeur plus grande (plus juste si le processus était vraiment Gaussien)
  # la valeur produite par dist = 1 semble largement assez grande
  M <- maxGaussianGenome(1e5, lambda, sqrt(s2), dist = 1, len = LE, seed = 2^32 * runif(1), nThreads = 4)
  
  # et voilà le seuil pour le test à droite...
  threshold <- mu + quantile(M, 0.95)
  z.max <- max(as$z.value)
  c(z.max = z.max, mu = mu, s2 = s2, lambda = lambda, threshold = threshold, signif = (z.max > threshold))
}