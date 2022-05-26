require(Fantasio2)

if( !require("HGDP.CEPH") ) {
  install.packages("HGDP.CEPH", repos="https://genostats.github.io/R/")
  require("HGDP.CEPH")
}

if( !require("HumanGeneticMap") ) {
  install.packages("HumanGeneticMap", repos="https://genostats.github.io/R/")
  require("HumanGeneticMap")
}

filepath <- system.file("extdata", "hgdp_ceph.bed", package="HGDP.CEPH")
x <- read.bed.matrix(filepath)
x <- set.stats(x)
x <- set.dist(x, HumanGeneticMap::genetic.map.b36) # vérifier le build ?

x.be <- select.inds(x, population == "Bedouin")

set.seed(1); sx <- new("sub.bed.matrix", x.be, submap = sort(sample.int(ncol(x.be), 100)))

A <- Fantasio2:::m4_logEmiss( sx@bed, sx@p, sx@submap, 1e-5 )

cat("test log emission\n")
for(k in 0:47) {
  B <- Fantasio2:::testLogEmiss( sx@bed, sx@p, sx@submap, 1e-5, k )
  if(!all( A[ k*2 + 1:2, ] == B)) stop(k)
}

cat("test log likelihood\n")
a <- 2
for(k in 0:47) {
  R1 <- Fantasio2:::testLikelihood( sx@bed, sx@p, sx@submap, delta.dist(sx), 1e-5, k, a, 0.1)
  R2 <- Fantasio2:::logLikelihood_gradient( A[k*2 + 1:2, ],  delta.dist(sx), a, 0.1)
  if(!all(R1 == R2)) stop(k)

  R1 <- Fantasio2:::testLikelihood( sx@bed, sx@p, sx@submap, delta.dist(sx), 1e-5, k, a, 1)
  R2 <- Fantasio2:::logLikelihood_gradient( A[k*2 + 1:2, ],  delta.dist(sx), a, 1)
  if(!all(R1 == R2)) stop(k)

  R1 <- Fantasio2:::testLikelihood( sx@bed, sx@p, sx@submap, delta.dist(sx), 1e-5, k, a, 0)
  R2 <- Fantasio2:::logLikelihood_gradient( A[k*2 + 1:2, ],  delta.dist(sx), a, 0)
  if(!all(R1 == R2)) stop(k)
}


## Il y a un souci de discontinuité du gradient en a = 0
## -> revoir les calculs papier...

k <- 17
Fantasio2:::testLikelihood( sx@bed, sx@p, sx@submap, delta.dist(sx), 1e-5, k, a = 1e-22, 0.052)
# [1] -66.86574837 -51.74281869   0.02803073
Fantasio2:::testLikelihood( sx@bed, sx@p, sx@submap, delta.dist(sx), 1e-5, k, a = 0, 0.052)
# [1]  -66.86574837 1153.61990432    0.02803073
k <- 2
Fantasio2:::testOptimLikelihood( sx@bed, sx@p, sx@submap, delta.dist(sx), 1e-5, k)




f <- function(x) -Fantasio2:::logLikelihood_gradient( A[k*2 + 1:2, ],  delta.dist(sx), x[1], x[2] )[1]
grad <- function(x) -Fantasio2:::logLikelihood_gradient( A[k*2 + 1:2, ],  delta.dist(sx), x[1], x[2] )[-1]

optim( c(0.1, 0.1), fn = f, gr = grad, method = "L-BFGS-B", lower = c(0,0), upper = c(Inf, 1) )
