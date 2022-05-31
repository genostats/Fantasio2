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
}

# -- définition d'une sub bed matrix
s <- segments.list.by.hotspots(x.be)
set.seed(1); sx <- sub.bed.matrix(x.be, s)

# ------------------


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




R <- NULL
for(k in 0:47) {
  f <- function(x) -Fantasio2:::logLikelihood_gradient( A[k*2 + 1:2, ],  delta.dist(sx), x[1], x[2] )[1]
  grad <- function(x) -Fantasio2:::logLikelihood_gradient( A[k*2 + 1:2, ],  delta.dist(sx), x[1], x[2] )[-1]

  o <- optim( c(0.1, 0.1), fn = f, gr = grad, method = "L-BFGS-B", lower = c(1e-10,0), upper = c(Inf, 1) )
  R <- rbind(R, c(o$par, o$value))
}
R <- as.data.frame(R)
names(R) <- c("a", "f", "fx")

S <- NULL
for(k in 0:47) {
  o <- Fantasio2:::testOptimLikelihood( sx@bed, sx@p, sx@submap, delta.dist(sx), 1e-5, k)
  S <- rbind(S, o)
}
S <- as.data.frame(S)

par(mfrow=c(1,2))
plot(R$a, S$a)
plot(R$f, S$f)

