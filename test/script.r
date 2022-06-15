require(gaston)

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
  x.be <- x.be[1:47,]
}

require(Fantasio2)

if(!exists("s")) {
  s <- segments.list.by.hotspots(x.be)
}

set.seed(4); sx <- sub.bed.matrix(x.be, s)

test.log.emission <- function() {
  cat("test log emission\n")
  A <- Fantasio2:::m4_logEmiss( sx@bed, sx@p, sx@submap, 1e-5 )
  for(k in 0:(nrow(x.be)-1)) {
    B <- Fantasio2:::testLogEmiss( sx@bed, sx@p, sx@submap, 1e-5, k )
    if(!all( A[ k*2 + 1:2, ] == B)) stop(k)
  }
}


test.log.likelihood <- function(a = 2) {
  cat("test log likelihood\n")
  A <- Fantasio2:::m4_logEmiss( sx@bed, sx@p, sx@submap, 1e-5 )
  for(k in 0:(nrow(x.be)-1)) {
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
}  
 
test.forward.backward <- function(k = 2, a = 0.025, f = 0.05) { 
  cat("test forward backward\n")
  
  A <- Fantasio2:::m4_logEmiss( sx@bed, sx@p, sx@submap, 1e-5 )
  R1 <- Fantasio2:::testForwardBackward( sx@bed, sx@p, sx@submap, delta.dist(sx), 1e-5, k, a, f)
  R2 <- Fantasio2:::forward_backward( A[k*2 + 1:2, ],  delta.dist(sx), a, f ) 
  list(new = R1, old = R2[2,])
}


## Il y a un souci de discontinuité du gradient en a = 0
## -> revoir les calculs papier...

k <- 17
Fantasio2:::testLikelihood( sx@bed, sx@p, sx@submap, delta.dist(sx), 1e-5, k, a = 1e-22, 0.052)
# [1] -66.86574837 -51.74281869   0.02803073
Fantasio2:::testLikelihood( sx@bed, sx@p, sx@submap, delta.dist(sx), 1e-5, k, a = 0, 0.052)
# [1]  -66.86574837 1153.61990432    0.02803073



test.optimization <- function() {
  cat("optimizing with optim\n")
  A <- Fantasio2:::m4_logEmiss( sx@bed, sx@p, sx@submap, 1e-5 )
  R <- NULL
  for(k in 0:(nrow(x.be)-1)) {
    f <- function(x) -Fantasio2:::logLikelihood_gradient( A[k*2 + 1:2, ],  delta.dist(sx), x[1], x[2] )[1]
    grad <- function(x) -Fantasio2:::logLikelihood_gradient( A[k*2 + 1:2, ],  delta.dist(sx), x[1], x[2] )[-1]
  
    o <- optim( c(0.1, 0.1), fn = f, gr = grad, method = "L-BFGS-B", lower = c(1e-10,0), upper = c(Inf, 1) )
    R <- rbind(R, c(o$par, o$value))
  }
  R <- as.data.frame(R)
  names(R) <- c("a", "f", "fx")
  
  cat("optimizing with LFBGSpp\n")
  # S <- NULL
  # for(k in 0:(nrow(x.be)-1)) {
  #   o <- Fantasio2:::testOptimLikelihood( sx@bed, sx@p, sx@submap, delta.dist(sx), 1e-5, k)
  #   S <- rbind(S, o)
  # }
  # S <- as.data.frame(S)
 
  S <- as.data.frame( Fantasio2:::festim(sx@bed, sx@p, sx@submap, delta.dist(sx), 1e-5) )
  par(mfrow=c(1,2))
  plot(R$a, S$a)
  plot(R$f, S$f)
  list(optim = R, fanta = S)
} 

# a = 0.337448, f = 0.0625532

