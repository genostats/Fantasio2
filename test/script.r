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
  # x.be <- x.be[1:47,]
}

require(Fantasio2)

if(!exists("s")) {
  s <- segments.list.by.hotspots(x.be)
}

# set.seed(4); submap <- Fantasio2:::rsubmap(s)
set.seed(1); submap <- Fantasio2:::rsubmap(s)
d.dist <- Fantasio2:::delta.dist(x.be, submap)

test.log.emission <- function() {
  cat("test log emission\n")
  A <- Fantasio2:::m4_logEmiss( x.be@bed, x.be@p, submap, 1e-5 )
  for(k in 0:(nrow(x.be)-1)) {
    B <- Fantasio2:::testLogEmiss( x.be@bed, x.be@p, submap, 1e-5, k )
    if(!all( A[ k*2 + 1:2, ] == B)) stop(k)
  }
}


test.log.likelihood <- function(a = 0.1, f = 0.1, epsilon = 1e-5) {
  cat("test log likelihood\n")
  A <- Fantasio2:::m4_logEmiss( x.be@bed, x.be@p, submap, 1e-5 )
  new <- old <- NULL
  for(k in 0:(nrow(x.be)-1)) {
    R1 <- Fantasio2:::testLikelihood( x.be@bed, x.be@p, submap, d.dist, 1e-5, k, a, f)
    R2 <- Fantasio2:::logLikelihood_gradient( A[k*2 + 1:2, ],  d.dist, a, f)
 # R2 <- .Call('_Fantasio_logLikelihood_gradient', PACKAGE = "Fantasio", A[k*2 + 1:2, ],  d.dist, a, f)
    new <- rbind(new, R1)
    old <- rbind(old, R2)
  }
  cat("Range des différences :", range(new - old), "\n")
  list(new = new, old = old)
}  
 
test.forward.backward <- function() { 
  cat("test forward backward\n")
  S <- as.data.frame( Fantasio2:::festim(x.be@bed, x.be@p, submap, d.dist, 1e-5) )
  
  A <- Fantasio2:::m4_logEmiss( x.be@bed, x.be@p, submap, 1e-5 )
  R <- NULL
  for(k in 0:(nrow(x.be)-1)) {
    R <- cbind(R, Fantasio2:::forward_backward( A[k*2 + 1:2, ],  d.dist, S$a[k+1], S$f[k+1] )[2,]) 
  }

  T <- Fantasio2:::probaHBD(x.be@bed, x.be@p, submap, d.dist, rep(TRUE, nrow(x.be)), S$a, S$f, 1e-5)
  list(new = T, old = R)
}


## Il y a un souci de discontinuité du gradient en a = 0
## -> revoir les calculs papier...

k <- 17
Fantasio2:::testLikelihood( x.be@bed, x.be@p, submap, d.dist, 1e-5, k, a = 1e-22, 0.052)
# [1] -66.86574837 -51.74281869   0.02803073
Fantasio2:::testLikelihood( x.be@bed, x.be@p, submap, d.dist, 1e-5, k, a = 0, 0.052)
# [1]  -66.86574837 1153.61990432    0.02803073



test.optimization <- function(epsilon = 1e-5) {
  cat("optimizing with optim\n")
  A <- Fantasio2:::m4_logEmiss( x.be@bed, x.be@p, submap, epsilon )
  R <- NULL
  for(k in 0:(nrow(x.be)-1)) {
    f <- function(x) -Fantasio2:::logLikelihood_gradient( A[k*2 + 1:2, ],  d.dist, x[1], x[2] )[1]
    grad <- function(x) -Fantasio2:::logLikelihood_gradient( A[k*2 + 1:2, ],  d.dist, x[1], x[2] )[-1]
  
    o <- optim( c(0.05, 0.05), fn = f, gr = grad, method = "L-BFGS-B", lower = c(1e-2,0), upper = c(Inf, 0.999) )
    R <- rbind(R, c(o$par, o$value))
  }
  R <- as.data.frame(R)
  names(R) <- c("a", "f", "fx")
  
  cat("optimizing with LFBGSpp\n")
  # S <- NULL
  # for(k in 0:(nrow(x.be)-1)) {
  #   o <- Fantasio2:::testOptimLikelihood( x.be@bed, x.be@p, submap, d.dist, 1e-5, k)
  #   S <- rbind(S, o)
  # }
  # S <- as.data.frame(S)
 
  S <- as.data.frame( Fantasio2:::festim(x.be@bed, x.be@p, submap, d.dist, epsilon) )
  par(mfrow=c(1,2))
  plot(R$a, S$a)
  plot(R$f, S$f)
  list(optim = R, fanta = S)
} 

# a = 0.337448, f = 0.0625532

