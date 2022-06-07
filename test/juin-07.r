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
  x.be <- x.be[1:47,]
}

if(!exists("segment.list")) {
  segment.list <- segments.list.by.hotspots(x.be)
}


set.seed(4)
sub <- Fantasio2:::rsubmap(segment.list)
d.dist <- Fantasio2:::delta.dist.0(x.be, sub)

stop()

n.submaps <- 100
Seeds <- matrix( nrow = length( get(".Random.seed", envir = .GlobalEnv) ), ncol = n.submaps )
A <- matrix( nrow = nrow(x.be), ncol = n.submaps )
F <- matrix( nrow = nrow(x.be), ncol = n.submaps )


for(i in 4) { #1:n.submaps) {
set.seed(i)
cat(i, "\n")
  Seeds[,i] <- get(".Random.seed", envir = .GlobalEnv)
  sub <- Fantasio2:::rsubmap(segment.list)
  d.dist <- Fantasio2:::delta.dist.0(x.be, sub)
  S <- Fantasio2:::festim(x.be@bed, x.be@p, sub, d.dist, 1e-5)
  A[,i] <- S$a
  F[,i] <- S$f
}
