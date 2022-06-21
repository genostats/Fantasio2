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
Fantasio.parameters(n_threads = 4)
Fantasio.parameters(use_float = TRUE)
debug(Fantasio2:::recap.HBD.FLOD)

set.seed(1); a <- Fantasio(x.be, n = 5, phen.code = "plink")


