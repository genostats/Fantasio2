require(gaston)
require(HGDP.CEPH)

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

# pour avoir les memes bornes pour l'optimisation que dans la v1
# MARGOT : essaie de commenter cette ligne pour voir la différence
Fantasio2::Fantasio.parameters(lower = c(0.01, 0), upper = c(Inf, 0.999))

## Avec la même graine on obtient la même submap, j'ai vérifié
set.seed(1); F1 <- Fantasio::Fantasio(x.be, "Hotspots", list(hotspots = Fantasio::hotspot_hg19), n = 1, verbose = FALSE)
set.seed(1); F2 <- Fantasio2::Fantasio(x.be, list(hotspots = Fantasio::hotspot_hg19), n = 1, epsilon = 1e-3)
# MARGOT : attention à ce epsilon = 1e-3 ... la valeur par défaut 1e-5 différente de la v2 est source de différences également


#### Récupération de a, f, likelihood ####

## -- Pour la V1 c'est facile toutes les valeurs sont conservées --
R1 <- data.frame(a = F1@submaps_list[[1]]@a, f = F1@submaps_list[[1]]@f, likelihood = F1@submaps_list[[1]]@likelihood1)

## -- Pour la V2 on n'a que a et f il faut refaire le calcul de la log vraisemblance --

# les a f likelihood 
R2 <- as.data.frame(F2@estimations)[,1:2]  # a et f

# calcul de la vraisemblance...
submap <- F1@submaps_list[[1]]@submap
d.dist <- Fantasio2:::delta.dist(x.be, submap)
R2$likelihood <- 0
for(k in 1:nrow(R2)) {
  R2$likelihood[k] <- Fantasio2:::testLikelihood(x.be@bed, x.be@p, submap, d.dist, epsilon = 1e-3, i = k-1, a = R2$a[k], f = R2$f[k])[1]
}

#### Fin récupération ###

# un petit graphe pour voir
plot(R2$likelihood - R1$likelihood)
abline(h=0)


