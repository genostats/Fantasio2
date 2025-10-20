# ---- installation and loading procedure ----

# For gaston 1.6.1 (not yet on CRAN...)
if(require(gaston)) {
  gaston_161 <- packageVersion("gaston") >= "1.6.1"
} else {
  gaston_161 <- FALSE
}

if(!gaston_161) {
  install.packages(c("Rcpp", "RcppParallel, RcppEigen"))
  remotes::install_github("genostats/gaston")
}

# For Fantasio2
if(!require(Fantasio2)) {
  remotes::install_github("genostats/Fantasio2")
}

# For the genetic map
if(!require(HumanGeneticMap)) {
  install.packages("HumanGeneticMap", repos = "https://genostats.github.io/R/")
}

# ---- example of data analysis: consanguinity -----

# set number of threads and verbosity
Fantasio.parameters(n_threads = 8, verbose = TRUE)

# read example data
bm <- read.bed.matrix(system.file("extdata", "test_Fantasio.bed", package = "Fantasio2"))

# NOTE : this data set has only autosomes. 
# In general, to remove X, Y and mt, do:
bm <- select.snps(bm, is.autosome(chr))

# set distance in cM in bm@snps$dist (change b37 to b38 if needed !!)
bm <- set.dist(bm, HumanGeneticMap::genetic.map.b37)
                      
# estimate consanguinity (n = 100 submaps is the default, one can try with 20 to begin with)
F <- Fantasio(bm, n = 100)

# look at results in F@submap_summary !!!
# how many inbred individuals?
sum(F@submap_summary$inbred)


# ----- example of data analysis: association --------

# load covariates
cov <- readRDS(system.file("extdata", "cov.rds", package="Fantasio2"))

# HBD-GWAS method on phenotype included in bed.matrix
hbd.gwas <- HBD.gwas(F, phen.code = "plink")

# HBD-GWAS with another binary phenotype
phen2 <- sample(c(1,2), size=length(bm@ped$pheno), replace=TRUE)
hbd.gwas2 <- HBD.gwas(F, phen = phen2, phen.code = "plink")


# logistic regression bed.matrix phenotype ~ FLOD
assoc.FLOD <- HBD.glm(F, phen.code = "plink")

# logistic regression bed.matrix phenotype ~ pHBD
assoc.pHBD <- HBD.glm(F, expl_var = "pHBD", phen.code = "plink")

# logistic regression phen2 ~ FLOD
assoc.FLOD2 <- HBD.glm(F, phen = phen2, phen.code = "plink")

# logistic regression bed.matrix phenotype ~ FLOD+cov
assoc.FLOD <- HBD.glm(F, covar_df = cov, phen.code = "plink")

# Manhattan plot for bed.matrix phenotype ~ FLOD, bilateral
glm.HBD.plot(assoc.FLOD)

# Manhattan plot for bed.matrix phenotype ~ FLOD, right
glm.HBD.plot(assoc.FLOD, test="right")
