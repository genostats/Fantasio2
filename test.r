library(Fantasio2)

Fantasio.parameters(use_froh = FALSE)

fxp <- read.bed.matrix("~/COURS/SDS/logiciels/TP-Fantasio/DataXP/afxp")
fxp <- set.stats(fxp)
fxp

freqs <- read.table("~/COURS/SDS/logiciels/TP-Fantasio/DataXP/freqs.txt", header=T)
summary(freqs)

#Run Fantasio
set.seed(123)
F1.fxp.me <- Fantasio(bedmatrix=fxp, allele.freq = freqs$freq.A2.me, n=100)
F1.fxp.me@FLOD_recap[1:5, 1:5]
F1.fxp.me@HBD_recap[1:5, 1:5]
F1.fxp.me@submap_summary

if(file.exists("/tmp/coin.hbd")) file.remove("/tmp/coin.hbd")
if(file.exists("/tmp/coin.flod")) file.remove("/tmp/coin.flod")
set.seed(123)
F1.fxp.me2 <- Fantasio(bedmatrix=fxp, allele.freq = freqs$freq.A2.me, n=100, basename = "/tmp/coin")
F1.fxp.me2@FLOD_recap
F1.fxp.me2@HBD_recap
F1.fxp.me2@submap_summary

# set.seed(1); FF <- Fantasio(bedmatrix=fxp, allele.freq = freqs$freq.A2.me, n=100, dense = FALSE)

#HFLOD genome-wide
HFLOD.fxp.me <- HBD.gwas(F1.fxp.me, phen.code = "plink")
HFLOD.fxp.me2 <- HBD.gwas(F1.fxp.me2, phen.code = "plink")

par(mfrow = c(1,2))
HFLOD.manhattan.plot(HFLOD.fxp.me)
HFLOD.manhattan.plot(HFLOD.fxp.me2)

#HFLOD Chr3
HFLOD.plot.chr(HFLOD.fxp.me, chr=3, MA = FALSE)
HFLOD.plot.chr(HFLOD.fxp.me2, chr=3, MA = FALSE)

#position gene XPC
genepos <- data.frame(chr=3, start=32.9, end=33)
HFLOD.plot.chr(HFLOD.fxp.me, chr=3, MA=FALSE, regions=genepos)
HFLOD.plot.chr(HFLOD.fxp.me2, chr=3, MA=FALSE, regions=genepos)

#plot HBD segments on chr 3
HBD.plot.chr(F1.fxp.me, chr=3, region = genepos)
HBD.plot.chr(F1.fxp.me2, chr=3, region = genepos)

#All patients share an HBD segment except Indiv P1. Why?
#HBD.plot.id(F1.fxp.me, id = "P1", famid = "FXP1")
F1.fxp.me@submap_summary

#No HBD segment evaluation for P1 because of low QUALITY (<95)
#It could be due to the allele frequencies 
#that might not be from the best reference population for that individuala

head(F1.fxp.me@HBD_segments[[1]])
#Indiv P2
HBD.plot.id(F1.fxp.me, id = "P2", famid = "FXP2")
HBD.plot.id(F1.fxp.me2, id = "P2", famid = "FXP2")

######################################

library(Fantasio2)

simu <- read.bed.matrix("~/COURS/SDS/logiciels/TP-Fantasio/Data1006/simu_H1_haplo_chr5-10.bed")

set.seed(123)
F.1006 <- Fantasio(simu, n = 10)
F.1006@HBD_recap[1:5,1:5]

if(file.exists("/tmp/coin.hbd")) file.remove("/tmp/coin.hbd")
if(file.exists("/tmp/coin.flod")) file.remove("/tmp/coin.flod")
set.seed(123)
F.1006.2 <- Fantasio(simu, n = 10, basename = "/tmp/coin")
F.1006.2@HBD_recap

glm.1006 <- HBD.glm(F.1006, expl_var = "FLOD", phen.code = "plink", score = FALSE)
glm.1006.2 <- HBD.glm(F.1006.2, expl_var = "FLOD", phen.code = "plink", score = FALSE)

par(mfrow = c(1,2))
glm.HBD.plot(glm.1006)
glm.HBD.plot(glm.1006.2)

glm.1006.sc <- HBD.glm(F.1006, expl_var = "FLOD", phen.code = "plink", score = TRUE)
glm.1006.sc.2 <- HBD.glm(F.1006.2, expl_var = "FLOD", phen.code = "plink", score = TRUE)

par(mfrow = c(1,2))
glm.HBD.plot(glm.1006.sc)
glm.HBD.plot(glm.1006.sc.2)

set.sed(1); covar <- runif(nrow(simu))

glm.1006.sc.cov <- HBD.glm(F.1006, expl_var = "FLOD", phen.code = "plink", score = TRUE, covar = covar)
glm.1006.sc.cov.2 <- HBD.glm(F.1006.2, expl_var = "FLOD", phen.code = "plink", score = TRUE, covar = covar)

par(mfrow = c(1,2))
glm.HBD.plot(glm.1006.sc.cov)
glm.HBD.plot(glm.1006.sc.cov.2)

