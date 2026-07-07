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


#HFLOD genome-wide
HFLOD.fxp.me <- HBD.gwas(F1.fxp.me, phen.code = "plink")
HFLOD.manhattan.plot(HFLOD.fxp.me)

#HFLOD Chr3
HFLOD.plot.chr(HFLOD.fxp.me, chr=3, MA = FALSE)

#position gene XPC
genepos <- data.frame(chr=3, start=32.9, end=33)
HFLOD.plot.chr(HFLOD.fxp.me, chr=3, MA=FALSE, regions=genepos)

#plot HBD segments on chr 3
HBD.plot.chr(F1.fxp.me, chr=3, region = genepos)

#All patients share an HBD segment except Indiv P1. Why?
HBD.plot.id(F1.fxp.me, id = "P1", famid = "FXP1")
F1.fxp.me@submap_summary

#No HBD segment evaluation for P1 because of low QUALITY (<95)
#It could be due to the allele frequencies 
#that might not be from the best reference population for that individuala

head(F1.fxp.me@HBD_segments[[1]])
#Indiv P2
HBD.plot.id(F1.fxp.me, id = "P2", famid = "FXP2")

#Change Quality threshold to see HBD segments of P1 individual
set.seed(123)
F1.fxp.me.Q50 <- Fantasio(bedmatrix=fxp, allele.freq = freqs$freq.A2.me, n=100, min.quality = 50)
F1.fxp.me.Q50@submap_summary

#P1 does not seem to have HBD segments compatible with having ~6% of his/her genome HBD (1C)
HBD.plot.id(F1.fxp.me.Q50, id = "P1", famid = "FXP1")
HBD.plot.chr(F1.fxp.me.Q50, chr=3)
genepos <- data.frame(chr=3, start=32.9, end=33)

HFLOD.fxp.me.Q50 <- HBD.gwas(F1.fxp.me.Q50, phen.code = "plink")
HFLOD.plot.chr(HFLOD.fxp.me.Q50, chr=3, MA=FALSE, regions=genepos)

#2.3# Change allele frequencies to HGDP-CEPH European (pre-calculated)
set.seed(123)
F1.fxp.eur <- Fantasio(bedmatrix=fxp, allele.freq = freqs$freq.A2.eur, n=100)
F1.fxp.eur@submap_summary

#Plot inbreeding F
#Middle East allele frequencies vs Sample frequencies
plot(F1.fxp@submap_summary$f_median, F1.fxp.me@submap_summary$f_median)
abline(0,1,col="red")
points(F1.fxp@submap_summary$f_median[9], F1.fxp.me@submap_summary$f_median[9], pch=16) #FXP3  P3 (rond)
points(F1.fxp@submap_summary$f_median[3], F1.fxp.me@submap_summary$f_median[3], pch=15) #FXP1  P1 (square)
#legend

#Middle East allele freq vs European
plot(F1.fxp.eur@submap_summary$f_median, F1.fxp.me@submap_summary$f_median)
abline(0,1,col="red")
points(F1.fxp.eur@submap_summary$f_median[9], F1.fxp.me@submap_summary$f_median[9], pch=16) #FXP3  P3
points(F1.fxp.eur@submap_summary$f_median[3], F1.fxp.me@submap_summary$f_median[3], pch=15) #FXP1  P1
#legend


