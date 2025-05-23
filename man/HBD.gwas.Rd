% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/HBD.gwas.R
\name{HBD.gwas}
\alias{HBD.gwas}
\title{Computation of HFLOD scores for HBD-GWAS}
\usage{
HBD.gwas(atlas, phen, phen.code = c("plink", "R"))
}
\arguments{
\item{atlas}{a atlas object}

\item{phen}{the phenotype (default is the bed.matrix phenotype)}

\item{phen.code}{phenotype coding :
       - 'R' : 0:control ; 1:case ; NA:unknown (default)
       - 'plink' : 1:control ; 2:case ; 0/-9/NA:unknown
if 'plink' the function automatically convert it to 'R' to run logistic regression description}
}
\value{
the atlas object with its slot HFLOD completed
}
\description{
This function is used to compute HFLOD scores on individuals in a sample for the HBD-GWAS method
}
