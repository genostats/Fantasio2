% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/HBD.glm.r
\name{HBD.glm}
\alias{HBD.glm}
\title{Logistic regression on HBD probability or FLOD score}
\usage{
HBD.glm(
  x,
  expl_var = c("FLOD", "pHBD"),
  phen,
  covar_df,
  covar,
  phen.code = c("plink", "R"),
  threshold
)
}
\arguments{
\item{x}{an atlas object}

\item{expl_var}{the explanatory variable 'FLOD' or 'pHBD'}

\item{phen}{the phenotype (default is the bed.matrix phenotype)}

\item{covar_df}{a dataframe or a matrix containing covariates}

\item{covar}{covariates of interest such as 'age', 'sex' , ...
if missing, all covariates of the dataframe are considered}

\item{phen.code}{phenotype coding :
       - 'R' : 0:control ; 1:case ; NA:unknown (default)
       - 'plink' : 1:control ; 2:case ; 0/-9/NA:unknown
if 'plink' the function automatically convert it to 'R' to run logistic regression}

\item{threshold}{if missing, no thresholding. If set, the explanatory variable will be dichotomized accordingly. Suggested values: 0 for FLOD, 0.5 for pHBD.}

\item{n.cores}{number of cores for parallelization calculation (default = 1)}
}
\description{
Logistic regression on HBD probability or FLOD score
}
