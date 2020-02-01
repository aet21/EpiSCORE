\docType{data}
\name{dbSCM2}
\alias{dbSCM2}
\title{Database of matched DNAm and gene expression from SCM2}
\format{
  Four separate matrices containing the integrated matched
  DNAm and gene expression datasets, the estimated probabilities
  of expression matrix and a matrix listing the imputable genes.
}
\usage{
data(dbSCM2)
}
\description{
  The integrated and matched DNAm and gene expression datasets are
  of dimension 13290 genes and 34 samples. Original data derives from
  the Stem-Cell Matrix compendium-2 (SCM2). The probability matrix is
  of the same dimension and derives from a 2-state Bayesian Gamma
  mixture-model giving us the probability of each gene belonging to the
  expressed state. A final matrix lists statistics for the subset of
  imputable genes.
}
\details{
Upon loading the data object, we find the following:
\itemize{
  \item \code{intSCM2dnaM.m}: the DNAm data matrix
  \item \code{intSCM2expR.m}: the corresponding matched expression matrix
  \item \code{pEgXscm2.m}: the probability matrix of being expressed
  \item \code{topIntMRscm2.m}: a matrix listing the imputable genes
}
}
\references{
Zhu T, Breeze CE, Beck S, Teschendorff AE
\emph{Cell-type deconvolution of bulk tissue DNA methylomes 
from single-cell RNA-Seq data}
Submitted.
}
\keyword{datasets}
