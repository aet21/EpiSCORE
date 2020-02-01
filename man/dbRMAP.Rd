\docType{data}
\name{dbRMAP}
\alias{dbRMAP}
\title{Database of matched DNAm and gene expression from RMAP}
\format{
  Four separate matrices containing the integrated matched
  DNAm and gene expression datasets, the estimated probabilities
  of expression matrix and a matrix listing the imputable genes.
}
\usage{
data(dbRMAP)
}
\description{
  The integrated and matched DNAm and gene expression datasets are
  of dimension 18651 genes and 45 samples. Original data derives from
  the NIH Epigenomics Roadmap (RMAP). The probability matrix is
  of the same dimension and derives from a 2-state Bayesian Gamma
  mixture-model giving us the probability of each gene belonging to the
  expressed state. A final matrix lists statistics for the subset of
  imputable genes.
}
\details{
Upon loading the data object, we find the following:
\itemize{
  \item \code{betaRMAPint.m}: the DNAm data matrix
  \item \code{expRMAPint.m}: the corresponding matched expression matrix
  \item \code{pEgXrmap.m}: the probability matrix of being expressed
  \item \code{topIntMRrmap.m}: matrix listing the imputable genes
}
}
\references{
Zhu T, Breeze CE, Beck S, Teschendorff AE
\emph{Cell-type deconvolution of bulk tissue DNA methylomes 
from single-cell RNA-Seq data}
Submitted.
}
\keyword{datasets}
