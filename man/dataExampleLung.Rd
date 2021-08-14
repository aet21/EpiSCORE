\docType{data}
\name{dataExampleLung}
\alias{dataExampleLung}
\title{Data Object Files for Running Vignette}
\format{Contains a number of R-objects for running the example in the vignette
}

\usage{
data(dataExampleLung)
}
\description{
Contains data objects for the simulated dataset and true mixing
fractions, the processed lung squamous cell carcinoma DNAm TCGA dataset,
summarized at the gene promoter level, associated phenotype information
and the lists of lung-cancer associated DMCTs.
}
\details{
Upon loading the data object, we find the following:
\itemize{
  \item \code{avSIM.m}: 100 simulated DNAm mixtures with DNAm values
  summarized at the gene promoter level.
  \item \code{trueW.m}: true mixing fractions for the 100 mixtures.
  \item \code{avLUSCtss.m}: the TCGA LUSC DNAm dataset of 316 samples,
  summarized at the gene-promoter level.
  \item \code{phenoLUSC.lv}: list of phenotype values for the TCGA
  dataset.
}
}
\references{
Teschendorff AE, Zhu T, Breeze CE, Beck S.
\emph{Cell-type deconvolution of bulk tissue DNA methylomes 
from single-cell RNA-Seq data}
Genome Biol.2020
}
\keyword{datasets}
