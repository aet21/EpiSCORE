\docType{data}
\name{probeInfoEPICv2}
\alias{probeInfoEPICv2}
\title{Probe information for Illumina EPICv2 DNAm beadarray}
\format{
  A list containing required information for all probes
  on the Illumina EPICv2 (950k) platform.
}
\usage{
data(probeInfoEPICv2)
}
\description{
The list contains entries describing unique probeID, probe design-type,
gene regulatory region (1=TSS1500, 2=TSS200, 3=5UTR, 4=1st Exon, 5=Body,
6=3UTR), gene symbol, NCBI gene identifier (EID) and CpG identifier (not
unique as there are multiple probes mapping to same CpG).
}
\details{
A list object \code{probeInfoEPICv2.lv} with following entries:
\itemize{
  \item \code{probeID:} probe IDs
  \item \code{Design:} probe design type
  \item \code{GeneGroup:} gene regulatory region
  \item \code{Symbol:} gene symbol it maps to
  \item \code{EID:} NCBI (Entrez) gene ID
  \item \code{CpG:} CpG identifier
}
}
\references{
Zhu T, Breeze CE, Beck S, Teschendorff AE
\emph{A pan-tissue DNA methylation atlas enables in silico decomposition
of human tissue methylomes at cell-type resolution.} Nat Methods 2022
Mar;19(3):296-305.
}
\keyword{datasets}
