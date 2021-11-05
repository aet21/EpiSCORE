\docType{data}
\name{probeInfo450k}
\alias{probeInfo450k}
\title{Probe information for Illumina 450k DNAm beadarray}
\format{
  A list containing required information for all probes
  on the Illumina 450k platform.
}
\usage{
data(probeInfo450k)
}
\description{
The list contains entries describing the type of cytosine (CpG vs
non-CpG), design-type of the probe (type-1 vs type-2), gene regulatory
region (1=TSS1500, 2=TSS200, 3=5UTR, 4=1st Exon, 5=Body, 6=3UTR),
if in a CpG island or not, probe-ID and Entrez gene ID.
}
\details{
  A list \code{probeInfo450k.lv} with following entries:
\itemize{
  \item \code{typeC}: CpG or non-CpG
  \item \code{Design}: probe design type
  \item \code{GeneGroup}: gene regulatory region
  \item \code{CGI}: CpGisland or not
  \item \code{probeID}: probe ID
  \item \code{EID}: Entrez gene ID
}
}
\references{
Zhu T, Breeze CE, Beck S, Teschendorff AE
\emph{Cell-type deconvolution of bulk tissue DNA methylomes 
from single-cell RNA-Seq data}
Submitted.
}
\keyword{datasets}
