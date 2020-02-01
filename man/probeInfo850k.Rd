\docType{data}
\name{probeInfo850k}
\alias{probeInfo850k}
\title{Probe information for Illumina 850k DNAm beadarray}
\format{
  A list containing required information for all probes
  on the Illumina 850k platform.
}
\usage{
data(probeInfo850k)
}
\description{
The list contains entries describing probeID, probe design-type, gene regulatory
region (1=TSS1500, 2=TSS200, 3=5UTR, 4=1st Exon, 5=Body, 6=3UTR), if in
a CpG island or not, the type of cytosine (CpG vs non-CpG), and Entrez
gene ID.
}
\details{
A list object \code{probeInfoALL.lv} with following entries:
\itemize{
  \item \code{probeID:} probe IDs
  \item \code{design:} probe design type
  \item \code{GeneGroup:} gene regulatory region
  \item \code{CGIgroup:} relation to CpG island
  \item \code{CG:} CpG or non-CpG
  \item \code{EID:} Entrez gene ID
}
}
\references{
Zhu T, Breeze CE, Beck S, Teschendorff AE
\emph{Cell-type deconvolution of bulk tissue DNA methylomes 
from single-cell RNA-Seq data}
Submitted.
}
\keyword{datasets}
