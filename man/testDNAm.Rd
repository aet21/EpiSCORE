\docType{data}
\name{testDNAm}
\alias{testDNAm}
\title{Simulated in-silico mixture DNAm dataset}
\format{
  A matrix representing 2 simulated mixtures of epithelial, endothelial,
  fibroblast and immune-cell DNAm 450k profiles, with rows labeling CpGs
  and columns labeling mixtures, with mixture weights specified in
  another matrix.
}
\usage{
data(testDNAm)
}
\description{
Relatively pure epithelial, endothelial, fibroblast and immune-cell
samples, all for which Illumina 450k DNAm profiles are available, were
mixed in-silico in known proportions. The DNAm data matrix is defined over
483793 common CpGs and 2 mixture samples. The matrix of mixing
fractions defined over 2 samples and the 4 cell-types.
}
\details{
 Loading the data object, there are two matrices:
\itemize{
  \item \code{testDNAm.m}: The simulated DNAm data matrix.
  \item \code{trueFrac.m}: The true cell-type fractions for all samples.
}
}
\references{
Zhu T, Breeze CE, Beck S, Teschendorff AE
\emph{Cell-type deconvolution of bulk tissue DNA methylomes 
from single-cell RNA-Seq data}
Submitted.
}
\keyword{datasets}
