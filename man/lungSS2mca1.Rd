\docType{data}
\name{lungSS2mca1}
\alias{lungSS2mca1}
\title{Lung scRNA-Seq atlas SmartSeq2 from Mouse Cell Atlas-1}
\format{Contains the normalized SmartSeq2 lung atlas from MCA1 with rows
labeling 15845 genes and columns labeling 1641 cells, an integer cell-type index
vector, and the names of the cell-types.
}
\usage{
data(lungSS2mca1)
}
\description{
The data has been log2-normalized with a pseudocount of 1, so that
zeroes correspond to zero mapped reads. The integer cell-type index
vector annotates each cell/column to a cell-type, with 1,2,3 and 4
indicating epithelial, endothelial, fibroblast and immune-cell,
respectively. An additional vector specifies the cell-types.
}
\details{
Upon loading the data object, we find the following:
\itemize{
  \item \code{lungSS2mca1.m}: the scRNA-Seq data matrix
  \item \code{celltypeSS2.idx}: the integer cell-annotation index vector
  \item \code{celltypeSS2.v}: the names of unique cell-types
}
}
\references{
Teschendorff AE, Zhu T, Breeze CE, Beck S.
\emph{Cell-type deconvolution of bulk tissue DNA methylomes 
from single-cell RNA-Seq data}
Genome Biol.2020
}
\keyword{datasets}
