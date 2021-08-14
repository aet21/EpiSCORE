\docType{data}
\name{HeartRef}
\alias{HeartRef}
\title{Heart mRNA expression and DNAm reference matrix}
\format{Two matrices representing the mRNA and DNAm reference matrices
}

\usage{
data(HeartRef)
}
\description{
The expression reference matrix is defined over 1180 marker genes and 5
cell-types (CM=cardiomyocytes, EC=endothelial, FB=fibroblast,
MP=macrophage, SMC=smooth muscle). The DNAm reference matrix is defined
over 207 marker gene promoters,and the same 5 cell-types, with the last column labeling the confidence weight.
}
\details{
Upon loading the data object, we find the following:
\itemize{
  \item \code{exprefHeart.m}: the expression reference matrix
  \item \code{mrefHeart.m}: the DNAm reference matrix
}
}
\references{
}
\keyword{datasets}
