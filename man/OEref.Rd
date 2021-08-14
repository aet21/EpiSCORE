\docType{data}
\name{OERef}
\alias{OERef}
\title{Olfactory epithelium (OE) mRNA expression and DNAm reference matrix}
\format{Two matrices representing the mRNA and DNAm reference matrices
}

\usage{
data(OEref)
}
\description{
The expression reference matrix is defined over 1889 marker genes and 9
cell-types (Basal, Fib=fibroblast, Gland, Macrophages, Immature Neurons,
Mature Neurons, Pericytes, Plasma-cells, T-cells). The DNAm
reference matrix is defined over 327 marker gene promoters,and the same
9 cell-types, with the last column labeling the
confidence weight.
}
\details{
Upon loading the data object, we find the following:
\itemize{
  \item \code{exprefOE.m}: the expression reference matrix
  \item \code{mrefOE.m}: the DNAm reference matrix
}
}
\references{
}
\keyword{datasets}
