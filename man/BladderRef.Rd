\docType{data}
\name{BladderRef}
\alias{BladderRef}
\title{Bladder mRNA expression and DNAm reference matrix}
\format{Two matrices representing the mRNA and DNAm reference matrices
}

\usage{
data(BladderRef)
}
\description{
The expression reference matrix is defined over 522 marker genes and 4
cell-types (EC=endothelial, Epi=epithelial, Fib=fibroblast, IC=immune
cell). The DNAm reference matrix is defined over 163 marker gene
promoters,and the same 4 cell-types, with the last column labeling the
confidence weight.
}
\details{
Upon loading the data object, we find the following:
\itemize{
  \item \code{exprefBladder.m}: the expression reference matrix
  \item \code{mrefBladder.m}: the DNAm reference matrix
}
}
\references{
}
\keyword{datasets}
