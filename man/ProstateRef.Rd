\docType{data}
\name{ProstateRef}
\alias{ProstateRef}
\title{Prostate mRNA expression and DNAm reference matrix}
\format{Two matrices representing the mRNA and DNAm reference matrices
}

\usage{
data(ProstateRef)
}
\description{
The expression reference matrix is defined over 680 marker genes and 6
cell-types (BE=basal epithelium, EC=endothelial, Fib=fibroblast,
LE=luminal epithelium, Leu=leukocytes, SM=smooth muscle). The DNAm
reference matrix is defined over 163 marker gene promoters,and the same
6 cell-types, with the last column labeling the
confidence weight.
}
\details{
Upon loading the data object, we find the following:
\itemize{
  \item \code{exprefProstate.m}: the expression reference matrix
  \item \code{mrefProstate.m}: the DNAm reference matrix
}
}
\references{
}
\keyword{datasets}
