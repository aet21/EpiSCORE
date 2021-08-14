\docType{data}
\name{ColonRef}
\alias{ColonRef}
\title{Colon mRNA expression and DNAm reference matrix}
\format{Two matrices representing the mRNA and DNAm reference matrices
}

\usage{
data(ColonRef)
}
\description{
The expression reference matrix is defined over 1124 marker genes and 6
cell-types (Enteroendocrine, Enterocyte, Goblet, immune cell, Paneth,
Progenitor). The DNAm reference matrix is defined over 185 marker gene
promoters,and the same 6 cell-types, with the last column labeling the confidence weight.
}
\details{
Upon loading the data object, we find the following:
\itemize{
  \item \code{exprefColon.m}: the expression reference matrix
  \item \code{mrefColon.m}: the DNAm reference matrix
}
}
\references{
}
\keyword{datasets}
