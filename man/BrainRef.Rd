\docType{data}
\name{BrainRef}
\alias{BrainRef}
\title{Brain mRNA expression and DNAm reference matrix}
\format{Two matrices representing the mRNA and DNAm reference matrices
}

\usage{
data(BrainRef)
}
\description{
The expression reference matrix is defined over 1243 marker genes and 6
cell-types (Neurons, Oligodendrocytes, Astrocytes, OligoPrecursorCells,
Endothelial Cells, Microglia). The DNAm reference matrix is defined over
204 marker gene promoters,and the same 6 cell-types, with the last
column labeling the confidence weight.
}
\details{
Upon loading the data object, we find the following:
\itemize{
  \item \code{exprefBrain.m}: the expression reference matrix
  \item \code{mrefBrain.m}: the DNAm reference matrix
}
}
\references{
}
\keyword{datasets}
