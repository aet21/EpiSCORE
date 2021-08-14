\docType{data}
\name{LungRef}
\alias{LungRef}
\title{Lung mRNA expression and DNAm reference matrix}
\format{Two matrices representing the mRNA and DNAm reference matrices
}

\usage{
data(LungRef)
}
\description{
The expression reference matrix is defined over 1788 marker genes and 7
cell-types (endothelial, epithelial, granulocyte, lymphocyte,
macrophage, monocyte, stromal). The DNAm reference matrix is defined over 369 marker gene promoters,and the same 7 cell-types, with the last column labeling the confidence weight.
}
\details{
Upon loading the data object, we find the following:
\itemize{
  \item \code{exprefLung.m}: the expression reference matrix
  \item \code{mrefLung.m}: the DNAm reference matrix
}
}
\references{
}
\keyword{datasets}
