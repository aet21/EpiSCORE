\docType{data}
\name{BreastRef}
\alias{BreastRef}
\title{Breast mRNA expression and DNAm reference matrix}
\format{Two matrices representing the mRNA and DNAm reference matrices
}

\usage{
data(BreastRef)
}
\description{
The expression reference matrix is defined over 2907 marker genes and 7
cell-types (basal, endothelial, fat, fibroblast, luminal, lymphocyte,
macrophage). The DNAm reference matrix is defined over 468 marker gene promoters,and the same 7 cell-types, with the last
column labeling the confidence weight.
}
\details{
Upon loading the data object, we find the following:
\itemize{
  \item \code{exprefBreast.m}: the expression reference matrix
  \item \code{mrefBreast.m}: the DNAm reference matrix
}
}
\references{
}
\keyword{datasets}
