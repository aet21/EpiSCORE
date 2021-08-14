\docType{data}
\name{PancreasRef}
\alias{PancreasRef}
\title{Pancreas mRNA expression and DNAm reference matrix}
\format{Two matrices representing the mRNA and DNAm reference matrices
}

\usage{
data(PancreasRef)
}
\description{
The expression reference matrix is defined over 822 marker genes and 6
cell-types (acinar, endocrine, ductal, endothelial, immune-cell,stellate). The DNAm reference matrix is defined over 178 marker gene promoters,and the same 6 cell-types, with the last column labeling the confidence weight.
}
\details{
Upon loading the data object, we find the following:
\itemize{
  \item \code{exprefPancreas.m}: the expression reference matrix
  \item \code{mrefPancreas.m}: the DNAm reference matrix
}
}
\references{
}
\keyword{datasets}
