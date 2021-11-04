\docType{data}
\name{EsoRef}
\alias{EsoRef}
\title{Esophageal mRNA expression and DNAm reference matrix}
\format{Two matrices representing the mRNA and DNAm reference matrices
}

\usage{
data(EsoRef)
}
\description{
The expression reference matrix is defined over 1227 marker genes and 8
cell-types (EC=endothelial, basal, stratified, suprabasal and upper
epithelium, fibroblasts, glandular and immune cells. The DNAm reference
matrix is defined over 161 marker gene promoters, and the same 8
cell-types, with the last column labeling the confidence weight.
}
\details{
Upon loading the data object, we find the following:
\itemize{
  \item \code{Eso_Expref.m}: the expression reference matrix
  \item \code{Eso_Mref.m}: the DNAm reference matrix
}
}
\references{
}
\keyword{datasets}
