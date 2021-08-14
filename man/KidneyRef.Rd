\docType{data}
\name{KidneyRef}
\alias{KidneyRef}
\title{Kidney mRNA expression and DNAm reference matrix}
\format{Two matrices representing the mRNA and DNAm reference matrices
}

\usage{
data(KidneyRef)
}
\description{
The expression reference matrix is defined over 289 marker genes and 6
cell-types (CD=columnar ductal, EC=endothelial, Fib=fibroblast,
IC=immune cell, LOH=loop of Henle, PT=proximal tubule). The DNAm
reference matrix is defined over 74 marker gene promoters,and the same 6 cell-types, with the last column labeling the confidence weight.
}
\details{
Upon loading the data object, we find the following:
\itemize{
  \item \code{exprefKidney.m}: the expression reference matrix
  \item \code{mrefKidney.m}: the DNAm reference matrix
}
}
\references{
}
\keyword{datasets}
