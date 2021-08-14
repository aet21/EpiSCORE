\docType{data}
\name{SkinRef}
\alias{SkinRef}
\title{Skin mRNA expression and DNAm reference matrix}
\format{Two matrices representing the mRNA and DNAm reference matrices
}

\usage{
data(SkinRef)
}
\description{
The expression reference matrix is defined over 652 marker genes and 7
cell-types (EC=endothelial, Fib=fibroblast, Kera_diff=differentiated
keratinocytes, Kera_undiff=undifferentiated keratinocytes,
Macro=macrophages, Mela=melanocytes, Tcell=T-cells). The DNAm
reference matrix is defined over 151 marker gene promoters,and the same
7 cell-types, with the last column labeling the
confidence weight.
}
\details{
Upon loading the data object, we find the following:
\itemize{
  \item \code{exprefSkin.m}: the expression reference matrix
  \item \code{mrefSkin.m}: the DNAm reference matrix
}
}
\references{
}
\keyword{datasets}
