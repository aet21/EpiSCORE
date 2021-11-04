\name{ColonRef}
\alias{ColonRef}
\docType{data}
\title{
EpiSCORE References for Colon
}
\description{
The data contains a 381-marker expression and a 104-marker DNAm
reference matrices for colon, which can be used to estimate cell type
fractions of endothelial cells, epithelial cells, lymphocytes, myeloid
cells and stromal cells in colon tissue with EpiSCORE framework. 
}
\usage{data("ColonRef")}
\format{
  Two reference matrices. 
  \describe{
    \item{\code{Colon_Expref.m}}{A scRNAseq expression reference for Colon. Each row is a marker gene (in gene symbol) and each column is a cell type. For each marker gene, the value of the reference matrix is the median expression level for each type of cells in scRNAseq data.}
    \item{\code{Colon_Mref.m}}{A DNAm reference for Colon. Each row is a marker gene (in entrezID) and each column is a cell type. The last column is the weight for each marker. The value of the matrix is the imputed DNAm level from expression reference.}
  }
}
\details{
The EpiSCORE colon expression reference matrix is derived from a
scRNA-Seq study (Elmentaite R et al Nature 2021). Details can be found
in published paper. The expression reference matrix was constructed
running EpiSCORE with marker specificity scores MSS= (4,4,3,4,4) for
endothelial, epithelial, lymphocytes, myeloid and stromal cells,
respectively. The expression reference matrix was validated in HCL colon
and rectum data. The corresponding DNAm reference was imputed from expression reference with EpiSCORE. The weights in the last column indicates the ability for a marker to discriminate different cell types, which can be implemented in the weighted robust partial correlation (wRPC) algorithm in the EpiSCORE package to estimate cell type fractions in bulk colon DNAm data.
}

\references{
Elmentaite et al
\emph{Cells of the human intestinal tract mapped across space and time.} 
Nature 2021

Han X, Zhou Z, Fei L, Sun H, Wang R, Chen Y, et al. 
\emph{Construction of a human cell landscape at single-cell level.} 
Nature. 2020

Teschendorff AE, Zhu T, Breeze CE, Beck S.
\emph{Cell-type deconvolution of bulk tissue DNA methylomes from single-cell RNA-Seq data.}
Genome Biol. 2020
}

\keyword{datasets}