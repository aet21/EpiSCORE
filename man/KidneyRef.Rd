\name{KidneyRef}
\alias{KidneyRef}
\docType{data}
\title{
EpiSCORE References for Kidney
}
\description{
The data contains a 272-marker expression and a 32-marker DNAm reference
matrix for kidney, which can be used to estimate cell type fractions of
endothelial cells, epithelial cells, fibroblasts and immune cells withh EpiSCORE framework. 
}
\usage{data("KidneyRef")}
\format{
  Two reference matrices. 
  \describe{
    \item{\code{Kidney_Expref.m}}{A scRNAseq expression reference for Kidney. Each row is a marker gene (in gene symbol) and each column is a cell type. For each marker gene, the value of the reference matrix is the median expression level for each type of cells in scRNAseq data.}
    \item{\code{Kidney_Mref.m}}{A DNAm reference for Kidney. Each row is a marker gene (in entrezID) and each column is a cell type. The last column is the weight for each marker. The value of the matrix is the imputed DNAm level from expression reference.}
  }
}
\details{
The EpiSCORE kidney expression reference matrix is derived from a human
kidney 10X scRNA-seq dataset from Muto et al. The expression reference
matrix was constructed running EpiSCORE with marker specificity scores
MSS= (3,2,2,2) for endothelial cells, epithelial cells, fibroblasts and
immune cells. The expression reference matrix was validated in a human
kidney 10X dataset from Wu et al. The corresponding DNAm reference was imputed from expression reference with EpiSCORE. The weights in the last column indicates the ability for a marker to discriminate different cell types, which can be implemented in the weighted robust partial correlation (wRPC) algorithm in the EpiSCORE package to estimate cell type fractions in bulk kidney DNAm data. 
}

\references{
Muto Y et al
\emph{Single cell transcriptional and chromatin accessibility profiling
  redefine cellular heterogeneity in the adult kidney.}
Nat Commun. 2021

Wu H, Malone AF, Donnelly EL, Kirita Y, Uchimura K, Ramakrishnan SM, et al. 
\emph{Single-Cell Transcriptomics of a Human Kidney Allograft Biopsy Specimen Defines a Diverse Inflammatory Response.} 
J Am Soc Nephrol. 2018

Teschendorff AE, Zhu T, Breeze CE, Beck S.
\emph{Cell-type deconvolution of bulk tissue DNA methylomes from single-cell RNA-Seq data.}
Genome Biol. 2020
}

\keyword{datasets}