% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ImputeDNAmRef.R
\name{ImputeDNAmRef}
\alias{ImputeDNAmRef}
\title{Impute a DNAm reference matrix}
\usage{
ImputeDNAmRef(
  refexp.m,
  db = c("SCM2", "RMAP"),
  geneID = c("SYMBOL", "ENTREZID")
)
}
\arguments{
\item{refexp.m}{The gene expression reference matrix with rows labeling genes and columns
labeling cell-types. If gene idgentifier is gene symbol this will be converted
to Entrez gene ID.}

\item{db}{A character specifying the database of matched expression and DNAm data.
This has to be either SCM2 for Stem-Cell-Matrix Compendium-2
or RMAP for Epigenomics Roadmap.}

\item{geneID}{A character specifying gene identifier used in expression reference matrix input. 
Either SYMBOL or ENTREZ GENE ID.}
}
\value{
The imputed DNAm reference matrix defined over the same cell-types
as defined in the input expression reference matrix.
}
\description{
This function takes as input a mRNA expression reference matrix and will generate
an imputed DNAm reference using one of two databases (NIH Epigenomics Roadmap
or Stem-Cell Matrix Compendium) specified by the user.
}
\examples{
data(lungSS2mca1)
out.l <- ConstExpRef(lungSS2mca1.m,celltypeSS2.idx,celltypeSS2.v,markspecTH.v=rep(3,4));
refDNAm.m <- ImputeDNAmRef(out.l$ref$med,db="SCM2",geneID="SYMBOL");
print(head(refDNAm.m));



}
\references{
Teschendorff AE, Zhu T, Breeze CE, Beck S.
\emph{Cell-type deconvolution of bulk tissue DNA methylomes 
from single-cell RNA-Seq data}
Genome Biol.2020

Zhu T, Liu J, Beck S, Pan S, Capper D, Lechner M, Thirlwell C, Breeze CE, Teschendorff AE.
\emph{A pan-tissue DNA methylation atlas enables in silico decomposition of human tissue methylomes at cell-type resolution}
Nat Methods 2022
}
