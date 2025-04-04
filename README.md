---
title: "Introduction to EpiSCORE and DNAm-Atlas"
author:
- name: "Andrew E. Teschendorff"
  affiliation: 
  - CAS Key Lab of Computational Biology, PICB, SINH
date: "2025-04-03"
package: EpiSCORE
output:
  BiocStyle::html_document:
    toc_float: true
---

# Summary

EpiSCORE is an R-package for constructing a tissue-specific DNA methylation reference matrix that can be subsequently used in conjunction with a reference-based cell-type deconvolution algorithm to (i) obtain cell-type fraction estimates in a corresponding bulk-tissue sample for which a genome-wide DNAm profile exists, and (ii) to infer cell-type specific differential DNA methylation signals in the context of a general Epigenome-Wide-Association Study. EpiScore is aimed at complex solid tissues, for which experimentally generating appropriate DNAm reference matrices representing all the major cell-types within the tissue is challenging. EpiScore exploits the tissue-specific single-cell RNA-Sequencing atlases to impute corresponding tissue-specific DNA methylation references. The EpiSCORE R-package currently contains DNAm reference matrices for 13 solid tissues. Of note, the EpiSCORE DNAm reference panels are defined at the level of gene promoters. As such, the method is applicable to any DNAm-assay where DNAm can be summarized as beta-values at the level of gene promoters. It is therefore applicable to all Illumina arrays, WGBS or RRBS.

# Installation

To install:

```r
library(devtools)
devtools::install_github("aet21/EpiSCORE")
```

# References

Teschendorff AE, Zhu T, Breeze CE, Beck S. EPISCORE: cell type deconvolution of bulk tissue DNA methylomes from single-cell RNA-Seq data. Genome Biology 2020 Sep 4;21(1):221. doi: 10.1186/s13059-020-02126-9 .

Zhu T, Liu J, Beck S, Pan S, Capper D, Lechner M, Thirlwell C, Breeze CE, Teschendorff AE. A pan-tissue DNA methylation atlas enables in silico decomposition of human tissue methylomes at cell-type resolution. Nat Methods 2022 Mar;19(3):296. doi: 10.1038/s41592-022-01412-7 .