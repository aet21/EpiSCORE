#' @title 
#' Construct a DNAm matrix at gene-level
#' 
#' @aliases constAvBetaTSS
#'  
#' @description 
#' This function takes as input a DNAm data matrix generated with
#' Illumina 450k or 850k technology and outputs a DNAm data matrix
#' at the gene level, by computing the average DNAm over a window
#' 200bp upstream of TSS, or if not available over 1st exon probes.
#' 
#' 
#' @param beta
#' A DNAm beta-valued data matrix with columns labeling samples and rows labeling CpGs.
#' At present we only support Illumina 450k or 850k methylation beadarrays.
#' 
#' @param type
#' This specifies whether data is 450k or 850k.
#'
#' 
#' @return The DNAm data matrix at the gene-level.
#' 
#' 
#' @references 
#' Zhu T, Breeze CE, Beck S, Teschendorff AE.
#' \emph{Cell-type deconvolution of bulk tissue DNA methylomes 
#' from single-cell RNA-Seq data}
#' Submitted.
#'
#' Jiao Y, Widschwendter M, Teschendorff AE. 
#' \emph{A systems-level integrative framework for genome-wide DNA methylation
#' and gene expression data identifies differential gene expression modules
#' under epigenetic control.}
#' Bioinformatics (2014) 30(16).
#' doi:\href{https://doi.org/10.1093/bioinformatics/btu316}{
#' 10.1093/bioinformatics/btu316}.
#' 
#' @examples 
#' data(testDNAm);
#' avDNAm.m <- constAvBetaTSS(testDNAm.m,type="450k");
#' 
#'
#' 
#' @export
#'     

constAvBetaTSS <- function(beta.m,type=c("450k","850k")){
 if(type=="450k"){
  data("probeInfo450k");
 }
 else if (type=="850k"){
  data("probeInfo850k");
 }

map.idx <- match(rownames(beta.m), probeInfoALL.lv$probeID)
probeInfo.lv <- lapply(probeInfoALL.lv, function(tmp.v,ext.idx){return(tmp.v[ext.idx]);}, map.idx);
beta.lm <- list();
for (g in 1:6) {
    group.idx <- which(probeInfo.lv$GeneGroup == g)
    tmp.m <- beta.m[group.idx, ]
    rownames(tmp.m) <- probeInfo.lv$EID[group.idx];
    sel.idx <- which(is.na(rownames(tmp.m)) == FALSE);
    tmp.m <- tmp.m[sel.idx,];
    nL <- length(factor(rownames(tmp.m)));
    nspg.v <- summary(factor(rownames(tmp.m)),maxsum=nL);
    beta.lm[[g]] <- rowsum(tmp.m,group=rownames(tmp.m))/nspg.v;
    print(paste("Done for regional gene group ", g, sep = ""))
}
unqEID.v <- unique(c(rownames(beta.lm[[2]]), rownames(beta.lm[[4]])));
avbeta.m <- matrix(nrow = length(unqEID.v), ncol = ncol(beta.m))
colnames(avbeta.m) <- colnames(beta.m)
rownames(avbeta.m) <- unqEID.v
for (gr in c(4, 2)) {
      avbeta.m[match(rownames(beta.lm[[gr]]), rownames(avbeta.m)), ] <- beta.lm[[gr]]
}
return(avbeta.m);
}
