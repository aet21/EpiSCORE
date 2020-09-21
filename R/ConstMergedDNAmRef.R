#' @title 
#' Construct merged DNAm reference matrix
#' 
#' @aliases ConstMergedDNAmRef
#'  
#' @description 
#' This function takes as input two separate imputed DNAm reference matrices
#' and merges them to arrive at a unique final DNAm reference matrix. 
#' 
#' 
#' @param ref1.m
#' The first DNAm reference matrix, say generated from the SCM2 database
#' 
#' @param ref2.m
#' The second DNAm reference matrix, say generated from the RMAP database
#' 
#' @return The merged DNAm reference matrix
#' 
#' @references 
#' Teschendorff AE, Zhu T, Breeze CE, Beck S
#' \emph{Cell-type deconvolution of bulk tissue DNA methylomes 
#' from single-cell RNA-Seq data}
#' Genome Biol.2020
#' 
#' 
#' @examples 
#' data(lungSS2mca1)
#' out.l <- ConstExpRef(lungSS2mca1.m,celltypeSS2.idx,celltypeSS2.v,markspecTH.v=rep(3,4));
#' refDNAm1.m <- ImputeDNAmRef(out.l$ref$med,db="SCM2",geneID="SYMBOL");
#' refDNAm2.m <- ImputeDNAmRef(out.l$ref$med,db="RMAP",geneID="SYMBOL");
#' refDNAm.m <- ConstMergedDNAmRef(refDNAm1.m,refDNAm2.m);
#' print(head(refDNAm.m));
#'
#' 
#' 
#' @export
#'     
ConstMergedDNAmRef <- function(ref1.m,ref2.m){

common.v <- intersect(rownames(ref1.m),rownames(ref2.m));
union.v <- union(rownames(ref1.m),rownames(ref2.m));
unq.v <- setdiff(rownames(ref1.m),common.v);
unq2.v <- setdiff(rownames(ref2.m),common.v);

refMG.m <- matrix(NA,nrow=length(union.v),ncol=ncol(ref1.m));
colnames(refMG.m) <- colnames(ref1.m);
rownames(refMG.m) <- union.v;

refMG.m[match(unq.v,union.v),] <- ref1.m[match(unq.v,rownames(ref1.m)),];
refMG.m[match(unq2.v,union.v),] <- ref2.m[match(unq2.v,rownames(ref2.m)),];

tmp1.m <- ref1.m[match(common.v,rownames(ref1.m)),];
tmp2.m <- ref2.m[match(common.v,rownames(ref2.m)),];

nNA1.v <- apply(tmp1.m,1,function(v){return(length(which(is.na(v))));})
nNA2.v <- apply(tmp2.m,1,function(v){return(length(which(is.na(v))));})

refMG.m[match(common.v,union.v),] <-  0.5*(tmp1.m + tmp2.m); ## this deals with no NAs in both and with NAs in both

tmp2.idx <- intersect(which(nNA1.v>0),which(nNA2.v==0));
refMG.m[match(common.v[tmp2.idx],union.v),] <- tmp2.m[tmp2.idx,];

tmp1.idx <- intersect(which(nNA1.v==0),which(nNA2.v>0));
refMG.m[match(common.v[tmp1.idx],union.v),] <- tmp1.m[tmp1.idx,];

return(refMG.m);

} 



