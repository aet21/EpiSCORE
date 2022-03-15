#' @title 
#' Estimation of cell-type fractions using a DNAm reference matrix
#' 
#' @aliases wRPC
#'  
#' @description 
#' This function takes as input a number of genome-wide DNAm profiles from
#' purified or complex bulk tissue samples and estimates the proportions
#' of the main cell-types in the samples.
#' 
#' 
#' @param data
#' A DNAm data matrix with columns labeling samples and rows labeling genes (using
#' Entrez gene identifier). This assumes that DNAm has been summarized at the gene
#' level, for instance by averaging DNAm values over CpGs or probes that map to
#' 200bp upstream of the gene's TSS. data can also be a vector if there is only one
#' sample.
#' 
#' @param ref.m
#' The DNAm reference matrix to be used with rows labeling genes (using Entrez gene
#' identifier) and columns labeling the main cell-types, with last column labeling
#' the weight. Weights need to be between 0 and 1. 
#'
#' @param useW
#' A logical. If 'TRUE' weights are used in the regression, otherwise not.
#'
#' @param wth
#' A threshold on the weights to select most informative genes. Only used
#' if \code{useW} is \code{TRUE}.
#' 
#' @param maxit
#' The maximum number of iterations in the robust linear regression.
#' 
#' @return A list containing the following elements
#'
#' @return estF
#' A matrix of estimated cell-type fractions, with rows labeling samples and
#' columns labeling the cell-types in the reference matrix.
#'
#' @return ref
#' The reference matrix used, i.e for the genes overlapping with those in
#' the data matrix.
#' 
#' 
#' @references 
#' Teschendorff AE, Zhu T, Breeze CE, Beck S.
#' \emph{Cell-type deconvolution of bulk tissue DNA methylomes 
#' from single-cell RNA-Seq data}
#' Genome Biol.2020
#'
#' Zhu T, Liu J, Beck S, Pan S, Capper D, Lechner M, Thirlwell C, Breeze CE, Teschendorff AE.
#' \emph{A pan-tissue DNA methylation atlas enables in silico decomposition of
#' human tissue methylomes at cell-type resolution.} Nat Methods 2022
#'
#' Teschendorff AE, Breeze CE, Zheng SC, Beck S. 
#' \emph{A comparison of reference-based algorithms for correcting cell-type 
#' heterogeneity in Epigenome-Wide Association Studies.}
#' BMC Bioinformatics (2017) 18: 105.
#' doi:\href{https://doi.org/10.1186/s12859-017-1511-5}{
#' 10.1186/s12859-017-1511-5}.
#' 
#' @examples 
#' data(lungSS2mca1)
#' out.l <- ConstExpRef(lungSS2mca1.m,celltypeSS2.idx,celltypeSS2.v,markspecTH.v=rep(3,4));
#' refDNAm1.m <- ImputeDNAmRef(out.l$ref$med,db="SCM2",geneID="SYMBOL");
#' refDNAm2.m <- ImputeDNAmRef(out.l$ref$med,db="RMAP",geneID="SYMBOL");
#' refDNAm.m <- ConstMergedDNAmRef(refDNAm1.m,refDNAm2.m);
#' data(testDNAm);
#' avDNAm.m <- constAvBetaTSS(testDNAm.m,type="450k");
#' wRPC.o <- wRPC(avDNAm.m,refDNAm.m,useW=TRUE,wth=0.4,maxit=200);
#' print(head(wRPC.o$est));
#'
#'
#' @importFrom MASS rlm
#' 
#' @export
#'     

wRPC <- function(data,ref.m,useW=TRUE,wth=0.4,maxit=100){
    if(is.matrix(data)){
     common.v <- intersect(rownames(ref.m),rownames(data));
     map.idx <- match(common.v,rownames(data));
     rep.idx <- match(common.v,rownames(ref.m));     
     data2.m <- data[map.idx,];
     ref2.m <- ref.m[rep.idx,-ncol(ref.m)];
     wG2.v <- ref.m[rep.idx,ncol(ref.m)];
    }
    else if (is.vector(data)){
      common.v <- intersect(rownames(ref.m),names(data));
      map.idx <- match(common.v,names(data));
      rep.idx <- match(common.v,rownames(ref.m));     
      data2.m <- matrix(data[map.idx],nrow=length(map.idx),ncol=1);
      rownames(data2.m) <- names(data[map.idx]);
      colnames(data2.m) <- c("sample");
      ref2.m <- ref.m[rep.idx,-ncol(ref.m)];
      wG2.v <- ref.m[rep.idx,ncol(ref.m)];
    }
    if(useW){
      selG.idx <- which(wG2.v>wth);
      if(length(selG.idx)<10){
          stop("Insufficient number of informative genes, try relaxing weight threshold!");
      }
    }
    est.m <- matrix(nrow=ncol(data2.m),ncol=ncol(ref2.m));
    colnames(est.m) <- colnames(ref2.m);
    rownames(est.m) <- colnames(data2.m);
    if(useW){
      for(s in 1:ncol(data2.m)){
          y <- sqrt(wG2.v[selG.idx])*data2.m[selG.idx,s];
          X <- diag(sqrt(wG2.v[selG.idx])) %*% ref2.m[selG.idx,];
          lm.o <- rlm( y ~ X , maxit=maxit);
          coef.v <- summary(lm.o)$coef[2:(ncol(ref2.m)+1),1];
          coef.v[which(coef.v<0)] <- 0;
          total <- sum(coef.v);
          coef.v <- coef.v/total;
          est.m[s,] <- coef.v;
      }
    }
    else {
      for(s in 1:ncol(data2.m)){
          rlm.o <- rlm( data2.m[,s] ~ ref2.m ,maxit=maxit);
          coef.v <- summary(rlm.o)$coef[2:(ncol(ref2.m)+1),1];
          coef.v[which(coef.v<0)] <- 0;
          total <- sum(coef.v);
          coef.v <- coef.v/total;
          est.m[s,] <- coef.v;
      }
    }
      
    
    return(list(estF=est.m,ref=ref2.m));
}



