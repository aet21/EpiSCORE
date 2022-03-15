#' @title 
#' Construct mRNA expression reference matrix
#' 
#' @aliases ConstExpRef
#'  
#' @description 
#' This function takes as input a tissue-specific scRNA-Seq atlas and 
#' produces an expression reference matrix for the cell-types specified
#' by the user
#' 
#' 
#' @param exp.m
#' The normalized scRNA-Seq data matrix with rows labeling unique genes (either human
#' gene symbol or Entrez gene ID) and columns labeling cells. Missing values are not
#' allowed and values should be positive or zero. Importantly, the value corresponding
#' zero read counts should also be zero. Typically, this matrix will be obtained by
#' scaling the readcounts, adding a pseudocount of +1 and then taking the log2 transform.
#' 
#' @param celltype.idx
#' An integer index vector with as many elements as there are columns in 'exp.m',
#' specifying the cell-type of each cell. Each unique integer in the vector defines
#' one cell-type and integers must be sequential, i.e. no gaps are allowed. For instance,
#' for 4 cell-types, the integers would be 1 to 4.
#' 
#' @param namesCellT.v
#' A vector of unique names for the cell-types, of length equal to the number
#' of unique elements in celltype.idx. The order of the names must match the numerical
#' integers in celltype.idx.
#' 
#' @param markspecTH.v
#' This is the marker specificity score (MSS) threshold for each cell-type, i.e.
#' for a given marker gene and cell-type, the minimum number of other cell-types where
#' we require this gene not to be expressed. For n cell-types, the maximum MSS value is
#' n-1 and minimum value should be 1. The parameter must be a vector with one threshold
#' value for each cell-type. If 'NULL' (the default option), each entry will be set to n-1.
#' 
#' 
#' @return A list with two elements. The first element 'ref' is a list itself
#' also with 2 elements, giving the reference expression matrices defined
#' using averages ('av') or medians ('med'). The second entry 'markers' is
#' a list of matrices, one for each cell-type listing for each marker gene
#' in the reference, their Wilcoxon rank sum statistics, median value
#' in each cell-type and their marker specificity score.
#'
#' 
#' @references 
#' Teschendorff AE, Zhu T, Breeze CE, Beck S
#' \emph{Cell-type deconvolution of bulk tissue DNA methylomes 
#' from single-cell RNA-Seq data}
#' Genome Biol.2020
#'
#' Zhu T, Liu J, Beck S, Pan S, Capper D, Lechner M, Thirlwell C, Breeze CE, Teschendorff AE.
#' \emph{A pan-tissue DNA methylation atlas enables in silico decomposition of
#' human tissue methylomes at cell-type resolution.} Nat Methods 2022
#' 
#' @examples
#' data(lungSS2mca1)
#' out.l <- ConstExpRef(lungSS2mca1.m,celltypeSS2.idx,celltypeSS2.v,markspecTH.v=rep(3,4));
#' print(head(out.l$ref$med));
#' 
#' @importFrom presto wilcoxauc
#'
#' 
#' @export
#'     

ConstExpRef <- function(exp.m,celltype.idx,namesCellT.v,markspecTH.v=NULL){
nCT <- length(unique(celltype.idx));
if(nCT!=length(namesCellT.v)){
  print("PROBLEM WITH CELL-TYPE ANNOTATION");
}
medexpMCT.m <- matrix(NA,nrow=nrow(exp.m),ncol=nCT);
rownames(medexpMCT.m) <- rownames(exp.m);
colnames(medexpMCT.m) <- namesCellT.v;

sigWT.lm <- list();
statWT.lm <- list();
idx.l <- as.list(1:nrow(exp.m));

if(is.null(markspecTH.v)){
   markspecTH.v <- rep(nCT-1,nCT);
}

print("Finding marker genes");
for(ct in 1:nCT){
    ct.idx <- which(celltype.idx==ct);
    medexpMCT.m[,ct] <- apply(exp.m[,ct.idx],1,median);
    other.idx <- which(celltype.idx %in% setdiff(1:nCT,ct));
    pheno.v <- rep(NA,ncol(exp.m));
    pheno.v[ct.idx] <- 1;
    pheno.v[other.idx] <- 2;
    out.m <- wilcoxauc(exp.m,y=pheno.v,groups_use=1:2);
    nf <- nrow(out.m)/2;
    statWT.lm[[ct]] <- out.m[1:nf,c("logFC","auc","pval","padj")];
    colnames(statWT.lm[[ct]]) <- c("LFC","AUC","P","Padj(BH)");
    rownames(statWT.lm[[ct]]) <- rownames(exp.m);
    sig.idx <- which(statWT.lm[[ct]][,4] < 0.05);
    tmp.s <- sort(statWT.lm[[ct]][sig.idx,1],decreasing=TRUE,index.return=TRUE);
    sigWT.lm[[ct]] <- statWT.lm[[ct]][sig.idx[tmp.s$ix],];
}
names(sigWT.lm) <- namesCellT.v;
print(lapply(sigWT.lm,nrow));

print("Now compute marker specificity scores and filter markers");
markspecMCT.lv <- list();
markerMCT.lm <- list();
for(ct in 1:nCT){
    ### compute specificity score
    match(rownames(sigWT.lm[[ct]]),rownames(medexpMCT.m)) -> map.idx;
    tmp.v <- vector(); gi <- 1;
    for(g in map.idx){
       tmp.v[gi] <- length(which(medexpMCT.m[g,-ct]==0));
       gi <- gi+1;
    }
    names(tmp.v) <- rownames(sigWT.lm[[ct]]);
    markspecMCT.lv[[ct]] <- tmp.v;
    #### now filter
    match(rownames(sigWT.lm[[ct]]),rownames(medexpMCT.m)) -> map.idx;
    maxMCT.v <- unlist(apply(medexpMCT.m[map.idx,],1,which.max));
    sel.idx <- intersect(which(maxMCT.v==ct),intersect(which(medexpMCT.m[map.idx,ct]>0),which(markspecMCT.lv[[ct]]>=markspecTH.v[ct])));
    markerMCT.lm[[ct]] <- cbind(sigWT.lm[[ct]][sel.idx,],medexpMCT.m[map.idx[sel.idx],],markspecMCT.lv[[ct]][sel.idx]);
    rownames(markerMCT.lm[[ct]]) <- rownames(sigWT.lm[[ct]][sel.idx,]);
    colnames(markerMCT.lm[[ct]]) <- c(colnames(sigWT.lm[[ct]]),paste("Med-",colnames(medexpMCT.m),sep=""),"MARK.SPEC");
}
names(markspecMCT.lv) <- namesCellT.v;
names(markerMCT.lm) <- namesCellT.v;
print(lapply(markerMCT.lm,nrow));
    
print("Now construct reference");
common.v <- rownames(markerMCT.lm[[1]]);
for(ct in 2:nCT){
  common.v <- union(rownames(markerMCT.lm[[ct]]),common.v);
}

avexpMCT.m <- matrix(NA,nrow=length(common.v),ncol=nCT);
rownames(avexpMCT.m) <- common.v;
colnames(avexpMCT.m) <- namesCellT.v;
medexpMCT.m <- avexpMCT.m;
match(common.v,rownames(exp.m)) -> map.idx;
for(ct in 1:nCT){
    ct.idx <- which(celltype.idx==ct);
    avexpMCT.m[,ct] <- apply(exp.m[map.idx,ct.idx],1,mean);
    medexpMCT.m[,ct] <- apply(exp.m[map.idx,ct.idx],1,median);
}
refexpMCT.lm <- list(av=avexpMCT.m,med=medexpMCT.m);

return(list(ref=refexpMCT.lm,markers=markerMCT.lm));   

}

### Auxilliary function

### doWTprl
#doWTprl <- function(idx,data.m,group.li){
#    lfc <- mean(data.m[idx,group.li[[1]]]) - mean(data.m[idx,group.li[[2]]]);
#    wt.o <- wilcox.test(data.m[idx,group.li[[1]]],data.m[idx,group.li[[2]]]);
#    pv <- wt.o$p.value;
#    n.v <- unlist(lapply(group.li,length));
#    auc <- wt.o$stat/prod(n.v);
#    out.v <- c(lfc,auc,pv);
#    names(out.v) <- c("LFC","AUC","P");
#    return(out.v);
#}


