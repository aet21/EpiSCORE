#' @title 
#' Impute a DNAm reference matrix
#' 
#' @aliases ImputeDNAmRef
#'  
#' @description 
#' This function takes as input a mRNA expression reference matrix and will generate
#' an imputed DNAm reference using one of two databases (NIH Epigenomics Roadmap
#' or Stem-Cell Matrix Compendium) specified by the user.
#' 
#' 
#' @param refexp.m
#' The gene expression reference matrix with rows labeling genes and columns
#' labeling cell-types. If gene idgentifier is gene symbol this will be converted
#' to Entrez gene ID.
#' 
#' @param db
#' A character specifying the database of matched expression and DNAm data.
#' This has to be either SCM2 for Stem-Cell-Matrix Compendium-2
#' or RMAP for Epigenomics Roadmap.
#'
#' @param geneID
#' A character specifying gene identifier used in expression reference matrix input. 
#' Either SYMBOL or ENTREZ GENE ID.
#'
#' 
#' @return The imputed DNAm reference matrix defined over the same cell-types
#' as defined in the input expression reference matrix.
#' 
#' 
#' @references 
#' Teschendorff AE, Zhu T, Breeze CE, Beck S
#' \emph{Cell-type deconvolution of bulk tissue DNA methylomes 
#' from single-cell RNA-Seq data}
#' Genome Biol.2020
#' 
#' Zhu T, Liu J, Beck S, Pan S, Capper D, Lechner M, Thirlwell C, Breeze CE, Teschendorff AE.
#' \emph{A pan-tissue DNA methylation atlas enables in silico decomposition
#' of human tissue methylomes at cell-type resolution.} Nat Methods 2022
#' 
#' @examples 
#' data(lungSS2mca1)
#' out.l <- ConstExpRef(lungSS2mca1.m,celltypeSS2.idx,celltypeSS2.v,markspecTH.v=rep(3,4));
#' refDNAm.m <- ImputeDNAmRef(out.l$ref$med,db="SCM2",geneID="SYMBOL");
#' print(head(refDNAm.m));
#'
#' 
#' 
#' @export
#'     

ImputeDNAmRef <- function(refexp.m,db=c("SCM2","RMAP"),geneID=c("SYMBOL","ENTREZID")){

if(geneID=="SYMBOL"){
  eid.v <- convertIDs(rownames(refexp.m), "SYMBOL", "ENTREZID", org.Hs.eg.db,ifMultiple="useFirst");
  na.idx <- which(is.na(eid.v));
  xx <- as.list(org.Hs.egALIAS2EG)
  xx <- xx[!is.na(xx)];
  map.idx <- match(rownames(refexp.m)[na.idx],names(xx));
  eidNA.v <- vector();
  for(n in 1:length(na.idx)){
    if(!is.na(map.idx[n])){
        eidNA.v[n] <- xx[[map.idx[n]]][1];
    }
  }
  eid.v[na.idx] <- eidNA.v;
  refexpEID.v <- eid.v;
}
else if (geneID=="ENTREZID"){
   refexpEID.v <- rownames(refexp.m);
}

    if(db=="SCM2"){
        data("dbSCM2");
        topIntMR.m <- topIntMRscm2.m;
        intDNAm.m <- intSCM2dnaM.m;
        intExpR.m <- intSCM2expR.m;
        pEgX.m <- pEgXscm2.m;
    }
    else if (db=="RMAP"){
        data("dbRMAP");
        topIntMR.m <- topIntMRrmap.m;
        intDNAm.m <- betaRMAPint.m;
        intExpR.m <- expRMAPint.m;
        pEgX.m <- pEgXrmap.m;
    }

markersEID.v <- intersect(refexpEID.v,rownames(topIntMR.m));
map.idx <- match(markersEID.v,rownames(intDNAm.m));

tmpRefM.m <- matrix(NA,nrow=length(markersEID.v),ncol=ncol(refexp.m));
colnames(tmpRefM.m) <- colnames(refexp.m);
rownames(tmpRefM.m) <- markersEID.v;

mapRef.idx <- match(markersEID.v,refexpEID.v);
mapExp.idx <- match(markersEID.v,rownames(intExpR.m));

for(ct in 1:ncol(refexp.m)){
    mark.idx <- which(refexp.m[mapRef.idx,ct]>0);
    notE.idx <- which(refexp.m[mapRef.idx,ct]==0);
    tmpRefM.m[mark.idx,ct] <- 0; ### impute zero DNAm
    p.m <- pEgX.m[mapExp.idx[notE.idx],];
    m.m <- intDNAm.m[mapExp.idx[notE.idx],];
    for(g in 1:nrow(p.m)){
     notExpS.idx <- which(p.m[g,] < 0.2);
     if(length(notExpS.idx)>0){
         tmpRefM.m[notE.idx[g],ct] <- median(m.m[g,notExpS.idx]);
     }
    }
}

wMref.v <- vector();
for(g in 1:nrow(tmpRefM.m)){
  ct <- which(tmpRefM.m[g,]==0);
  wMref.v[g] <- mean(tmpRefM.m[g,-ct]);
}
wMref.v[which(is.na(wMref.v))] <- 0;

refMscbu.m <- cbind(tmpRefM.m,wMref.v);
colnames(refMscbu.m) <- c(colnames(tmpRefM.m),"weight");

return(refMscbu.m);
} 

### Auxilliary function
#' @import org.Hs.eg.db
#' @importFrom AnnotationDbi select
### convertIDs
convertIDs <- function( ids, from, to, db, ifMultiple=c("putNA", "useFirst")) {
  stopifnot( inherits( db, "AnnotationDb" ) )
  ifMultiple <- match.arg( ifMultiple )
  suppressWarnings( selRes <- AnnotationDbi::select(
    db, keys=ids, keytype=from, columns=c(from,to) ) )
 
  if ( ifMultiple == "putNA" ) {
    duplicatedIds <- selRes[ duplicated( selRes[,1] ), 1 ]
    selRes <- selRes[ ! selRes[,1] %in% duplicatedIds, ]
  }
 
  return( selRes[ match( ids, selRes[,1] ), 2 ] )
}



