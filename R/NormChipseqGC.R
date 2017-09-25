# This function adjusts the difference of ChIP-seq depth between library pairs with local GC content
# Given a matrix of ChIP-seq depth, the columns are ChIP-seq libraries and the rows are genomic regions of interest, such as TSS and enhancers.
# The function calculates the difference of all pairs of columns and adjusts it with GC contents within the regions.
# The goal of this function is to remove the impact of GC contents on differential ChIP-seq analysis between sample groups.

NormChipseqGC <- function(d, gc) {
  rmn <- rowMeans(d); 
  prs <- t(combn(1:ncol(d), 2)); 
  dff <- apply(prs, 1, function(ind) {
    cat(ind[1], '-', ind[2], '\n');
    dff <- d[, ind[2]]-d[, ind[1]]; 
    loess(dff ~ gc)$residuals;
  }); 
  adj <- sapply(1:ncol(d), function(i) {
    ind <- which(prs[,1]==i | prs[, 2]==i); 
    sub <- dff[, ind]; 
    sgn <- sapply(ind, function(j) if (prs[j, 1]==i) 1 else -1);
    rmn-colMeans(sgn*t(sub)); 
  }); 
  colnames(adj) <- colnames(d); 
  adj; 
}; 

