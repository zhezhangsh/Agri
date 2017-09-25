# Remove small gaps from GRanges object
# Strand information and elementMetadata will be ignored
RemoveGap <- function(gr, size=0) {
  # gr    GRanges object
  # size  Maximum size of gaps to be removed
  
  require(GenomicRanges);
  
  if (size < 0) size <- 0;
  
  gr <- GRanges(as.vector(seqnames(gr)), IRanges(start(gr), end(gr)));
  cov <- coverage(gr); 
  
  cov <- lapply(cov, function(c) {
    v <- c@values;
    l <- c@lengths;
    v[v==0 & l<=size] <- 1;
    Rle(v, l); 
  });
  
  rng <- lapply(cov, function(c) {
    rng <- cbind(start(c), end(c), runValue(c));
    rng[rng[, 3]>0, , drop=FALSE];
  });
  
  chr <- rep(names(cov), sapply(rng, nrow)); 
  stt <- unlist(lapply(rng, function(x) x[, 1]));
  end <- unlist(lapply(rng, function(x) x[, 2]));
  
  gr <- GRanges(as.vector(chr), IRanges(as.vector(stt), as.vector(end)));
  gr; 
}