
# This function uses a 2-step procedure to normalize sequencing depth on a set of genomic locations betweeen samples.
# The main input of this function is a set of numeric matrixes, each corresponding to one sample or sequencing library.
# Each matrix should have the same numbers of rows and columns: 
# - rows are a set of genomic regions with the same or different width.
# - columns are sub-regions of each region. The sub-regions could have fixed or not fixed width, but must be ordered.
# - each cell in the matrix is the average sequencing depth of that sub-region.
#
# The 2-step normalization will be performed:
# Step 1. subtract background from the sequencing depth
# Step 2. normalize data between samples

NormChipseq1 <- function(cov, log.transform=TRUE, bg=0, bg.start=c(), bg.end=c(), 
                         trim.low=0.05, trim.high=0.95, exclude=c(), max.seed=25000) {
  require(DEGandMore);
  
  # Transform data to log-scale
  if (log.transform) cov <- lapply(cov, function(c) log2(c+1)); 

  # Subtract background if known a general background of all samples or background of individual samples
  if (length(bg)==1) cov <- lapply(cov, function(c) c-bg) else 
    if (length(bg)==length(cov)) cov <- lapply(1:length(cov), function(i) cov[[i]]-bg[i]); 
  
  # For trimmed mean
  mn <- min(max(0, trim.low), min(1, trim.high));
  mx <- max(max(0, trim.low), min(1, trim.high));
  
  if (length(bg.start)>0 & length(bg.start)==length(bg.end)) { # When use local region as background
    cov <- lapply(cov, function(c) {
      cv  <- c[!((1:nrow(c)) %in% exclude), , drop=FALSE];
      bg <- sapply(1:length(bg.start), function(i) {  # For each local region used for calculating background
        m <- rowMeans(cv[, bg.start[i]:bg.end[i], drop=FALSE], na.rm=TRUE);
        m <- m[!is.na(m)];
        m <- sort(m); 
        l <- length(m); 
        m <- m[max(1, round(mn*l)):min(l, round(mx*l))]; 
        mean(m); 
      });
      c - mean(bg); 
    });
  }; 
  
  nr <- nrow(cov[[1]]);
  nc <- ncol(cov[[1]]); 
  if (max.seed > nr) seed <- 1:nr else seed <- sort(sample(1:nr, max.seed)); 
  
  for (i in 1:nc) {
    d <- sapply(cov, function(c) c[, i]); 
    d <- NormLoess(d, seed = seed); 
    for (j in 1:length(cov)) cov[[j]][, i] <- d[, j]; 
  }; 
  
  if (log.transform) cov <- lapply(cov, function(c) exp(c*log(2))-1); 

  cov; 
}


