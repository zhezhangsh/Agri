# Given a set of CNVs with copy number values, map them to a set of genomic features
# The genomic features could have 0 to 2 levels of parent features (ex. fea->transcript->gene)
# See the MapCNV2Gene function for more specifically mapping CNVs to genes
MapCNV2Feature<-function(cnv, fea, copy='copy', parent1=NA, parent2=NA) {
  # cnv         GRanges object of CNV locations
  # copy        Column name of <cnv> metadata corresponding to copy number of each CNV
  # fea         GRanges object of location of feas, must have the same chromosome name and length as <cnv>. All feas must be named.
  # parent1     The parent level of the basic feature units, such as transcript vs. feas
  # parent2     The parent level of <parent1>, such as gene vs. transcripts
  
  require('GenomicRanges'); 
  require('CHOPseq'); 
  
  if (is.null(names(cnv))) names(cnv) <- 1:length(cnv);
  if (is.null(names(fea))) names(fea) <- 1:length(fea); 
  
  copy <- copy[1]; 
  if (!(copy %in% names(elementMetadata(cnv)))) stop('Error: metadata column, ', copy, ' is not available.\n'); 
  cnv <- cnv[!is.na(as.numeric(elementMetadata(cnv)[, copy]))]; 
  
  chr <- unique(as.vector(seqnames(fea))); 
  seqlevels(fea) <- chr;
  len <- sapply(split(end(fea), as.vector(seqnames(fea))), max);
  seqlengths(fea) <- base::pmax(seqlengths(fea), len[seqlevels(fea)], na.rm=TRUE); 
  
  # Match chromosome name and length
  cnv <- cnv[as.vector(seqnames(cnv)) %in% seqlevels(fea)]; 
  if (length(cnv) == 0) warning('None of the CNVs was overlapped to any genomic features; make sure their chromosome names match.\n'); 
  seqlevels(cnv) <- seqlevels(fea); 
  len <- sapply(split(end(cnv), as.vector(seqnames(cnv))), max)[seqlevels(cnv)];
  seqlengths(fea)[seqlevels(cnv)] <- seqlengths(cnv) <- pmax(seqlengths(cnv)[seqlevels(cnv)], seqlengths(fea)[seqlevels(cnv)], len, na.rm=TRUE)

  # Fix potential overlapping of CNV regions
  cnv<-cnv[order(as.vector(seqnames(cnv)), start(cnv))]; 
  if (length(cnv)>1) {
    ind<-(2:length(cnv))[which(as.vector(seqnames(cnv)[-1]) == as.vector(seqnames(cnv)[-length(cnv)]))]; 
    start(cnv)[ind]<-pmax(start(cnv)[ind], end(cnv)[ind-1]+1); 
  } 
  
  nm <- names(cnv); 
  if (is.null(nm)) names(cnv) <- 1:length(cnv) else names(cnv)[which(is.na(nm))] <- which(is.na(nm)); 
  nm <- names(fea); 
  if (is.null(nm)) names(fea) <- 1:length(fea) else names(fea)[which(is.na(nm))] <- which(is.na(nm)); 
  
  copy <- as.numeric(elementMetadata(cnv)[, copy]); 
  
  # Copy number of the whole genome
  cov <- coverage(cnv, weight=copy+1); 
  cov[cov==0] <- 3;
  cov <- cov-1;
  
  cat('Mapping', length(cnv), 'CNVs to', length(fea), 'feas\n'); 
  
  ########################################################################
  # Average copy number of each genomic region
  cat('Summarizing copy number at each given genomic region\n'); 
  fea.copy<-rep(2, length(fea));
  names(fea.copy)<-names(fea); 
  
  ct<-countOverlaps(fea, cnv); # feature overlapped to CNVs?
  fea.olap<-fea[ct>0]; # feature overlapped to CNVs
  seqlevels(fea.olap)<-unique(as.vector(seqnames(fea.olap)));
  
  c<-data.frame(chr=as.vector(seqnames(fea.olap)), start=start(fea.olap), end=end(fea.olap), stringsAsFactors = FALSE); 
  fea.cov<-lapply(1:nrow(c), function(i) cov[[c[i, 1]]][c[i, 2]:c[i, 3]]); 
  names(fea.cov)<-names(fea.olap); 
  fea.copy[names(fea.olap)]<-sapply(fea.cov, mean); 
  fea$copy<-fea.copy;
  
  olap <- findOverlaps(fea.olap, cnv);
  nm1 <- names(fea.olap)[olap@queryHits]; 
  nm2 <- names(cnv)[olap@subjectHits];

  out<-list(cnv=cnv, map2cnv=split(nm2, nm1), feature=fea); 

  ########################################################################
  
  parents <- c(parent1[1], parent2[1]); 
  parents <- parents[!is.na(parents)]; 
  if (length(setdiff(parents, names(elementMetadata(fea)))) > 0) {
    x <- setdiff(parents, names(elementMetadata(fea))); 
    warning('Warning: parent level name(s) ', paste(x, collapse=' and '), ' not column of feature metadata\n'); 
    parents <- parents[!(parents %in% x)]; 
  }
  
  ########################################################################
  if (length(parents) > 0) {
  
    ########################################################################
    # summarize parent level 1
    cat('Summarizing the first parent level:', parents[1], '\n'); 
    id0 <- as.vector(elementMetadata(fea)[, parents[1]]); 
    id1 <- unique(id0[fea$copy!=2]); 
    fea1 <- fea[as.vector(elementMetadata(fea)[, parents[1]]) %in% id1]; 
    cpy1 <- sapply(id1, function(i) {
      x <- fea1[as.vector(elementMetadata(fea1)[, parents[1]]) %in% i]; 
      weighted.mean(x$copy, w = width(x)); 
    }); 
    id <- unique(id0); 
    cpy <- rep(2, length(id));
    names(cpy) <- id; 
    cpy[id1] <- cpy1;
    
    map2cnv <- out$map2cnv; 
    names(map2cnv) <- elementMetadata(fea[names(out$map2cnv)])[, parents[1]]; 
    x <- unlist(map2cnv, use.names=FALSE);
    y <- rep(names(map2cnv), sapply(map2cnv, length)); 
    map2cnv <- lapply(split(x, y), unique); 
    
    out[[parents[1]]]<-list(copy=cpy, length=sapply(split(width(fea), id0), sum), mapping=split(names(fea), id0)[names(cpy)], map2cnv=map2cnv); 

    if (length(parents) > 1) {
      ########################################################################
      # summarize parent level 2
      cat('Summarizing the seconde parent level:', parents[2], '\n');  
      mp <- split(id0, as.vector(elementMetadata(fea)[, parents[2]])); 
      mp <- lapply(mp, unique); 
      p1 <- unlist(mp, use.names=FALSE);
      p2 <- rep(names(mp), sapply(mp, length)); 
      l <- out[[length(out)]]$length[p1];
      c <- out[[length(out)]]$copy[p1]; 
      lc <- cbind(split(as.vector(l), as.vector(p2)), split(as.vector(c), as.vector(p2))); 
      cp <- apply(lc, 1, function(x) weighted.mean(x[[2]], w=x[[1]]));
      
      names(p2) <- p1; 
      names(map2cnv) <- p2[names(map2cnv)]; 
      map2cnv <- lapply(split(unlist(map2cnv, use.names=FALSE), rep(names(map2cnv), sapply(map2cnv, length))), unique); 
        
      out[[parents[2]]] <- list(copy=cp, mapping=mp, map2cnv=map2cnv); 
    }
  }
  
  out; 
}


