# Given a set of CNVs with copy number values, map them to genes and calculate the gene-level copy numbers
MapCNV2Gene<-function(cnv, exon, copy='copy', parent='parent', gene='gene') {
  # cnv         GRanges object of CNV locations
  # copy        If integer, must have the same length as <cnv> and correspond to copy number of each CNV
  #             If character, column name of <cnv> metadata, corresponding to copy number of each CNV
  # exon        GRanges object of location of exons, must have the same chromosome name and length as <cnv>. All exons must be named.
  # parent      The parent level of the exon, such as a transcript
  # gene        The gene where the exon/transcript is located in
  
  require('GenomicRanges'); 
  require('CHOPseq'); 
  
  if (length(copy)!=length(cnv)) copy<-elementMetadata(cnv)[, copy]; 
  copy[copy<0]<-0; 
  
  if (length(parent) != length(exon)) parent<-elementMetadata(exon)[, parent[1]];
  if (length(gene) != length(exon)) gene<-elementMetadata(exon)[, gene[1]];
  elementMetadata(exon)<-NULL;
  exon$parent<-parent;
  exon$gene<-gene;
  if (is.null(names(exon))) names(exon)<-1:length(exon); 
  
  cnv<-cnv[as.vector(seqnames(cnv)) %in% as.vector(seqnames(exon))]; 
  seqlengths(cnv)<-seqlengths(exon)[names(seqlengths(cnv))];
  
  # Fix potential overlapping of CNV regions
  cnv<-cnv[order(as.vector(seqnames(cnv)), start(cnv))]; 
  if (length(cnv)>1) {
    ind<-(2:length(cnv))[which(seqnames(cnv)[-1] == seqnames(cnv)[-length(cnv)])]; 
    start(cnv)[ind]<-pmax(start(cnv)[ind], end(cnv)[ind-1]+1); 
  } 
  
  # Copy number of the whole genome
  cov<-coverage(cnv, weight=copy+1); 
  cov[cov==0]<-3;
  cov<-cov-1;
  
  cat('Mapping', length(cnv), 'CNVs to', length(exon), 'exons\n'); 
  
  ########################################################################
  # Average exon copy number
  cat('Summarizing copy number at exon level\n'); 
  exon.copy<-rep(2, length(exon));
  names(exon.copy)<-names(exon); 
  
  ct<-countOverlaps(exon, cnv); # exon overlapped to CNVs?
  exon.olap<-exon[ct>0]; # exon overlapped to CNVs
  seqlevels(exon.olap)<-unique(as.vector(seqnames(exon.olap)));
  
  exon.cov<-cov[exon.olap]; 
  names(exon.cov)<-names(exon.olap); 
  exon.copy[names(exon.cov)]<-mean(exon.cov); 
  exon$copy<-exon.copy;

  ########################################################################
  # Break exons into intervals with the same copy number
  nrun<-sapply(exon.cov, nrun);
  ids<-rep(names(exon.cov), nrun);
  chr<-rep(as.vector(seqnames(exon.olap)), nrun); 
  str<-rep(as.vector(strand(exon.olap)), nrun); 
  stt<-rep(start(exon.olap), nrun); 

  stts<-unlist(lapply(exon.cov, start), use.names=FALSE); 
  ends<-unlist(lapply(exon.cov, end), use.names=FALSE); 
  cpys<-unlist(lapply(exon.cov, runValue), use.names=FALSE); 
  
  exon.intv<-GRanges(chr, IRanges(stt+stts-1, stt+ends-1), str, exon=ids, copy=cpys); 
  exon.dpld<-exon[!(names(exon) %in% ids)];
  elementMetadata(exon.dpld)<-NULL;
  exon.dpld$exon<-names(exon.dpld); 
  exon.dpld$copy<-2; 
  exon.intv<-c(exon.intv, exon.dpld);
  names(exon.intv)<-paste(exon.intv$exon, exon.intv$copy, sep='_'); 
  
  ########################################################################
  # summarize transcript/parent level
  cat('Summarizing copy number at transcript/exon set level\n'); 
  cpy<-split(as.vector(exon$copy), exon$parent); 
  wid<-split(width(exon), exon$parent); 
  str<-as.vector(sapply(split(as.vector(strand(exon)), exon$parent), function(x) unique(x)[1]));
  chr<-as.vector(sapply(split(as.vector(seqnames(exon)), exon$parent), function(x) unique(x)[1]));
  stt<-as.vector(sapply(split(start(exon), exon$parent), min));
  end<-as.vector(sapply(split(end(exon), exon$parent), max)); 
  len<-as.vector(sapply(wid, sum)); 
  cpy<-as.vector(sapply(names(cpy), function(nm) weighted.mean(cpy[[nm]], wid[[nm]]))); 
  tx<-GRanges(chr, IRanges(stt, end), str, exon=as.vector(sapply(wid, length)), length=len, copy=cpy); 
  names(tx)<-names(wid); 
  
  ########################################################################
  # summarize gene level
  cat('Summarizing copy number at gene level\n'); 
  mp<-gene;
  names(mp)<-parent; 
  mp<-mp[names(wid)]; 
  cpy<-split(as.vector(tx$copy), mp); 
  wid<-split(tx$length, mp); 
  str<-as.vector(sapply(split(as.vector(strand(tx)), mp), function(x) unique(x)[1]));
  chr<-as.vector(sapply(split(as.vector(seqnames(tx)), mp), function(x) unique(x)[1]));
  stt<-as.vector(sapply(split(start(tx), mp), min));
  end<-as.vector(sapply(split(end(tx), mp), max)); 
  len<-as.vector(sapply(split(tx$length, mp), mean)); 
  cpy<-as.vector(sapply(names(cpy), function(nm) weighted.mean(cpy[[nm]], wid[[nm]]))); 
  gn<-GRanges(chr, IRanges(stt, end), str, transcript=as.vector(sapply(wid, length)), length=len, copy=cpy); 
  names(gn)<-names(wid); 
  
  ########################################################################  
  # TSS upstream
  tss<-resize(tx, 1); 
  tss10k<-resize(tss, 10001, fix='end');
  tss.up<-resize(tss10k, 1, fix='center'); 
  cov.up<-apply(RetrieveSegmentList(cov, tss.up, 5000), 1, Rle); 
  nrun<-sapply(cov.up, nrun);
  ids<-rep(names(tss), nrun);
  chr<-rep(as.vector(seqnames(tss)), nrun); 
  str<-rep(as.vector(strand(tss)), nrun); 
  stt<-rep(start(tss10k), nrun); 
  
  stts<-unlist(lapply(cov.up, start), use.names=FALSE); 
  ends<-unlist(lapply(cov.up, end), use.names=FALSE); 
  cpys<-unlist(lapply(cov.up, runValue), use.names=FALSE); 
  
  up<-GRanges(chr, IRanges(stt+stts-1, stt+ends-1), str, copy=cpys, transcript=ids, tss=start(tss[ids])); 
  loc<-start(tss[ids]);
  dist1<--1*abs(start(up)-loc);
  dist2<--1*abs(end(up)-loc); 
  up$near<-pmax(dist1, dist2);
  up$distant<-pmin(dist1, dist2); 
  
  list(genome=cov, exon=exon, interval=exon.intv, transcript=tx, gene=gn, upstream=up); 
}