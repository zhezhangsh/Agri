ScorePWMWrapper <- function(pwm.set, subject.set, num.cluster=4, both.strand=TRUE, local.max=TRUE, sampling=10000, pvalue=0.01, fdr=FALSE) {
  # pwm          A list of PWM objects to be passed to the ScorePWMMatch() function: "For PWM: a rectangular character vector or rectangular DNAStringSet object ("rectangular" means that all elements have the same number of characters) with no IUPAC ambiguity letters, or a Position Frequency Matrix represented as an integer matrix with row names containing at least A, C, G and T (typically the result of a call to consensusMatrix)."
  # subject      A DNAStringSet or character vector of DNA sequences
  # num.cluster  Number of clusters for parallele computing
  # both.strand  If TRUE, the reverse-complement sequence of <subject> will also be matched
  # local.max    If TRUE, matches overlapped to any matches with higher scores will be removed.
  # sampling     Number of random sampling to get background distribution of matching scores
  # pvalue       Cutoff of empirical p value based on matching scores of random sampling
  # fdr          If TRUE, calculate FDR of the p values based on total length of the sequences

  require("Biostrings");
  require("parallel");

  # PWM matrix
  pwm.set <- lapply(pwm.set, function(pwm) if (is.character(pwm) | is.factor(pwm) | class(pwm)=='DNAStringSet') ConvertSeq2PWM(pwm) else pwm);

  # Union sequences
  if (is.character(subject.set)) subject.set <- DNAStringSet(subject.set);
  seq.all <- BiocGenerics::unlist(subject.set);

  # Make sure all sequences have names
  nms <- names(subject.set);
  if (is.null(nms)) nms <- paste('Seq', 1:length(subject.set), sep='_');
  ind <- which(is.na(nms) | nms=='');
  if (length(ind) > 0) nms[ind] <- paste('Seq', ind, sep='_');
  names(subject.set) <- nms;

  # Sequence break points
  len <- width(subject.set);
  end <- cumsum(len);
  stt <- c(1, 1+end[-length(end)]);
  rng <- Views(seq.all, start = stt, end = end);

  # Minimal matching scores
  frq <- colSums(alphabetFrequency(subject.set));
  frq <- frq[1:4]/sum(frq);
  rsc <- mclapply(pwm.set, function(pwm) {
    s <- ScorePWMDist(pwm, frq=frq, num=sampling);
    c <- s[min(ceiling(pvalue*length(s)), length(s))];
    c <- floor(100*(c/maxScore(pwm)));
    list(c, s);
  }, mc.cores = num.cluster);
  min.score <- as.vector(sapply(rsc, function(x) x[[1]]));
  rsc <- lapply(rsc, function(x) x[[2]]);

  ########################################################################################
  match.all <- mclapply(1:length(pwm.set), function(i) {
    ScorePWMMatch(pwm.set[[i]], seq.all, both.strand = both.strand, local.max = local.max,
                  min.score = min.score[i], min.match = 1);
  }, mc.cores = num.cluster);
  ########################################################################################

  ################################################################
  # Split matches
  tbls <- lapply(1:length(match.all), function(i) {
    m <- match.all[[i]];
    s <- rsc[[i]];

    olap  <- findOverlaps(m, rng, type='within');
    # qry   <- olap@queryHits;
    # sub   <- olap@subjectHits;
    qry   <- olap@from;
    sub   <- olap@to;
    stts  <- start(rng)[sub];
    mcol  <- mcols(m)[qry, ];

    seqID <- names(subject.set)[sub];
    stt   <- start(m)[qry] - stts + 1;
    end   <- end(m)[qry] - stts + 1;
    str   <- mcols(m)$strand[qry];
    score <- mcols(m)$percent[qry];
    pval  <- ScorePWMPvalue(mcols(m)$percent[qry]/100*maxScore(pwm.set[[i]]), s);
    seqs  <- as.character(m)[qry];
    seqs[str==-1] <- as.character(reverseComplement(DNAStringSet(seqs[str==-1])));

    t <- data.frame(seqID=seqID, start=stt, end=end, strand=str, score=score, pvalue=pval, seq=seqs, stringsAsFactors = FALSE);
    t <- t[rev(order(t$score)), , drop=FALSE];
    rownames(t) <- 1:nrow(t);

    # Calculate FDR
    if (fdr) {
      n  <- sum(len) - ncol(pwm.set[[i]]) + length(len);
      p0 <- seq(1/n/2, 1-1/n/2, 1/n);
      p0 <- t$pvalue;
      p1 <- c(p0, rep(1, n-nrow(t)));
      t$fdr  <- round(p.adjust(p1, method='BH')[1:nrow(t)], 4);
      t <- t[, c(1, 2, 3, 4, 5, 6, 8, 7)];
    }

    t;
  });
  names(tbls) <- names(pwm.set);

  tbls;
}

# A more efficient way of PWM to DNA sequence match of 1 or more PWMs and 1 or more sequences
ScorePWMMatchSet <- function(pwm.set, subject.set, num.cluster=4, both.strand=TRUE, local.max=TRUE, min.score=85, min.match=0) {
  # pwm          A list of PWM objects to be passed to the ScorePWMMatch() function: "For PWM: a rectangular character vector or rectangular DNAStringSet object ("rectangular" means that all elements have the same number of characters) with no IUPAC ambiguity letters, or a Position Frequency Matrix represented as an integer matrix with row names containing at least A, C, G and T (typically the result of a call to consensusMatrix)."
  # subject      A DNAStringSet or character vector of DNA sequences
  # num.cluster  Number of clusters for parallele computing
  # both.strand  If TRUE, the reverse-complement sequence of <subject> will also be matched
  # local.max    If TRUE, matches overlapped to any matches with higher scores will be removed.
  # min.score    Minimal match score, between 0 and 100
  # min.match    Minimal number of matches within all sequences to return; reduce matching score if getting less matches using <min.score>

  require("Biostrings");
  require("parallel");

  # Union sequences
  if (is.character(subject.set)) subject.set <- DNAStringSet(subject.set);
  seq.all <- BiocGenerics::unlist(subject.set);

  # Sequence break points
  len <- width(subject.set);
  end <- cumsum(len);
  stt <- c(1, 1+end[-length(end)]);
  rng <- Views(seq.all, start = stt, end = end);

  # parameters
  if (length(min.score) == 1) min.score <- rep(min.score, length(pwm.set));
  if (length(min.match) == 1) min.match <- rep(min.match, length(pwm.set));
  num.cluster <- min(length(pwm.set), num.cluster);

  ########################################################################################
  match.all <- mclapply(1:length(pwm.set), function(i) {
    ScorePWMMatch(pwm.set[[i]], seq.all, both.strand = both.strand, local.max = local.max,
                  min.score = min.score[i], min.match = min.match[i]);
  }, mc.cores = num.cluster);
  ########################################################################################

  ################################################################
  # Split matches
  match.all <- lapply(match.all, function(m) {
    olap <- findOverlaps(m, rng, type='within');
    mtch <- lapply(1:length(subject.set), function(i) {
      # ind <- olap@queryHits[olap@subjectHits==i];
      ind <- olap@from[olap@to==i];
      if (length(ind) == 0) Views(subject.set[[i]]) else {
        x <- m[ind];
        v <- Views(subject.set[[i]], start = start(x) - stt[i] + 1, end = end(x) - stt[i] + 1);
        v@elementMetadata <- x@elementMetadata;
        v;
      }
    });
    names(mtch) <- names(subject.set);
    mtch;
  });
  names(match.all) <- names(pwm.set);
  ################################################################

  match.all;
}

# PWM to DNA sequence match with matching scores.
ScorePWMMatch <- function(pwm, subject, both.strand=TRUE, local.max=TRUE, min.score=85, min.match=0) {
  # pwm          Same as the first parameter of the matchPWM() function: "For PWM: a rectangular character vector or rectangular DNAStringSet object ("rectangular" means that all elements have the same number of characters) with no IUPAC ambiguity letters, or a Position Frequency Matrix represented as an integer matrix with row names containing at least A, C, G and T (typically the result of a call to consensusMatrix)."
  # subject      Same as the second parameter of the matchPWM() function: "Typically a DNAString object. A Views object on a DNAString subject, a MaskedDNAString object, or a single character string, are also supported. IUPAC ambiguity letters in subject are ignored (i.e. assigned weight 0) with a warning."
  # both.strand  If TRUE, the reverse-complement sequence of <subject> will also be matched
  # local.max    If TRUE, matches overlapped to any matches with higher scores will be removed.
  # min.score    Minimal match score, between 0 and 100
  # min.match    Minimal number of matches to return; reduce matching score if getting less matches using <min.score>

  require("Biostrings");

  # Convert PWM from original format
  if (class(pwm) == "DNAStringSet" | is.character(pwm) | is.factor(pwm)) { # PWM is provided as DNAStringSet of seed sequences
    pwm <- ConvertSeq2PWM(pwm);
  } else if (is.matrix(pwm) | is.data.frame(pwm) | is.table(pwm)) { # PWM is provided as integer matrix
    pwm <- as.matrix(pwm)[1:4, , drop=FALSE];
    mx  <- max(pwm, na.rm=TRUE);
    if (mx <= 1) { # PWM is given as proportion, guess the total number of seeds and transform the matrix
      mn  <- min(pwm[pwm>0 & !is.na(pwm)]);  # non-zero minimum
      adj <- 1/mn;
      pwm <- round(pwm*adj)
    }
  } else {
    stop('Unknown format of PWM: ', class(pwm), '\n');
  }
  rownames(pwm) <- c('A', 'C', 'G', 'T');
  pwm[is.na(pwm)] <- 0;
  if (!both.strand) pwm <- list(forward=pwm) else {
    rev <- pwm[4:1, ncol(pwm):1];
    rownames(rev) <- rownames(pwm) <- c('A', 'C', 'G', 'T');
    pwm <- list(forward=pwm, reverse=rev);
  }

  # Max possible score of match
  mx <- maxScore(pwm[[1]]);
  sc <- paste(min(100, max(0, min.score)), '%', sep='');

  # Reduce score cutoff if less matches than required
  count <- sum(sapply(pwm, function(pwm) countPWM(pwm, subject, min.score = sc)));
  if (count < min.match) sc <- '0%'; # This is a very aggresive step to match all bases, instead, try to guess a better score cutoff,

  ###################################################################
  # Match PWM
  match <- lapply(pwm, function(pwm) matchPWM(pwm, subject, min.score = sc, with.score = TRUE));
  ###################################################################

  # Organize matches
  str   <- rep(c(1, -1)[1:length(match)], sapply(match, length));
  if (length(match) == 1) match <- match[[1]] else match <- c(match[[1]], match[[2]]);
  mcols(match)$strand = str;
  match <- match[order(start(match))];

  # Make sure minimum number of matches, by lower cutoff score
  if (min.match > 0) {
    s <- 100*mcols(match)$score/mx;
    if (min(s) < min.score) {
      o <- rev(sort(s));
      c <- o[min(length(o), min.match)];
      c <- min(c, min.score);
      match <- match[s>=c];
    };
  }

  # Remove match if it overlaps to another match with higher score
  if (local.max & length(match)>0) {
    olap <- findOverlaps(match, match);
    # olap <- cbind(olap@queryHits, olap@subjectHits);
    olap <- cbind(olap@from, olap@to);
    olap <- olap[olap[, 1]!=olap[, 2], , drop=FALSE];
    if (nrow(olap) > 0) {
      s0 <- mcols(match)$score;
      s1 <- s0[olap[, 1]];
      s2 <- s0[olap[, 2]];
      ls <- olap[s1 < s2, , drop=FALSE];
      if (nrow(ls) > 0) {
        ind <- unique(ls[, 1]);
        match <- match[-ind];
        score <- mcols(match)$score;
      }
    }
  }

  # Percent of maximum score
  mcols(match)$percent <- round(100*mcols(match)$score/mx, 1);
  mcols(match) <- mcols(match)[c(2, 1, 3)];

  match;
}

# Calculate distribution of PWM match scores by random sampling
ScorePWMDist <- function(pwm, frq=NA, seq=NA, num=4096) {
  # pwm   A Position Weight Matrix represented as a numeric matrix with row names A, C, G and T.
  # frq   Known base frequency. Will be overwritten if background sequence is provided by <seq>
  # seq   A character vector or DNAString or DNAStringSet, used to calculate background base frequency. Use 0.25 for all 4 bases if not given.
  # num   Number of random sampling to draw a random sequence, recommended value: 1000 - 100,000

  # Base frequency
  if(identical(NA, frq) | length(frq)!=4) frq <- rep(0.25, 4);
  if (!identical(seq, NA)) {
    if (is.character(seq) | is.factor(seq)) seq <- DNAStringSet(seq);
    if (class(seq) == 'DNAString') seq <- DNAStringSet(seq);
    if (class(seq) == 'DNAStringSet') {
      cnt <- colSums(alphabetFrequency(seq)[, 1:4]);
      frq <- cnt/sum(cnt);
    }
  };

  # Random scores
  s <- sapply(1:num, function(i) sum( apply(pwm, 2, function(x) sample(x, 1, prob=frq)) ));

  rev(sort(s));
};

# Calculate the empiracal p values of matching score, based on random background scores from function <ScorePWMDist>
ScorePWMPvalue <- function(score, bg) {
  # score   A set of match scores
  # bg      The random matching scores generated by the <ScorePWMDist> function

  bg <- sort(bg[!is.na(bg)]);
  qt <- (score-mean(bg)) / sd(bg);
  pv <- pt(qt, df=length(bg)-1, lower.tail=FALSE, log.p=TRUE);

  exp(pv)*2;
};

# Convert a set of DNA seed sequence to a position weight matrix
ConvertSeq2PWM <- function(seq) {
  # seq   A vector or DNAStringSet including seed sequences of the PWM

  if (is.character(seq) | is.factor(seq)) seq <- DNAStringSet(as.vector(seq))
  nch <- width(seq);
  seq <- seq[nch==max(nch)];
  pwm <- sapply(1:max(nch), function(i) table(substr(seq, i, i))[c('A', 'C', 'G', 'T')]);
  pwm[is.na(pwm)] <- 0;
  rownames(pwm) <- c('A', 'C', 'G', 'T');

  pwm;
}

# Plot the distribution of the relative position of PWM matches within sequences.
# If the matching p value is available, also plot the average -log10(p) within each interval
PlotPWMMatchPosition <- function(match.position, sequence.length, p.value=NA, motif.length=NA, min.count=5, max.interval=20) {
  # match.position  Numeric vector of the start positions of the matches.
  # sequence.length Numeric vector with the same length as <motif.position>. The length of the corresponding sequence.
  # p.value         The matching p value of each match
  # motif.length    The length of the motif. If available, motif.length-1 bases will be truncated from the end of the sequences
  # min.count       The minimal average number of matches per interval
  # max.interval    The maximal number of intervals

  # Adjust sequence length
  if (!identical(motif.length[1], NA)) sequence.length <- sequence.length - motif.length[1] + 1;

  # Number of intervals
  num.interval <- length(match.position)/max(5, min.count);
  num.interval <- min(num.interval, max.interval);
  num.interval <- max(5, 10*floor(num.interval/10));

  match.position <- pmin(match.position, sequence.length);

  pos <- ceiling(num.interval * match.position/sequence.length);
  cnt <- sapply(1:num.interval, function(i) length(pos[pos==i]));

  if (length(p.value) == length(match.position)) margin <- 5 else margin <- 2;
  par(mar=c(5,5,3,margin));
  barplot(cnt, width=1, space=0, xlab='Position in sequences', ylab='Number of matches', cex.lab=2,
          xlim=c(0, num.interval), xaxs='i', ylim=c(0, 1.1*max(cnt, na.rm=TRUE)));
  axis(1, at=c(0, num.interval/2, num.interval), label=c('Start', 'Center', 'End'));

  if (length(p.value) == length(pos)) {
    p <- -log10(p.value);
    p[p==Inf] <- max(p[p<Inf]);
    m <- sapply(1:num.interval, function(i) mean(p[pos==i]));
    par(new=TRUE);
    plot((1:num.interval)-0.5, m, axes=FALSE, type='l', col='lightgreen', xlab='', ylab='');
    points((1:num.interval)-0.5, m, col='darkblue', pch='*', cex=2);
    mtext("Mean of -Log10(P)", side=4, line=3, cex=1, col='darkblue');
    axis(4);
  }

  abline(v=num.interval/2, lty=2);
  box();

  invisible();
}
