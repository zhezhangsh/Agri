# A more efficient way of PWM to DNA sequence match of 1 or more PWMs and 1 or more sequences
ScorePWMMatchSet <- function(pwm.set, subject.set, num.cluster=4, both.strand=TRUE, local.max=TRUE, min.score=85, min.match=0) {
  # pwm          A list of PWM objects to be passed to the ScorePWMMatch() function: "For PWM: a rectangular character vector or rectangular DNAStringSet object ("rectangular" means that all elements have the same number of characters) with no IUPAC ambiguity letters, or a Position Frequency Matrix represented as an integer matrix with row names containing at least A, C, G and T (typically the result of a call to consensusMatrix)."
  # subject      A DNAStringSet or character vector of DNA sequences
  # both.strand  If TRUE, the reverse-complement sequence of <subject> will also be matched
  # local.max    If TRUE, matches overlapped to any matches with higher scores will be removed.
  # min.score    Minimal match score, between 0 and 100
  # min.match    Minimal number of matches within all sequences to return; reduce matching score if getting less matches using <min.score>
  # num.cluster  Number of clusters for parallele computing
  
  require("Biostrings");
  require("parallel");
  
  # Union sequences
  if (is.character(subject.set)) subject.set <- DNAStringSet(subject.set); 
  seq.all <- unlist(subject.set); 
  
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
      ind <- olap@queryHits[olap@subjectHits==i];
      if (length(ind) == 0) Views(subject.set[[i]]) else {
        x <- m[ind];
        Views(subject.set[[i]], start = start(x) - stt[i] + 1, end = end(x) - stt[i] + 1); 
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
    if (is.character(pwm) | is.factor(pwm)) pwm <- DNAStringSet(as.vector(pwm)); 
    nch <- width(pwm); 
    pwm <- pwm[nch==max(nch)]; 
    pwm <- sapply(1:max(nch), function(i) table(substr(pwm, i, i))[c('A', 'C', 'G', 'T')]); 
    pwm[is.na(pwm)] <- 0; 
  } else if (is.matrix(pwm) | is.data.frame(pwm) | is.table(pwm)) { # PWM is provided as integer matrix
    pwm <- as.matrix(pwm)[1:4, , drop=FALSE];
    rownames(pwm) <- c('A', 'C', 'G', 'T');
    mx  <- max(pwm, na.rm=TRUE); 
    if (mx <= 1) { # PWM is given as proportion, guess the total number of seeds and transform the matrix
      mn  <- min(pwm[pwm>0 & !is.na(pwm)]);  # non-zero minimum
      adj <- 1/mn; 
      pwm <- round(pwm*adj)
    }
  } else {
    stop('Unknown format of PWM: ', class(pwm), '\n');
  }
  pwm[is.na(pwm)] <- 0;
  if (!both.strand) pwm <- list(pwm) else {
    rev <- pwm[4:1, ncol(pwm):1];
    rownames(rev) <- rownames(pwm);
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
  str   <- rep(c(1, -1), sapply(match, length)); 
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
    olap <- cbind(olap@queryHits, olap@subjectHits); 
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