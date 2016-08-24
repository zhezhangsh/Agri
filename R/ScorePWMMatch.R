# Match PWM to a sequence
ScorePWMMatch <- function(pwm, subject, both.strand=TRUE, min.score=85, min.match=1, local.max=TRUE) {
  # pwm          Same as the first parameter of the matchPWM() function: "For PWM: a rectangular character vector or rectangular DNAStringSet object ("rectangular" means that all elements have the same number of characters) with no IUPAC ambiguity letters, or a Position Frequency Matrix represented as an integer matrix with row names containing at least A, C, G and T (typically the result of a call to consensusMatrix)."
  # subject      Same as the first parameter of the matchPWM() function: "Typically a DNAString object. A Views object on a DNAString subject, a MaskedDNAString object, or a single character string, are also supported. IUPAC ambiguity letters in subject are ignored (i.e. assigned weight 0) with a warning."
  # both.strand  If TRUE, the reverse-complement sequence of <subject> will also be matched
  # min.score    Minimal match score, between 0 and 100
  # min.match    Minimal number of matches to return; reduce matching score if getting less matches using <min.score>
  # local.max    If TRUE, matches overlapped to any matches with higher scores will be removed.
  
  require("Biostrings");
  
  # Convert PWM from original format
  if (class(pwm) == "DNAStringSet" | is.character(pwm)) { # PWM is provided as DNAStringSet of seed sequences
    if (is.character(pwm)) pwm <- DNAStringSet(pwm); 
    nch <- nchar(pwm); 
    pwm <- pwm[nch==max(nch)]; 
    pwm <- sapply(1:max(nch), function(i) table(substr(pwm, i, i))[c('A', 'C', 'G', 'T')]); 
    pwm[is.na(pwm)] <- 0; 
  } if (is.matrix(pwm) | is.data.frame(pwm) | is.table(pwm)) { # PWM is provided as integer matrix
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
  
  # Max possible score of match
  mx <- maxScore(pwm);
  sc <- paste(min(100, max(0, min.score)), '%', sep='');

  # Reduce score cutoff if less matches than required
  count <- countPWM(pwm, subject, min.score = sc);
  if (min.match > count) {
    match <- matchPWM(pwm, subject, min.score = '0%', with.score = TRUE);
    score <- rev(sort(mcols(match)$score));
    sc    <- paste(floor(100*(score[min(min.match, length(score))]/mx)), '%', sep='');
  }
  
  ###################################################################
  # Match PWM
  match <- matchPWM(pwm, subject, min.score = sc, with.score = TRUE)
  mcols(match)$strand <- 1;
  
  # Match to reverse-complement sequence
  if (both.strand) { 
    match0 <- matchPWM(pwm, reverseComplement(subject), min.score = sc, with.score = TRUE); 
    range0 <- IRanges(length(subject)-end(match0)+1, length(subject)-start(match0)+1);
    match0@ranges <- range0;
    mcols(match0)$strand <- -1
    match <- c(match, match0);
    match <- match[order(start(match))];
  }
  
  # Remove match if it overlaps to another match with higher score
  if (local.max) {
    olap <- as.matrix(findOverlaps(match, match)); 
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