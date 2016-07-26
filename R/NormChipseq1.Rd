\name{NormChipseq1{Agr}}
\alias{Normalize ChIP-seq depth}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Normalization of ChIP-seq depth, Method 1
}
\description{
This function uses a 2-step procedure to normalize sequencing depth on a set of genomic locations betweeen ChIP-seq libraries. The first step subtracts background from the sequencing depth, and the second step normalizes data between samples. The inputs to this function require the processed ChIP-seq data to give the sequencing depth at a set of genomic regions, such as TSSs and enhancers. The regions are split into the same number of sub-regions and the average sequencing depth, at each sub-region, within each region, from each library, should be provided. 
}
\usage{     
 NormChipseq1(cov, col.start, col.end, log.transform = TRUE, 
              bg = 0, bg.start = c(), bg.end = c(), 
              trim.low = 0.05, trim.high = 0.95, exclude = c()) 
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
%%     ~~Describe \code{x} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (x) 
{
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
