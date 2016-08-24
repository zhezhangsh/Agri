\name{NormChipseq1}
\alias{Normalize ChIP-seq depth}
\title{
Normalization of ChIP-seq depth, Method 1
}
\description{
This function uses a 2-step procedure to normalize sequencing depth on a set of genomic locations across ChIP-seq libraries. The first step subtracts background from the sequencing depth, and the second step normalizes data between samples. The inputs to this function require the processed ChIP-seq data to give the sequencing depth at a set of genomic regions, such as TSSs and enhancers. Each region is split into the same number of sub-regions and the average sequencing depth of each sub-region within each region from each library should be provided. 
}
\usage{     
NormChipseq1 <- function(cov, log.transform=TRUE, bg=0, 
                         bg.start=c(), bg.end=c(), trim.low=0.05, 
                         trim.high=0.95, max.seed=25000, exclude=c())
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{cov}{A list of numeric matrix with matching rows and columns. Each matrix corresponds to a ChIP-seq library to be normalized. Each row of the matrixes is a genomic region with the same number of sub-regions as columns. Each cell is the average sequencing depth at each sub-region.}
  \item{log.transform}{If TRUE, the normalization will be performed after log2-transformation of the matrixes, as log2(depth+1).}
  \item{bg}{If the a single number, the universal background level to be subtracted from all matrixes. If a numeric vector with the same length as the length of <cov> matrxies, the background levels to be subtracted from corresponding matrixes. Ignored otherwise.}
  \item{bg.start, bg.end}{The start/end column indexes of sub-regions to be used to calculate background level. Each parameter could be empty, one single number, or a vector of number, but both must have the same length.}
  \item{trim.low, trim.high}{To calculate the background level, the top and bottom regions/rows to be trimmed.}
  \item{exclude}{The indexes of the regions/rows to be excluded from calculating background.}
  \item{max.seed}{The maximum number of rows to be used as seed of normalization, to save computational time. Rows will be selected randomly if there are more rows/regions in the matrixes than this number. Use all rows/regions if empty.}
}
\details{
This is a two-step procedure, respectively addressing the difference of total throughput (number of sequencing reads) and ChIP efficiency between libraries.

The frist step is to adjust the matrixes by their background sequencing depth that reflects the overall throughput of the libraries. The background depth can be provided by the <bg> parameter, as a univeral number for all libraries or a vector of the same length as the length of the <cov> list of matrixes. It can also be calculated using a subset of columns in the matrixes. The global mean of these columns of each matrix will be used as background level of each library. To calculate the mean, the rows/regions can be trimed to remove certain percents of top/bottom values and filtered by excluding a subset. 

The second step is to perform a loess normalization between the samples on background-adjusted data. The normalization will be performed column by column on all matrixes. Random rows/columns can be selected as seeds of normalization to save computational time.

The <cov> parameter should include a list of non-negative matrixes with exactly the same rows. Rows of all matrixes correspond to the same genomic regions. The columns of these matrixes could correspond to sub-regions of each row with the same or different length. 
}
\value{
The normalized <cov> matrixes as a named list. 
}
\references{

}
\author{
Zhe Zhang
}
\note{
}
\seealso{
}
\examples{

## Create random matrix with log-normal distribution
# Create 10 pairs of random numbers
random.number <- cbind(sample(seq(-0.2, 0.2, 0.01), 10), 
                       sample(seq(0.75, 1.25, 0.01), 10)); 
# Create 10 random vectors based on the random numbers
random.vector <- apply(random.number, 1, function(num) 
                        num[2]*rlnorm(5000, num[1]));
# Create column adjustment factor
column.adjust <- exp(c(seq(0, 2, 0.1), seq(2, 0, -0.1)[-1]))
# Create 10 random matrixes
random.matrix <- lapply(1:10, function(i) {v <- random.vector[,i];
                         sapply(column.adjust, function(a) 
                          v*a+rlnorm(5000, -5)); }); 

names(random.matrix) <- paste('S', 1:10, sep='');

## Normalize by known background
norm.matrix1 <- NormChipseq1(random.matrix, log.transform = TRUE, 
                             bg = random.number[, 1]);
                             
## Normalize by reference columns
norm.matrix2 <- NormChipseq1(random.matrix, log.transform = TRUE, 
                             bg.start = c(1, 36), bg.end = c(5, 41));
                             
## Effect of normalization
x0 <- lapply(random.matrix, function(x) x[, 21]); 
x1 <- lapply(norm.matrix1,  function(x) x[, 21]); 
x2 <- lapply(norm.matrix2,  function(x) x[, 21]); 

boxplot(x0, outline=FALSE); # No normalization
boxplot(x1, outline=FALSE); # Normalize by known background
boxplot(x2, outline=FALSE); # Normlaize by reference columns as background
}
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
