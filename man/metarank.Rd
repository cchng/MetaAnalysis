\name{metarank}
\alias{metarank}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Meta-rank Analysis
}
\description{
  A meta-analysis approach using meta-ranks.
}
\usage{
metarank(merged)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{merged}{
    A data frame with raw p-values of genes (rows) in each study (columns).
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
function (merged) 
{
    k <- apply(merged[, c(2:ncol(merged))], 1, function(x) {
        length(which(!is.na(x)))
    })
    merged <- merged[which(k >= 3), ]
    within.rank <- apply(merged[, c(2:ncol(merged))], 2, function(x) {
        rank(x, na.last = "keep")
    })
    merged$Ave.Rank <- rowMeans(within.rank, na.rm = TRUE)
    sqindex <- (within.rank - merged$Ave.Rank)^2
    merged$Heterogeneity <- rowSums(sqindex, na.rm = TRUE)
    merged$Meta.Rank <- rank(merged$Ave.Rank)
    merged.ord <- merged[order(merged$Meta.Rank), ]
    return(merged.ord)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
