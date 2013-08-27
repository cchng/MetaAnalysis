\name{Fisher}
\alias{Fisher}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Fisher's combined probability test
}
\description{
A meta-analysis approach using combined significance.
}
\usage{
Fisher(merged, k.threshold = 3)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{merged}{
    A data frame with raw p-values of genes (rows) in each study (columns).
}
  \item{k.threshold}{
    Minimum number of data sets included in the meta-analysis. Default = 3.
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
%%  Carolyn Ch'ng
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
function (merged, k.threshold = 3) 
{
    F.stat <- function(dat) {
        pvals <- as.numeric(dat[2:length(dat)])
        k <- length(which(!is.na(pvals)))
        F <- -2 * sum(log(pvals), na.rm = TRUE)
        return(c(dat["Group.1"], F = F, k = k))
    }
    meta <- data.frame(t(apply(merged, 1, F.stat)))
    meta$F <- as.numeric(levels(meta$F))[meta$F]
    meta$k <- as.numeric(levels(meta$k))[meta$k]
    meta$df <- 2 * meta$k
    meta <- meta[which(meta$k >= k.threshold & meta$Group.1 != 
        ""), ]
    meta$P <- pchisq(meta$F, df = meta$df, lower.tail = FALSE)
    meta$Q <- qvalue(meta$P)$qvalues
    meta$BH <- p.adjust(meta$P, method = "fdr")
    meta.ord <- meta[order(meta$BH), ]
    return(meta.ord)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
