\name{weightedRankSumTestWithCorrelation}
\alias{weightedRankSumTestWithCorrelation}

\title{
Weighted Rank Sum Test with Correlation
}

\description{
This is an extension of the two-sample Mann-Whitney U test (a.k.a. Wilcoxon rank sum test), which incorporates different weighting for the correlated elements in the test group.
}

\usage{
weightedRankSumTestWithCorrelation(index, statistics, weight=NULL, correlation=0, df=Inf)
}

\arguments{
  \item{index}{vector of indices of the elements belonging to the test group.}
  \item{statistics}{numeric vector contains all the elements (summary statistics or observed values) from both groups.}
  \item{weight}{numeric vector, of the same length as \code{index}, giving weights for the elements in the test group. If \code{NULL}, all elements are weighted equally with unit weight.}
  \item{correlation}{correlation in the test group. Elements in the second group are assumed independent of each other and of the test group.} 
  \item{df}{degrees of freedom, based on which the correlation is estimated.} 
}

%\details{
%%  ~~ If necessary, more details than the description above ~~
%}

\value{
P-values (\code{less}, \code{greater}) for one-sided left and right tails respectively.
}

\references{
Barry, W.T., Nobel, A.B., and Wright, F.A. (2008).
A statistical framework for testing functional categories in microarray data.
\emph{Annals of Applied Statistics}, 286-315.

Wu, D, and Smyth, GK (2012). Camera: a competitive gene set test accounting for inter-gene correlation.
\emph{Nucleic Acids Research}, 40(17), e133-e133.
}

\author{
Yun Zhang
}

%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{wilcox.test}} performs the standard Wilcoxon rank sum test.

\code{\link[limma]{rankSumTestWithCorrelation}} in the \code{limma} package performs the correlation extention of the rank sum test.
}

\examples{
set.seed(1)
stat <- rnorm(100)
index <- 1:10
stat[index] <- stat[1:10]+1
weight <- runif(length(index), 0, 1)

weightedRankSumTestWithCorrelation(index, stat)
weightedRankSumTestWithCorrelation(index, stat, correlation=0.1)
weightedRankSumTestWithCorrelation(index, stat, weight=weight, correlation=0.1)

## Wilcoxon rank sum test
group <- rep(2,100)
group[index] <- 1
group <- factor(group)
wilcox.test(stat ~ group, alternative = "less")
wilcox.test(stat ~ group, alternative = "greater")

## rank sum test with correlation
library(limma)
rankSumTestWithCorrelation(index, stat, correlation=0.1)
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ weightedRankSumTestWithCorrelation }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
