\name{wMWUTest}
\alias{wMWUTest}

\title{
    An extened Mann-Whitney U Rank Sum Test That Incorporates Pre-computed
    Weights and Correlation.
}

\description{
    This is an extension of the two-sample Mann-Whitney U test which
    incorporates pre-calculated weights for the correlated inputs in the
    test group. Note that the pre-calculated correlation only applies to the
    test group (the gene set of interest).  The correlation of the
    background genes is assumed to be zero.  Pre-calculated weights are
    typically computed by function \code{FUNNEL.GSEA}.
}

\usage{
    wMWUTest(test.index, statistics, weight=NULL, correlation=0, df=Inf)
}

\arguments{
    \item{test.index}{A vector of indices (or names) of genes that belong
        to the test group (the gene set to be tested).}
    \item{statistics}{A (named) numeric vector contains all the elements
    (summary statistics or observed values) from *both* groups, such
    that statitics[test.index] is a vector of statistics of the test set.}
    \item{weight}{A numeric vector, of the same length as
      \code{test.index}, giving weights for the elements in the test
      group. If \code{NULL}, all elements are weighted equally with unit
      weight.}
    \item{correlation}{An estimate of the correlation in the test
      group. Genes in the second group are assumed to be independent
      of each other and of the genes in the test group.}
    \item{df}{Degrees of freedom, based on which the correlation is
    estimated. For \code{funnel.gsea}, we define \code{df==N.timepoints-1}. }
}

%\details{
    %%  ~~ If necessary, more details than the description above ~~
                                                               %}

\value{
    P-values (\code{less}, \code{greater}) for one-sided left- and
    right-tail tests, respectively.
}

\references{
  \cite{Barry, W.T., Nobel, A.B., and Wright, F.A. (2008).
    A statistical framework for testing functional categories in microarray data. \emph{Annals of Applied Statistics}, 286-315.}
  
  \cite{Wu, D, and Smyth, GK (2012). Camera: a competitive gene set test
    accounting for inter-gene correlation. \emph{Nucleic Acids Research},
    40(17), e133-e133.}
  
  \cite{Yun Zhang, Juilee Thakar, Xing Qiu (2016) FUNNEL-GSEA:
    FUNctioNal ELastic-net Regression in time-course Gene Set Enrichment
    Analysis.  \emph{submitted to Nucleic Acids Research}.
  }
}

\author{
  Yun Zhang, Juilee Thakar, Xing Qiu
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

## We define the first 20 "genes" to be a gene set
test.index <- 1:20

## Add some true signal (>) to the first 10 test genes
stat[1:8] <- stat[1:8]+1 

pL <- wilcox.test(stat[test.index], stat[-test.index], alternative = "less")$p.value
pU <- wilcox.test(stat[test.index], stat[-test.index], alternative = "greater")$p.value

wMWUTest(test.index, stat)
## compare it with the following. pU is what we are looking for
c("less"=pL, "greater"=pU)

## With just 0.1 correlation, p-values are not significant anymore
wMWUTest(test.index, stat, correlation=0.1)

## Our results are equivalent to the implementation provided by limma
library(limma)
rankSumTestWithCorrelation(test.index, stat, correlation=0.1)

## First set of weight: attenuates the signal. With weights<1 for all
## signal-carrying genes, the test is less significant
ww1 <- runif(length(test.index), 0, 1)
wMWUTest(test.index, stat, weight=ww1, correlation=0.1)

## Second set of weight: All the signal-carrying genes have weight==1;
## the rest 15 genes have less weights.  Now the p-value is
## signficant again!
ww2 <- c(rep(1, 8), runif(12, 0, 1))
wMWUTest(test.index, stat, weight=ww2, correlation=0.1)

## In the context of FUNNEL.GSEA
## Load the sample data
data("H3N2-Subj1")

## It takes about 10 minutes to run on my Laptop; YMMV.
\dontrun{t1 <- system.time(results1 <- FUNNEL.GSEA(X, tt, geneset=geneset))

    geneset2 <- lapply(geneset, function(z) { intersect(z, rownames(X)) })
    gg1 <- geneset2[[1]]
    ww1 <- rr1[["weight.list"]][["GLYCOLYSIS_GLUCONEOGENESIS"]]
    rho <- rr1$correlation

    ## The test
    o1 <- wMWUTest(gg1, rr1$Fstats, ww1, rho, df=15)
    o1["greater"]                       #p-value for gene test

    ## should be the same as the p-value below
    rr1$pvals["GLYCOLYSIS_GLUCONEOGENESIS"]

}
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ htest }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
