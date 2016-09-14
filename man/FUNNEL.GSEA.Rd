\name{FUNNEL.GSEA}
\alias{FUNNEL.GSEA}

\title{FUNNEL-GSEA: FUNctioNal ELastic-net Regression in time-course
  Gene Set Enrichment Analysis
}

\description{
  This is the main function of this package.  It takes gene expression
  matrix (\code{X}), time points (\code{tt}), a pre-defined list of gene sets
  (\code{genesets}), and several tuning parameters (see below) as input. It then
  apply functional data analysis techniques and elastic-net regression
  to define weights for genes that are shared by more than one gene sets
  (the overlapping genes). Finally, it applies an extended Mann-Whitney
  U test that incorporates both intergene correlation and weights
  computed from the elastic-net regression to pre-defined gene sets
  to obtain significance of non-constant temporal trend for each gene set.
}

\usage{
FUNNEL.GSEA(X, tt, genesets, lambda=10^-3.5, rr=rep(1,length(tt)),
            selection_k="FVE", FVE_threshold=0.9, nharm=3, centerfns=FALSE,
            equiv.threshold=0.01, lam1=0.4, lam2=0.01, alpha.level=0.05)
}

\arguments{
  \item{X}{A gene expression matrix. Genes must be ordered as rows and
  each column represent a time point. Row names must be gene names used
  in \code{genesets}. Columns must be ordered by \code{tt}. }
  \item{tt}{A vector of ordered unique time points.}
  \item{genesets}{A list of pre-defined gene sets.  Each element of this
  list is a vector of gene names that are rownames of the \code{X}
  matrix. Each element should have gene set name as element name.
  Note that due to data pre-processsing, many genes in genesets
  may not be present in \code{rownames(X)}. Internally, we define a new
  genesets object that is a collection of the intersection of each
  pre-defined pathway and \code{rownames(X)} for the actual computation. }
  \item{lambda}{Roughness penalty (smoothing) parameter used by function
  \code{smooth.basis} in package \code{fda}.  Empirical evidences show
  that the default value works well for a large collection of
  time-course gene expression data, as long as \code{X} is z-transformed and
  \code{tt} is normalized to [0,1].  In this function, we carry out
  these standardizations automatically.} 
  \item{rr}{Number of repetitions at each unique time points.} 
  \item{selection_k}{Method of choosing the number of principal
    components. Choices are: \code{"FVE"} (fraction of variance explained),
    which use scree plot approach to select number of principal
    components. See \code{FVE_threshold} below; a positive integer \code{K}:
    user-specified number of principal components. Default value:
    \code{"FVE"}.} 
  \item{FVE_threshold}{A positive number between 0 and 1. It is
    used with the option \code{selection_k = "FVE"} to select the number
    of principal components that explain at least \code{FVE_threshold}
    of total variation. Default value: 0.9. } 
  \item{nharm}{Number of eigenfunctions to be used to represent the
    overall temporal trend of each gene set. Default to 3.} 
  \item{centerfns}{An indicator if to center functions in \code{pca.fd()}. 
    Default: \code{FALSE}.}
  \item{equiv.threshold}{ This is a threshold to be used in function
    \code{equiv.regression}. Any eigenvalue that explains less than this
    proportion of total variance will be discarded.  Default to 0.01 (one
    percent of total variance). } 
  \item{lam1}{The L1 (LASSO) penalty used in elastic-net regression. } 
  \item{lam2}{The L2 (ridge) penalty used in elastic-ent regression. } 
  \item{alpha.level}{Significance level for gene set test (extended
  Mann-Whitney U test). Default to 0.05.} 
}

\details{
  Technical details of this procedure is documented in Yun Zhang,
  Juilee Thakar, Xing Qiu (2016) FUNNEL-GSEA: FUNctioNal ELastic-net
  Regression in time-course Gene Set Enrichment Analysis, submitted to
  Bioinformatics.
}

\value{
  A list with the following objects
  \item{pvals}{A vector of unadjusted p-values for pre-defined gene sets. }
  \item{weight.list}{A list of weights (a.k.a. empirical gene set
  membership) computed from the elastic-net regression. }
  \item{correlation}{Mean intergene correlation coefficient estimated
  from genes within each gene set.}
  \item{sig.genesets}{A vector of names (or indices) of significant gene sets.}
  \item{Fstats}{A vector of F-statistics, which are gene-level summary
  statistics used in extended Mann-Whitney U test.}
}

\references{
  Yun Zhang, Juilee Thakar, Xing Qiu (2016) FUNNEL-GSEA: FUNctioNal ELastic-net
  Regression in time-course Gene Set Enrichment Analysis, submitted to
  Bioinformatics.
}

\author{
  Yun Zhang, Juilee Thakar, Xing Qiu
}

\examples{

## Load the sample data
data("H3N2-Subj1")

## It takes about 10 minutes to run on my Laptop; YMMV.
\dontrun{t1 <- system.time(results1 <- FUNNEL.GSEA(X, tt, genesets=genesets))}

}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ htest }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
