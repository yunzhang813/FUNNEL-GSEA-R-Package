\name{FPCA.Fstats}
\alias{FPCA.Fstats}

\title{Functional F-statistics for time-course gene expression data}

\description{
This function takes time-course expression data and time points as input values. By using Functional Principal Component Analysis (FPCA) techniques, it returns the functional F-statistic for each gene.
}

\usage{
FPCA.Fstats(X, tt, rr=rep(1,length(tt)), selection_k="FVE", FVE_threshold=0.9)
}

\arguments{
  \item{X}{A gene expression matrix. Genes must be ordered as rows and
  each column represent a time point. Columns must be ordered by \code{tt}. }
  \item{tt}{A vector of ordered unique time points.}
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
}

\details{
  Technical details of FPCA procedure on time-course expression data is documented in 
  \cite{Wu, S. and Wu, H., 2013. BMC bioinformatics.}
}

\value{
  A vector of F-statistics. The larger the F-statistic is, the more evidence there is for the temporally differentially expressed gene (TDEG).
  
  Note: although the functional F-statistics are defined as ratios of variance between the constant expression model (null model) and the alternative model; due to the nature of functional data analysis, the functional F-statistics may not follow a standard F-distribution under the null.  To select significant TDEGs, we must resort to a resampling method.  Please consult the following paper for more details.
}

\references{
  \cite{Wu, S. and Wu, H., 2013. More powerful significant testing for time course gene expression data using functional principal component analysis approaches. BMC bioinformatics, 14(1), p.1.}
}

\author{
  Shuang Wu
}

\examples{

## Load the sample data
data("H3N2-Subj1")

## Standardize time and X as in FUNNEL.GSEA
tt2 <- (tt - min(tt))/diff(range(tt))
X2 <- t(scale(t(X)))

## Get functional F-statistics
Fstats <- FPCA.Fstats(X2, tt2)

## Sort the F-statistics from the largest to the smallest. 
## The top genes are more "significant".
sorted.genes <- names(sort(Fstats, decreasing = TRUE))
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ htest }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
