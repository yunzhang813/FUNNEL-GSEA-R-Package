\name{plotWeight}
\alias{plotWeight}

\title{Plot the weights (a.k.a. empirical gene set membership) obtained
  from \code{FUNNEL.GSEA}.
}

\description{
  This function takes a list of weights (a.k.a. empirical membership),
  an index (or name) of a specific gene set, then it produces a heatmap
  (by calling function \code{pheatmap}) of these weights.  
}

\usage{
plotWeight(weight.list, geneset.index, main="Weight", ...)
}

\arguments{
  \item{weight.list}{A list of weights (a.k.a. empirical membership) of
  genes in each gene set.  This object is typically obtained from
  \code{FUNNEL.GSEA}. }
  \item{geneset.index}{A specific index or name of the geneset to be plotted. }
  \item{main}{The title of the figure. See \code{pheatmap}.}
  \item{...}{Other graphical parameters to be passed to function \code{pheatmap()}.}
}

\references{
  \cite{Zhang, Y., Topham, D.J., Thakar, J. and Qiu, X. (2017). FUNNEL-GSEA:
    FUNctioNal ELastic-net Regression in time-course Gene Set Enrichment
    Analysis. \emph{Bioinformatics}, btx104.}
}

\author{
  Yun Zhang, Juilee Thakar, Xing Qiu
}

\seealso{\code{\link{pheatmap}}
}

\examples{

## Load the sample data
data("H3N2-Subj1")

## It takes about 10 minutes to run on my Laptop; YMMV.
\dontrun{t1 <- system.time(results1 <- FUNNEL.GSEA(X, tt, genesets=genesets))
## Generate a heatmap of weights for a particular gene set
    plotWeight(results1$weight.list, "PRIMARY_IMMUNODEFICIENCY")
}
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{hplot}
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line

