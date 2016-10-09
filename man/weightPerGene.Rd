\name{weightPerGene}
\alias{weightPerGene}

\title{Get weights (a.k.a. empirical gene set membership) of a list of
  genes of interest from a computed weight list.
}

\description{
  This is a convenience function that extracts weights (a.k.a. empirical
  gene set membership) of a vector of genes of interest from a weight
  list, which is typically computed by \code{FUNNEL.GSEA()}.
}

\usage{
weightPerGene(weight.list, genesOfInterest)
}

\arguments{
  \item{weight.list}{A list of weights (a.k.a. empirical gene set
  membership) computed from the elastic-net regression. }
  \item{genesOfInterest}{A vector of genes of interest. }
}

%% \details{
%%   more details than the description above
%% }

\value{
  A list of weights for the pre-specified genes.
}

\references{
  Yun Zhang, Juilee Thakar, Xing Qiu (2016) FUNNEL-GSEA: FUNctioNal ELastic-net
  Regression in time-course Gene Set Enrichment Analysis, submitted to
  Bioinformatics.
}

\author{
  Yun Zhang, Juilee Thakar, Xing Qiu
}

%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

%% \seealso{
%% }

\examples{

## Load the sample data
data("H3N2-Subj1")

## It takes about 10 minutes to run on my Laptop; YMMV.
\dontrun{t1 <- system.time(results1 <- FUNNEL.GSEA(X, tt, genesets=genesets))
    weightPerGene(results1$weight.list, c("TLR5", "C3AR1"))
}

}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ misc }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line

