\name{equiv.regression}
\alias{equiv.regression}

\title{An equivalence between functional and multivariate regression.}

\description{
  This function takes functional covariates (xfd) and response (yfd) as
  input, and returns a list of equivalent multivariate covariates (Xmat)
  and response (y).  Running penalized (or ordinary least square)
  regression on y and Xmat is equivalent to the corresponding concurrent
  functional regression with constant beta.
}

\usage{
equiv.regression(yfd, xfd, threshold=0.01)
}

\arguments{
  \item{yfd}{A list of functional objects as response variables.}
  \item{xfd}{A list of functional objects as covariates.}
  \item{threshold}{Any eigenvalue that explains less than this
    proportion of total variance will be discarded.  Default to 0.01 (one
    percent of total variance). }
}

\details{
  Technical details of this transformation is documented in Yun Zhang,
  Juilee Thakar, Xing Qiu (2016) FUNNEL-GSEA: FUNctioNal ELastic-net
  Regression in time-course Gene Set Enrichment Analysis, submitted to
  Nucleic Acids Research.
}

\value{
  A list with the following objects
  \describe{
    \item{y}{A matrix of dimension KxM, where K is the number of
      nontrivial eigenvalues kept in the transformation and M is the
      dimension of yfd (number of response curves). }
    \item{Xmat}{A matrix of dimension KxP, where K is the number of
      nontrivial eigenvalues kept in the transformation and P is the
      dimension of xfd (number of functional covariates).}
    \item{cut.value}{ \code{threshold} x total variance. We only keep
      those eigenvalues that explain greater or equal to this amount of
      total variance in the transformation. }
  }
}

\references{
  Yun Zhang, Juilee Thakar, Xing Qiu (2016) FUNNEL-GSEA: FUNctioNal ELastic-net
  Regression in time-course Gene Set Enrichment Analysis, submitted to
  Nucleic Acids Research.
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
## library(FUNNEL)
library(fda)
library(quadrupen)

## Let us first create some functional objects
data("H3N2-Subj1")

## Remove genes in predefined genesetss that are not present in X the
## filtered input data
genenames <- rownames(X)
newGenesets <- lapply(genesets, function(z) { intersect(z, genenames) } )

## Standardize timepoint and X so that the optimum roughness/L1/L2
## penality parameters are applicable
tt2 <- round((tt - min(tt))/diff(range(tt)), 2)
X2 <- t(scale(t(X)))

## smoothing
mybasis <- create.bspline.basis(range(tt2), length(tt2)+4-2, 4, tt2)
mypar <- fdPar(mybasis, 2, lambda=10^-3.5)
fdexpr <- smooth.basis(tt2, t(X2), mypar)$fd

## take genes C5AR1 and C3AR1 as two examples
gene.i <- c("C5AR1", "C3AR1")


## They belong to the following two pathways
newGeneset.i <- newGenesets[c("NEUROACTIVE_LIGAND_RECEPTOR_INTERACTION", "COMPLEMENT_AND_COAGULATION_CASCADES")]

## The response is just the smoothed curves of these two genes C3AR1
yfd <- fdexpr[gene.i]

## Let us use the first 3 eigen-functions of both pathways as covariates
xfd <- FUNNEL:::PCA.genesets(fdexpr, newGeneset.i, nharm = 3, centerfns = FALSE)$harmonics

## Calculate the equivalent multivariate regression datasets
equiv <- equiv.regression(yfd, xfd, threshold = 0.01)
Y <- equiv$y; colnames(Y) <- gene.i     #3x7 matrix
X <- equiv$X                            #3x6 matrix
colnames(X) <- paste(rep(names(newGeneset.i), each=3), rep(paste0("eigfun", 1:3), length(newGeneset.i)), sep=".")

## Now we can run multivariate elastinet regression on X and Y, as
## implemented in package quadrupen.
en <- elastic.net(X, Y[, "C3AR1"], lambda1 = 0.4, lambda2 = 0.1,
                             intercept = FALSE, normalize = FALSE)

## beta.en are the regression coefficientsx
beta.en <- as.numeric(attributes(en)$coef)
names(beta.en) <- colnames(X)

}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%% \keyword{models}
