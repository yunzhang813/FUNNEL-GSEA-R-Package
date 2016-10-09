\name{equiv.regression}
\alias{equiv.regression}

\title{An equivalence between functional and multivariate regression}

\description{
  This function takes functional covariates (\code{xfd}) and response (\code{yfd}) as
  input, and returns a list of equivalent multivariate covariates (\code{Xmat})
  and response (\code{y}).  Running penalized (or ordinary least square)
  regression on \code{y} and \code{Xmat} is equivalent to the corresponding concurrent
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
  Technical details of this transformation is documented in \cite{Yun Zhang,
  Juilee Thakar, Xing Qiu (2016)}.
}

\value{
  A list with the following objects
  \item{y}{A matrix of dimension KxM, where K is the number of
    nontrivial eigenvalues kept in the transformation and M is the
    dimension of \code{yfd} (number of response curves). }
  \item{Xmat}{A matrix of dimension KxP, where K is the number of
    nontrivial eigenvalues kept in the transformation and P is the
    dimension of \code{xfd} (number of functional covariates).}
  \item{cut.value}{ \code{threshold} x total variance. We only keep
    those eigenvalues that explain greater or equal to this amount of
    total variance in the transformation. }
}

\references{
  \cite{Yun Zhang, Juilee Thakar, Xing Qiu (2016) FUNNEL-GSEA: FUNctioNal ELastic-net
  Regression in Time-course Gene Set Enrichment Analysis, submitted to
  Bioinformatics.}
}

\author{
  Yun Zhang, Juilee Thakar, Xing Qiu
}

\examples{
## library(FUNNEL)
library(fda)
library(quadrupen)

## Let us first create some functional objects
data("H3N2-Subj1")

## Remove genes in predefined gene sets that are not present in X the
## filtered input data
genenames <- rownames(X)
newGenesets <- lapply(genesets, function(z) { intersect(z, genenames) } )

## Standardize time and X so that the smoothing penality is applicable 
## to any real data with comparable signal-to-noise ratio
tt2 <- (tt - min(tt))/diff(range(tt))
X2 <- t(scale(t(X)))

## Smoothing
mybasis <- create.bspline.basis(range(tt2), length(tt2)+4-2, 4, tt2)
mypar <- fdPar(mybasis, 2, lambda=10^-3.5)
fdexpr <- smooth.basis(tt2, t(X2), mypar)$fd

## Take genes C5AR1 and C3AR1 as two examples
gene.i <- c("C5AR1", "C3AR1")

## They belong to the following two pathways
newGeneset.i <- newGenesets[c("NEUROACTIVE_LIGAND_RECEPTOR_INTERACTION", "COMPLEMENT_AND_COAGULATION_CASCADES")]

## The response is just the smoothed curves of these two genes
yfd <- fdexpr[gene.i]

##############
## Method 1 ##
##############

## Let us use the first 3 eigen-functions of both pathways as covariates
xfd <- FUNNEL:::PCA.genesets(fdexpr, newGeneset.i, nharm = 3, centerfns = FALSE)$harmonics

## Calculate the equivalent multivariate regression datasets
equiv <- equiv.regression(yfd, xfd, threshold = 0.01)
equiv.Y <- equiv$y; colnames(equiv.Y) <- gene.i     #3x2 matrix
equiv.X <- equiv$Xmat                               #3x6 matrix
colnames(equiv.X) <- paste(rep(names(newGeneset.i), each=3), rep(paste0("eigfun", 1:3), length(newGeneset.i)), sep=".")

## Now we can run multivariate elastinet regression on the equivalent X and Y, as
## implemented in package quadrupen
en <- elastic.net(equiv.X, equiv.Y[, "C3AR1"], lambda1 = 0.4, lambda2 = 0.01, intercept = FALSE, normalize = FALSE)

## beta.en are the regression coefficients
beta.en <- as.numeric(attributes(en)$coef)
names(beta.en) <- colnames(equiv.X)

##############
## Method 2 ##
##############

## Alternatively, let us try using the first 2 eigen-functions and increase the variance 
## threshold in the equivalence rotated regression to 0.1.
xfd2 <- FUNNEL:::PCA.genesets(fdexpr, newGeneset.i, nharm = 2, centerfns = FALSE)$harmonics

## Calculate the equivalent multivariate regression datasets
equiv2 <- equiv.regression(yfd, xfd2, threshold = 0.1)
equiv2.Y <- equiv2$y; colnames(equiv2.Y) <- gene.i     #3x2 matrix
equiv2.X <- equiv2$Xmat                                #3x6 matrix
colnames(equiv2.X) <- paste(rep(names(newGeneset.i), each=2), rep(paste0("eigfun", 1:2), length(newGeneset.i)), sep=".")

## Now we can run multivariate elastinet regression on the equivalent X and Y, as
## implemented in package quadrupen
en2 <- elastic.net(equiv2.X, equiv2.Y[, "C3AR1"], lambda1 = 0.4, lambda2 = 0.01, intercept = FALSE, normalize = FALSE)

## beta.en are the regression coefficients
beta.en2 <- as.numeric(attributes(en2)$coef)
names(beta.en2) <- colnames(equiv2.X)
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%% \keyword{models}
