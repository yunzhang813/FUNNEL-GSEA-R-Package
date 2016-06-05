\name{H3N2-Subj1}
\alias{H3N2-Subj1}
\docType{data}
\title{
  Time-course gene expression data collected from a H3N2 influenza
  challenge study (Subject 1).
}

\description{
  This is the processed gene expression data extracted from Huang et al
  2011 PloS Comp. Genet. study, with GEO accession number
  GSE52428. Expressions are log2-transformed and filtered by the
  criterion IQR>=0.3.
}

\usage{data("H3N2-Subj1")}

\format{
  The data set contains the following:
  \describe{
    \item{\code{X}}{A 11189x16 pre-processed gene expression matrix taken
      from a H3N2 human influenza challenge study (Subject1). Genes are
      ordered as rows and each column represent a time point. Row names must
      be gene names used in \code{genesets}. Columns must be ordered
      by \code{tt}.}
    \item{\code{tt}}{A vector of ordered time points.}
    \item{\code{genesets}}{List of 186 MSigDB:KEGG gene sets.
      Each gene set is an element of the list. Genes in \code{genesets}
      should be annotated in the same way rownames of the \code{X} matrix does.}
  }
}

\examples{
data(H3N2-Subj1)
}
\keyword{datasets}
