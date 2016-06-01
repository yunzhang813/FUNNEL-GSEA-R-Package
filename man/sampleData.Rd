\name{sampleData}
\alias{sampleData}
\docType{data}
\title{
Simulated and Preprocessed Data
}
\description{
The sample data are simulated and preprocessed as required. The users may clean their raw data as they prefer, and make the clean data in the format this package requires by functions such as \code{scaleX, scaleTime, modifyGeneset}.
}

\usage{data("sampleData")}

\format{
The data set contains the following:
\describe{
\item{\code{X}}{Gene expression matrix. Genes in rows and time points in columns. Row names must be in gene names, with the same annotation used in the \code{geneset} later on. Columns must be ordered corresponding to \code{time}. The whole matrix is z-transformed gene-wise, i.e. row-wise.}
\item{\code{time}}{Vector of ordered time points. They should be scaled to the range of [0,1].}
\item{\code{geneset}}{List of gene sets. Each gene set is an element of the list. Names of the gene sets should be recorded in the names of each corresponding element. Genes in \code{geneset} should be annotated in the same way rownames of the \code{X} matrix.}
\item{\code{fdexpr}}{Smoothed gene expression curves obtained through the function \code{smoothExpr(X, time)}.}
}
}

\examples{
data(sampleData)
}
\keyword{datasets}


