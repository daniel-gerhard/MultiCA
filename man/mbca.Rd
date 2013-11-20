\name{mbca}
\alias{mbca}

\title{Between Group Analysis}

\description{Performs a principle component analysis with subject to an instrumental grouping structure.}

\usage{
mbca(mat, group)
}

\arguments{
  \item{mat}{An object of class matrix with response vectors in columns and observations per subject in rows.}
  \item{group}{A factor allocating the rows of the response matrix into several clusters.}
}


\details{The main part of the code is taken from functions dudi.pca and bca of package ade4.}

\author{Daniel Gerhard}

\seealso{\code{\link[ade4]{dudi.pca}}, \code{\link[ade4]{bca}}}

\keyword{ misc }
