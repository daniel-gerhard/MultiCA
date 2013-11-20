\name{bcatest}
\alias{bcatest}

\title{Multiple hypothesis testing in between component analysis.}

\description{Permutation testing of cluster distances in between component analysis using a maximum statistic.}

\usage{
bcatest(x, Klist, nf=NULL, nperm=999)
}

\arguments{
  \item{x}{An object of class mbca.}
  \item{Klist}{A list of matrices defining the hypotheses of interest.}
  \item{nf}{Number of factors used for comparisons. If NULL (default) the full rank is used.}
  \item{nperm}{Number of permutations. default=999}
}


\details{Adaption of the function randtest of package ade4 to perform a maximum test of multiple partial hypotheses.}

\author{Daniel Gerhard}

\seealso{\code{\link[ade4]{randtest}}}

\keyword{ misc }
