\name{contrList}
\alias{contrList}

\title{Construct a list of pseudo contrast matrices.}

\description{Constructs a list of contrasts to perform predefined multiple comparisons between groups in the bcatest function.}

\usage{
contrList(n, type=c("Dunnett", "Tukey", "AVE"), base=1)
}

\arguments{
  \item{n}{A (possibly named) vector of sample sizes for each group.}
  \item{type}{Type of contrast.}
  \item{base}{an integer specifying which group is considered the baseline group for Dunnett contrasts.}
}


\details{Computes the requested list matrices for comparisons of multivariate mean levels.}

\author{Daniel Gerhard}

\seealso{\code{\link{bcatest}}, \code{\link[multcomp]{contrMat}}}

\keyword{ misc }
