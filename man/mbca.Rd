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

\seealso{\code{\link[ade4]{dudi.pca}}, \code{\link[ade4]{bca}}, \code{\link{bcatest}}}

\examples{
### see ?coagulation in package SimComp
data("coagulation", package = "SimComp")

### level "S" is the standard, "H" and "B" are novel procedures
coagulation$Group <- relevel(coagulation$Group, ref = "S")

(bcaobj <- mbca(mat=as.matrix(coagulation[,c(2,3,4)]), coagulation$Group))
plot(bcaobj)

bcatest(bcaobj, Klist="Tukey", nperm=999)
}

\keyword{ misc }
