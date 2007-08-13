\name{[.circular}
\title{Extract or Replace Parts of a Circular Object}
\alias{[.circular}

\description{
  Operators act on vectors and matrices to extract or
  replace subsets, methods for Circular Data.
}

\usage{
\method{[}{circular}(x, i, \dots)
}

\arguments{
  \item{x}{object from which to extract elements.}
  \item{i, \dots}{elements to extract or replace.}
}

\author{Claudio Agostinelli}

\examples{
x <- circular(matrix(rwrappednormal(n=100, mu=circular(0)), nrow=5))
dim(x)
x[1,]
x[,1]
x[,1, drop=FALSE]
}

\keyword{array}
