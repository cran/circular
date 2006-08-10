\name{[.circular}
\title{Extract or Replace Parts of a Circular Object}
\alias{[.circular}

\description{
  Operators act on vectors and matrices to extract or
  replace subsets, methods for Circular Data.
}

\usage{
x[i]
x[i, j, \dots]
x[i, j, \dots , drop=TRUE]
}

\arguments{
  \item{x}{object from which to extract elements.}
  \item{i, j, \dots}{elements to extract or replace.}
  \item{drop}{For matrices.  If \code{TRUE} the
    result is coerced to the lowest possible dimension (see examples
    below).}
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
