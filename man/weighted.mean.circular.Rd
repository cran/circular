\name{weighted.mean.circular}
\title{Weighted Mean Direction}
\alias{weighted.mean.circular}

\description{
Returns the weighetd mean direction of a vector of circular data. 
}
 
\usage{
\method{weighted.mean}{circular}(x, w, na.rm=FALSE,
  control.circular=list(), \dots)
}

\arguments{
  \item{x}{a vector. The object is coerced to class
    \code{\link{circular}}.}
  \item{w}{a numerical vector of weights the same length as \code{x}
    giving the weights to use for elements of \code{x}.}
  \item{na.rm}{logical, indicating if \code{\link{NA}}'s should be
	omitted.}
  \item{control.circular}{the attribute of the resulting object.}
  \item{\dots}{further arguments passed to or from other methods.}
}

\value{
Returns the weighted mean direction of the data as an object of class \code{circular} with the attribute given by \code{control.circular} or from \code{x} if missed in \code{control.circular}.
}

\details{
Each observation is treated as a unit vector, or point  on the unit
circle.  The resultant vector of the observations is found, and the
direction of the resultant vector is returned. An \code{\link{NA}} is
returned if the weighted resultant length is less than
\code{\link{.Machine}}.

If \code{w} is missing then all elements of \code{x} are given the same
weight, otherwise the weights coerced to numeric by \code{as.numeric} and normalized to sum to one.

Missing values in \code{w} are not handled specially and so give a
missing value as the result. However, zero weights are handled
specially and the corresponding \code{x} values are omitted from the
computation.
}

\author{Claudio Agostinelli}

\seealso{
\code{\link{mean.circular}}
}

\examples{
# Compute the weighted mean direction of a random sample of observations.
x <- circular(runif(50, circular(0), pi))
w <- runif(50, 0, 1)
weighted.mean(x, w)
}

\keyword{univar}
