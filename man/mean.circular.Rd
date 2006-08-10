\name{mean.circular}
\title{Mean Direction}
\alias{mean.circular}

\description{
Returns the mean direction of a vector of circular data. 
}
 
\usage{
\method{mean}{circular}(x, na.rm=FALSE, control.circular=list(), \dots)
}

\arguments{
  \item{x}{a vector. The object is coerced to class
    \code{\link{circular}}.}
  \item{na.rm}{logical, indicating if \code{\link{NA}}'s should be
	omitted.}
  \item{control.circular}{the attribute of the resulting object.}
  \item{\dots}{further arguments passed to or from other methods.}
}

\value{
Returns the mean direction of the data as an object of class \code{circular} with the attribute given by \code{control.circular} or from \code{x} if missed in \code{control.circular}.
}

\details{
Each observation is treated as a unit vector, or point  on the unit
circle.  The resultant vector of the observations is found, and the
direction of the resultant vector is returned. An \code{\link{NA}} is
returned if the resultant length (see \code{\link{rho.circular}}) is
less than \code{\link{.Machine}} 
}

\author{Claudio Agostinelli and Ulric Lund}

\references{
Jammalamadaka, S. Rao and SenGupta, A. (2001). Topics in Circular Statistics, Section 1.3, World Scientific Press, Singapore.
}

\seealso{
\code{\link{var.circular}}, \code{\link{summary.circular}},
\code{\link{mle.vonmises}}, \code{\link{rho.circular}} and \code{\link{.Machine}}.
}

\examples{
# Compute the mean direction of a random sample of observations.
data <- circular(runif(50, circular(0), pi))
mean(data)
}

\keyword{univar}
