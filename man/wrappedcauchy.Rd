\name{wrappedCauchy}
\title{Wrapped Cauchy Density Function}
\alias{dwrappedcauchy}
\alias{rwrappedcauchy}

\description{
Density, and random generation for the wrapped Cauchy circular distribution.   
}

\usage{
dwrappedcauchy(x, mu = 0, rho = exp(-1))
rwrappedcauchy(n, mu = 0, rho = exp(-1), units = c("radians", "degrees"), \dots)
}

\arguments{
  \item{x}{a vector. The object is coerced to class
    \code{\link{circular}}.}
  \item{n}{number of observations.}  
  \item{mu}{mean direction of the distribution.}
  \item{rho}{concentration parameter of the distribution.  \code{rho}
    must be in the interval from 0 to 1.}
  \item{units}{units of measure desired.}
  \item{\dots}{parameters passed to \code{\link{circular}}.}

}

\value{
\code{dwrappedcauchy} gives the density and \code{rwrappedcauchy} generates random deviates.
}

\author{Claudio Agostinelli and Ulric Lund}

\references{
Jammalamadaka, S. Rao and SenGupta, A. (2001). Topics in Circular Statistics, Section 2.2.7, World Scientific Press, Singapore.
}

\keyword{distribution}
