\name{wrappedNormal}
\title{Wrapped Normal Density Function}
\alias{dwrappednormal}
\alias{rwrappednormal}

\description{
Density, and random generation for the wrapped normal circular distribution.   
}

\usage{
dwrappednormal(x, mu = 0, rho, sd = 1, K, min.k = 10)
rwrappednormal(n, mu = 0, rho, sd = 1, units = c("radians", "degrees"), \dots)
}

\arguments{
  \item{x}{a vector. The object is coerced to class
    \code{\link{circular}}.}
  \item{n}{number of observations.}  

  \item{mu}{mean direction of the distribution.}
  \item{rho}{concentration parameter of the distribution.  \code{rho}
    must be in the interval from 0 to 1.}
  \item{sd}{standard deviation of the (unwrapped) normal distribution.}
  \item{K}{number of terms to be used in approximating the density.}
  \item{min.k}{minimum number of terms used in approximating the density.}
  \item{units}{units of measure desired.}
  \item{\dots}{parameters passed to \code{\link{circular}}.}
}

\value{
\code{dwrappednormal} gives the density and \code{rwrappednormal} generates random deviates.
}

\author{Claudio Agostinelli and Ulric Lund}

\references{
Jammalamadaka, S. Rao and SenGupta, A. (2001). Topics in Circular Statistics, Section 2.2.7, World Scientific Press, Singapore.
}

\keyword{distribution}
