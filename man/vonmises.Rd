\name{vonMises}
\title{von Mises Density Function}
\alias{dvonmises}
\alias{rvonmises}
\alias{pvonmises}

\description{
  Density, distribution function and random generation for the von Mises
  circular distribution.  
}

\usage{
dvonmises(x, mu, kappa)
pvonmises(q, mu, kappa, tol = 1e-20)
rvonmises(n, mu, kappa, units = c("radians", "degrees"), \dots)
}

\arguments{
  \item{x, q}{a vector. The object is coerced to class
    \code{\link{circular}}.}
\item{n}{number of observations.}  
\item{mu}{mean direction of the distribution.}
\item{kappa}{non-negative numeric value for the concentration parameter
  of the distribution.}
\item{units}{units of measure desired.}
\item{tol}{the precision in evaluating the distributon function.}
\item{\dots}{parameters passed to \code{\link{circular}}.}
}

\value{
\code{dvonmises} gives the density, \code{pvonmises} gives the
distribution function and \code{rvonmises} generates random deviates.
}

\author{Claudio Agostinelli and Ulric Lund}

\references{
Jammalamadaka, S. Rao and SenGupta, A. (2001). Topics in Circular Statistics, Section 2.2.4, World Scientific Press, Singapore.
}

\keyword{distribution}
