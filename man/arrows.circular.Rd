\name{arrows.circular}
\alias{arrows.circular}
\title{Add Arrows to a Circular Plot}
\description{
Draw arrows in a circular plot. 
}
\usage{
arrows.circular(x, y = NULL, x0 = 0, y0 = 0, na.rm = FALSE, 
  shrink = 1, plot.info = NULL, zero = NULL, rotation = NULL, ...)
}
\arguments{
  \item{x}{a vector. The object is coerced to class \code{\link{circular}}.}
  \item{y}{a vector with the same length as \code{x}.}
  \item{x0}{a vector of origins (x axis).}
  \item{y0}{a vector of origins (y axis).}
  \item{na.rm}{logical, indicating if \code{\link{NA}}'s should be omitted.}
  \item{shrink}{parameter that controls the size of the plotted circle.
  Default is 1.  Larger values shrink the circle, while smaller values
  enlarge the circle.}
  \item{plot.info}{an object from \code{\link{plot.circular}} that
  contains information on the \code{zero}, the \code{rotation} and \code{next.points}.}
  \item{zero}{the zero used in the plot. Ignored if \code{plot.info} is provided.}
  \item{rotation}{the rotation used in the plot. Ignored if \code{plot.info} is provided.}
  \item{\dots}{futher parameters passed to \code{\link{arrows}}.}
}

\author{Claudio Agostinelli}

\note{The function call \code{\link{arrows}} and it is not a method of \code{\link{arrows}}.}

\seealso{\code{\link{arrows}}}

\examples{
  plot(rvonmises(10, circular(0), kappa=1))
  arrows.circular(rvonmises(10, circular(0), kappa=1))
  arrows.circular(rvonmises(10, circular(0), kappa=1), y=runif(10), col=2)
  arrows.circular(rvonmises(10, circular(0), kappa=1), y=runif(10), 
    x0=runif(10, -1, 1), y0=runif(10, -1, 1), col=3)
}

\keyword{aplot}
\keyword{hplot}
