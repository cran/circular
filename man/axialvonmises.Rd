\name{axialvonMises}
\title{Axial von Mises Density Function}
\alias{daxialvonmises}

\description{
  Density for the axial von Mises circular distribution.  
}

\usage{
daxialvonmises(x, mu, kappa, l = 2)
}

\arguments{
  \item{x}{a vector. The object is coerced to class \code{\link{circular}}.}
 \item{mu}{mean direction of the distribution. The object is coerced to class \code{\link{circular}}}
 \item{kappa}{non-negative numeric value for the concentration parameter of the distribution.}
 \item{l}{a positive number. \code{l=2} provide the axial distribution in the range [0, pi].}
}

\value{
\code{daxialvonmises} gives the density.
}

\author{Claudio Agostinelli}

\references{
Jammalamadaka, S. Rao and SenGupta, A. (2001). Topics in Circular Statistics, Section 2.2.4, World Scientific Press, Singapore.
}

\keyword{distribution}

