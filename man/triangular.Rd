\name{triangular}    
\title{Triangular Density Function}
\alias{dtriangular}
\alias{rtriangular}

\description{
Density and random generation for the Triangular circular distribution. 
}

\usage{
dtriangular(x, rho)
rtriangular(n, rho, units=c("radians", "degrees"), \dots)
}

\arguments{
\item{x}{a vector. The object is coerced to class \code{\link{circular}}.}
\item{n}{number of observations.}  
\item{rho}{concentration parameter of the distribution.  rho must be
  between 0 and \eqn{4/pi^2}.}
\item{units}{units of measure desired.}
\item{\dots}{parameters passed to \code{\link{circular}}.}
}

\value{
\code{dtriangular} gives the density and \code{rtriangular} generates
random deviates.
}

\author{Claudio Agostinelli and Ulric Lund}

\references{
Jammalamadaka, S. Rao and SenGupta, A. (2001). Topics in Circular Statistics, Section 2.2.3, World Scientific Press, Singapore.
}

\keyword{distribution}
