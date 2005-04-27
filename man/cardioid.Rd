\name{Cardioid}
\title{Cardioid Density Function}
\alias{dcardioid}
\alias{rcardioid}

\description{
Density and random generation for the Cardioid circular distribution.  
}

\usage{
dcardioid(x, mu = 0, rho = 0)
rcardioid(n, mu = 0, rho = 0, units = c("radians", "degrees"), \dots)
}

\arguments{
  \item{x}{a vector. The object is coerced to class \code{\link{circular}}.}
\item{n}{number of observations.}  
\item{mu}{mean direction of the distribution (in the same 'units' of 'x').}
\item{rho}{concentration parameter of the distribution. Absolute value of \code{rho} must be less than 0.5.}
\item{units}{units of measure desired.}
\item{\dots}{parameters passed to \code{\link{circular}}.}
}

\value{
\code{dcardioid} gives the density and \code{rcardioid} generates random deviates.
}

\author{Claudio Agostinelli and Ulric Lund}

\references{
Jammalamadaka, S. Rao and SenGupta, A. (2001). Topics in Circular Statistics, Section 2.2.2, World Scientific Press, Singapore.
}

\keyword{distribution}
