\name{Cardioid}
\title{Cardioid Density Function}
\alias{dcardioid}
\alias{rcardioid}

\description{
Density and random generation for the Cardioid circular distribution.  
}

\usage{
dcardioid(x, mu = circular(0), rho = 0)
rcardioid(n, mu = circular(0), rho = 0, control.circular=list())
}

\arguments{
  \item{x}{a vector. The object is coerced to class \code{\link{circular}}.}
 \item{n}{number of observations.}  
 \item{mu}{mean direction of the distribution. The object is coerced to class \code{\link{circular}}.}
 \item{rho}{concentration parameter of the distribution. Absolute value of \code{rho} must be less than 0.5.}
 \item{control.circular}{the coordinate system used in the output of \code{rcardioid}. See \code{\link{circular}} for details.}
}

\value{
\code{dcardioid} gives the density and \code{rcardioid} generates random deviates.
}

\author{Claudio Agostinelli and Ulric Lund}

\references{
Jammalamadaka, S. Rao and SenGupta, A. (2001). Topics in Circular Statistics, Section 2.2.2, World Scientific Press, Singapore.
}

\examples{
  set.seed(1234) 
  resrad <- rcardioid(n=10)
  set.seed(1234)
  resdeg <- rcardioid(n=10, control.circular=list(units="radians", zero=pi))  
  max(abs(resrad - conversion.circular(resdeg, zero=0)))
}

\keyword{distribution}

