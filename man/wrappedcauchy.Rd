\name{wrappedCauchy}
\title{Wrapped Cauchy Density Function}
\alias{dwrappedcauchy}
\alias{rwrappedcauchy}

\description{
Density, and random generation for the wrapped Cauchy circular distribution.   
}

\usage{
dwrappedcauchy(x, mu = circular(0), rho = exp(-1))
rwrappedcauchy(n, mu = circular(0), rho = exp(-1), control.circular=list())
}

\arguments{
  \item{x}{a vector. The object is coerced to class
    \code{\link{circular}}.}
  \item{n}{number of observations.}  
  \item{mu}{mean direction of the distribution as a \code{circular} object.}
  \item{rho}{concentration parameter of the distribution.  \code{rho}
    must be in the interval from 0 to 1.}
  \item{control.circular}{the attribute of the resulting object.}
}

\value{
\code{dwrappedcauchy} gives the density and \code{rwrappedcauchy} generates random deviates.
}

\author{Claudio Agostinelli and Ulric Lund}

\references{
Jammalamadaka, S. Rao and SenGupta, A. (2001). Topics in Circular Statistics, Section 2.2.7, World Scientific Press, Singapore.
}

\examples{

data1 <- rwrappedcauchy(100, mu=circular(0), rho=0.7, 
  control.circular=list(units="degrees"))
plot(data1)

ff <- function(x) dwrappedcauchy(x, mu=circular(pi), rho=0.7)
curve.circular(ff, join=TRUE, xlim=c(-2, 1), 
  main="Density of a Wrapped Cauchy Distribution \n mu=pi, rho=0.7")

}

\keyword{distribution}
