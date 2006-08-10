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
pvonmises(q, mu, kappa, tol = 1e-020)
rvonmises(n, mu, kappa, control.circular=list())
}

\arguments{
  \item{x, q}{a vector. The object is coerced to class
    \code{\link{circular}}.}
\item{n}{number of observations.}  
\item{mu}{mean direction of the distribution. The object is coerced to class \code{\link{circular}}.}
\item{kappa}{non-negative numeric value for the concentration parameter
  of the distribution.}
\item{control.circular}{the attribute of the resulting object.}
\item{tol}{the precision in evaluating the distributon function.}
}

\value{
\code{dvonmises} gives the density, \code{pvonmises} gives the
distribution function and \code{rvonmises} generates random deviates.

Since version 0.3-5 the random deviates are generated using a C code.
}

\author{Claudio Agostinelli, Ulric Lund and Harry Southworth}

\references{
Jammalamadaka, S. Rao and SenGupta, A. (2001). Topics in Circular Statistics, Section 2.2.4, World Scientific Press, Singapore.
}

\examples{

data1 <- rvonmises(100, circular(0), 10, control.circular=list(units="degrees"))
plot(data1)

ff <- function(x) dvonmises(x, mu=circular(pi), kappa=10)
curve.circular(ff, join=TRUE, xlim=c(-2.3, 1), main="Density of a VonMises Distribution \n mu=pi, kappa=10")

ff <- function(x) pvonmises(x, mu=circular(pi), kappa=10)
curve.circular(ff, join=FALSE, xlim=c(-2, 2), ylim=c(-2, 1), to=(2*pi-3*.Machine$double.eps), modulo="asis", nosort=TRUE, main="Probability of a VonMises Distribution \n mu=pi, kappa=10")

}


\keyword{distribution}
