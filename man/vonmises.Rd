\name{vonMises}
\title{von Mises Density Function}
\alias{dvonmises}
\alias{rvonmises}
\alias{pvonmises}
\alias{qvonmises}

\description{
  Density, distribution function, random generation and quantiles for
  the von Mises circular distribution. 
}

\usage{
rvonmises(n, mu, kappa, control.circular=list())
dvonmises(x, mu, kappa)
pvonmises(q, mu, kappa, from=NULL, tol = 1e-020)
qvonmises(p, mu = circular(0), kappa=NULL, from=NULL, tol = .Machine$double.eps^(0.6), 
  control.circular = list(), \dots)
}

\arguments{
  \item{x, q, p}{a vector. The \code{x} and \code{q} objects are coerced to class
    \code{\link{circular}}.}
  \item{n}{number of observations.}  
  \item{mu}{mean direction of the distribution. The object is coerced to class \code{\link{circular}}.}
  \item{kappa}{non-negative numeric value for the concentration parameter
    of the distribution.}
  \item{from}{if \code{NULL}  is set to \eqn{mu-pi}. This is the value from which the pvonmises and qvonmises are evaluated. It should be a \code{circular} object.}
  \item{tol}{the precision in evaluating the distributon function or the
  quantile.}
  \item{control.circular}{the attribute of the resulting object.}
  \item{\dots}{parameters passed to \code{\link{integrate}}.}
}

\value{
\code{dvonmises} gives the density, \code{pvonmises} gives the
distribution function, \code{rvonmises} generates random deviates and
\code{qvonmises} provides quantiles.

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
curve.circular(ff, join=TRUE, xlim=c(-2.3, 1),
  main="Density of a VonMises Distribution \n mu=pi, kappa=10")

ff <- function(x) pvonmises(x, mu=circular(pi), kappa=10)
curve.circular(ff, join=FALSE, xlim=c(-2, 2), ylim=c(-2, 1), 
  to=(2*pi-3*.Machine$double.eps), modulo="asis", nosort=TRUE, 
  main="Probability of a VonMises Distribution \n mu=pi, kappa=10")

plot(function(x) qvonmises(x, mu=circular(0), kappa=5), from=0, to=1)
##curve do not work!
plot(function(x) qvonmises(x, mu=circular(pi), kappa=5), from=0, to=1)
plot(function(x) qvonmises(x, mu=circular(pi), kappa=5, from=circular(pi/2)), from=0, to=1)
}


\keyword{distribution}
