\name{wrappedNormal}
\title{Wrapped Normal Density Function}
\alias{dwrappednormal}
\alias{rwrappednormal}
\alias{pwrappednormal}
\alias{qwrappednormal}

\description{
Density, and random generation for the wrapped normal circular distribution.   
}

\usage{
rwrappednormal(n, mu = circular(0), rho = NULL, sd = 1, 
  control.circular = list())
dwrappednormal(x, mu = circular(0), rho = NULL, sd = 1, 
  K = NULL, min.k = 10)
pwrappednormal(q, mu = circular(0), rho = NULL, sd = 1, 
  from = NULL, K = NULL, min.k = 10, \dots)
qwrappednormal(p, mu = circular(0), rho = NULL, sd = 1, 
  from = NULL, K = NULL, min.k = 10, tol = .Machine$double.eps^(0.6), 
  control.circular = list(), \dots)
}

\arguments{
  \item{x, q}{vector of quantiles. The object is coerced to class
    \code{\link{circular}}.}
  \item{p}{vector of probabilities.}
  \item{n}{number of observations.} 
  \item{mu}{mean direction of the distribution as a \code{circular} object.}
  \item{rho}{concentration parameter of the distribution.  \code{rho}
    must be in the interval from 0 to 1.}
  \item{sd}{standard deviation of the (unwrapped) normal distribution.}
  \item{from}{if \code{NULL}  is set to \eqn{mu-pi}. This is the value from which the pwrappednormal and qwrappednormal are evaluated. It should be a \code{circular} object.}
  \item{K}{number of terms to be used in approximating the density.}
  \item{min.k}{minimum number of terms used in approximating the density.}
  \item{tol}{passed to \code{\link{uniroot}}.}
  \item{control.circular}{the attribute of the resulting object.}
  \item{\dots}{parameters passed to \code{\link{integrate}}.}
}

\value{
\code{dwrappednormal} gives the density and \code{rwrappednormal}
generates random deviates,  \code{pwrappednormal} gives the
distribution function, and \code{qwrappednormal} provides quantiles.
}

\author{Claudio Agostinelli and Ulric Lund}

\references{
Jammalamadaka, S. Rao and SenGupta, A. (2001). Topics in Circular Statistics, Section 2.2.7, World Scientific Press, Singapore.
}

\examples{

data1 <- rwrappednormal(100, mu=circular(0), rho=0.7, 
  control.circular=list(units="degrees"))
plot(data1)

ff <- function(x) dwrappednormal(x, mu=circular(pi), rho=0.7)
curve.circular(ff, join=TRUE, xlim=c(-1.5, 1), 
  main="Density of a Wrapped Normal Distribution \n mu=pi, rho=0.7")

ff <- function(x) pwrappednormal(x, mu=circular(pi), rho=0.7)
curve.circular(ff, join=FALSE, xlim=c(-2, 2), ylim=c(-2, 2), 
  to=(2*pi-3*.Machine$double.eps), modulo="asis", nosort=TRUE, 
  main="Probability of a Wrapped Normal Distribution \n mu=pi, 
  rho=0.7, from=0")

ff <- function(x) pwrappednormal(x, mu=circular(pi), rho=0.7, from=circular(pi))
curve.circular(ff, join=FALSE, xlim=c(-2, 2), ylim=c(-2, 2), from=-pi, 
  to=(pi-3*.Machine$double.eps), modulo="asis", nosort=TRUE, 
  main="Probability of a Wrapped Normal Distribution \n mu=pi, 
  rho=0.7, from=pi")

plot(qwrappednormal, from=0, to=1)
plot(function(x) qwrappednormal(p=x, mu=circular(pi)), from=0, to=1)
}


\keyword{distribution}
