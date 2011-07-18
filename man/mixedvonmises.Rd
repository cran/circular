\name{mixedvonmises}    
\title{Mixture of von Mises Distributions}
\alias{dmixedvonmises}
\alias{rmixedvonmises}
\alias{pmixedvonmises}

\description{
  Density and random generation for the mixed von Mises circular distribution.  
}

\usage{
dmixedvonmises(x, mu1, mu2, kappa1, kappa2, prop)
rmixedvonmises(n, mu1, mu2, kappa1, kappa2, prop, control.circular = list())
pmixedvonmises(q, mu1, mu2, kappa1, kappa2, prop, from=NULL, tol = 1e-020)
}

\arguments{
  \item{x, q}{a vector. The object is coerced to class \code{\link{circular}}.}
  \item{n}{number of observations.}
  \item{mu1}{mean direction of one of the two von Mises distributions as a \code{circular} object.}
  \item{mu2}{mean direction of the other von Mises distribution as a \code{circular} object.}
  \item{kappa1}{concentration parameter of one of the two von Mises distributions.}
  \item{kappa2}{concentration parameter of the other von Mises distribution.}
  \item{prop}{mixing proportion.}
  \item{from}{if \code{NULL} is set equal to \eqn{0} (Notice the
  difference from the correspoding vonmises distribution). This is the value from which the pmixedvonmises is evaluated. It should be a \code{circular} object.}
  \item{tol}{the precision in evaluating the distributon function or the
  quantile.}
  \item{control.circular}{the attribute of the resulting object.}
}

\value{
\code{dmixedvonmises} gives the density, \code{pmixedvonmises} gives the
distribution function and \code{rmixedvonmises} generates random deviates. 
}

\author{Claudio Agostinelli and Ulric Lund}

\seealso{
\code{\link{dvonmises}} \code{\link{pvonmises}} and \code{\link{rvonmises}}
}

\examples{

x <- rmixedvonmises(n=100, mu1=circular(0), mu2=circular(pi), kappa1=15, 
  kappa2=15, prop=0.5)
plot(x)


}


\keyword{distribution}
