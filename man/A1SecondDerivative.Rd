\name{A1SecondDerivative}
\title{Second derivative of the Ratio of First and Zeroth Order Bessel Functions.}
\alias{A1SecondDerivative}

\description{Evaluates the second derivative of the second derivative of the Ratio of First and Zeroth Order Bessel Functions.}

\usage{
  A1SecondDerivative(kappa)
}

\arguments{
  \item{kappa}{non-negative numeric value at which to evaluate the second derivative of A1 function.}
}

\value{
  The value of the second derivative of A1 function in the point \code{kappa}.
}

\details{
  Formula (3.49) of  Fisher (1993), pag. 52 is implemented.
  The function uses \code{\link{A1}}, \code{\link{A1FirstDerivative}} and \code{\link{besselI}}.
}

\author{Claudio Agostinelli and Alessandro Gagliardi.}

\references{
  N.I. Fisher (1993) Statistical Analysis of Circular Data, Cambridge University Press.
}

\seealso{
  \code{\link{A1}}, \code{\link{A1FirstDerivative}}, \code{\link{besselI}}, \code{\link{A1inv}}.
}

\keyword{math}
