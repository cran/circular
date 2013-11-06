\name{A1FirstDerivative}
\title{First derivative of the Ratio of First and Zeroth Order Bessel Functions.}
\alias{A1FirstDerivative}

\description{Evaluates the first derivative of the Ratio of First and Zeroth Order Bessel Functions}

\usage{
  A1FirstDerivative(kappa)
}

\arguments{
  \item{kappa}{non-negative numeric value at which to evaluate the first derivative of A1 function.}
}

\value{
  The value of the first derivative of A1 function in the point \code{kappa}.
}

\details{
  The formula (3.48) of  Fisher (1993), pag. 52 is implemented.
  The function uses \code{\link{A1}} and \code{\link{besselI}}.
}

\author{Claudio Agostinelli and Alessandro Gagliardi.}

\references{
  N.I. Fisher (1993) Statistical Analysis of Circular Data, Cambridge University Press.
}

\seealso{
  \code{\link{A1}}, \code{\link{besselI}}, \code{\link{A1inv}}.
}

\keyword{math}
